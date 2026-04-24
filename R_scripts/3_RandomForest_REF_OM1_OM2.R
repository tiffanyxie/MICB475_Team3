# Adapted from "Random Forest" by Avril Metcalfe-Roach
# R Script for multi-class RF model to predict REF, OM1, OM2

#### Part 1: Load Libraries and Data ####
library(randomForest)
library(caret)
library(ranger)
library(pROC)
library(boot)
library(ggplot2)
library(phyloseq)
library(tidyverse)



#Load dataset and aggregate reads to Genus level
load("phylo_soil_genus.RData")

# Load ASVs of interest

#DeSeq results
deseq_om1_ref<-read.delim("output/DeSEQ_OM1_vs_REF.csv",sep=",")
deseq_om2_ref<-read.delim("output/DeSEQ_OM2_vs_REF.csv",sep=",")
deseq_om2_om1<-read.delim("output/DeSEQ_OM2_vs_OM1.csv",sep=",")


deseq_results<-rbind(deseq_om1_ref,deseq_om2_ref,deseq_om2_om1) %>%
  filter(!grepl("Incertae_Sedis",Genus))
unique_asvs<-deseq_results %>% pull(ASV) %>% unique()
unique_genus<-deseq_results %>% pull(Genus) %>% unique()




#### Part 2: Format Data ####

# Normalize microbiome data and select only the significant taxa



# CLR transformation
ps_clr = phylo_soil_genus %>% microbiome::transform('clr') 

#Filter taxa
ps_filt = prune_taxa(unique_asvs,ps_clr)

# Melt the dataset
df = psmelt(ps_filt) 
# Add Z transformation to each Genus individually
# (mean of zero, standard deviation of 1)
df_transformed = df %>% 
  group_by(Genus) %>% 
  mutate(Abundance = as.numeric(scale(Abundance))) %>% #Used GenAI to troubleshoot and added as.numeric
  ungroup()

# Final table should ONLY contain the outcome and explanatory variables, each as their own column. For now we'll also include the sample id.
df_pivot = df_transformed %>% 
  select(Sample,pH,Total.Carbon,Total.Nitrogen, pH, CN.Ratio,LTSP.Treatment,Genus,Abundance) %>% 
  # Turn each Genus into its own column
  pivot_wider(names_from = Genus, values_from = Abundance) 

# Remove rows with NA values in the metadata
df_noNA = df_pivot %>% na.omit()

# Remove the sample ID column - otherwise the code will try to use it as an explanatory variable (just like the microbial genera).
# We do this after pivoting (try doing it before pivoting, see what happens!)
df_final = df_noNA %>% select(-Sample)


#### Part 3: Set up RF ####


# Predictors and outcome
predictors = df_final %>% select(-LTSP.Treatment) %>% as.data.frame()


outcome = df_final %>% pull(LTSP.Treatment) %>%
  factor(levels=c("REF","OM1","OM2"))

#Set up k-fold cross-validation
# Randomly subsets the rows into k equal bins.
k = 5
set.seed(421)
folds = createFolds(outcome, k = k, list = TRUE)

# Each of these folds will take a turn being the test dataset.
str(folds)


#Hyperparameters

# mtry: number of variables that will be used per forest. 
#       High = overfitting, low = uninformative

# splitrule: affects how decision trees are calculated.
#            Use gini or extratrees for boolean outcomes (ex. subject)
#            Use variance for continuous outcomes (ex. age)

# min.node.size: Related to tree complexity. Larger = simpler tree.
# Often best as a proportional fraction of your sample size.

# These are generic values. Depending on your dataset, you may need to adjust the numeric ranges up or down.
tune_grid = expand.grid(mtry = c(3,6,10), 
                        splitrule = c("gini","extratrees"),
                        min.node.size = c(2,3,4))

#### Part 4: Run RF ####
source('3_randomforest_functions_modified.R')

soil_model = run_rf(X = predictors, y = outcome, 
                  fold_list = folds,
                  hyper = tune_grid, 
                  rngseed = 421)

names(soil_model)

#Save model
save(soil_model,file = "output/soil_model_OM1_OM2_REF.Rdata")

#### Part 5: Interpret Results ####

#Test results with confusion matrix
test_true<-soil_model$test_labels %>% pull(true_labels) %>%
  factor(levels=c("REF","OM1","OM2"))

#For confusion matrix need factor with eventual result
# Use max probability to determine result
test_result<-data.frame(Result = soil_model$test_labels %>% 
                          select(REF,OM1,OM2) %>% max.col()) %>%
  mutate(Result = case_when(Result == 1 ~ "REF",
                            Result == 2 ~ "OM1",
                            Result == 3 ~ "OM2")) %>%
  pull(Result) %>% 
  factor(levels=c("REF","OM1","OM2"))

test_conf_mat<-confusionMatrix(test_result,test_true)
test_conf_mat$table

#Training results with confusion matrix
train_true<-soil_model$train_labels %>% pull(true_labels) %>%
  factor(levels=c("REF","OM1","OM2"))

#For confusion matrix need factor with eventual result
# Use max probability to determine result
train_result<-data.frame(Result = soil_model$test_labels %>% 
                          select(REF,OM1,OM2) %>% max.col()) %>%
  mutate(Result = case_when(Result == 1 ~ "REF",
                            Result == 2 ~ "OM1",
                            Result == 3 ~ "OM2")) %>%
  pull(Result) %>% 
  factor(levels=c("REF","OM1","OM2"))

train_conf_mat<-confusionMatrix(train_result,train_true)
train_conf_mat$table


#Plot importance values
soil_model$importance %>% 
  # Data are automatically arranged by decreasing importance - turn it into a factor.
  # Otherwise the features will show up alphabetically in the plot.
  mutate(Feature = factor(.$Feature,levels = .$Feature)) %>% 
  ggplot(aes(Feature,MeanDecreaseGini,fill=MeanDecreaseGini)) +
  geom_col() +
  theme_classic(base_size=18) +
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1)) +
  ylab('Importance (Gini)') + xlab(NULL)
ggsave("output/importance_plot_OM1_OM2_REF.png",width = 14, height=4, units=c("in"))

