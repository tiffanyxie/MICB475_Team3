# Adapted from "Random Forest" by Avril Metcalfe-Roach


#### Part 1: Load Libraries and Data ####
library(randomForest)
library(caret)
library(ranger)
library(pROC)
library(boot)
library(ggplot2)
library(phyloseq)
library(tidyverse)
library(patchwork)



#Load dataset and aggregate reads to Genus level

load("phylo_soil_genus.RData")

# Load ASVs of interest

#DeSeq results

deseq_om2_om1<-read.delim("output/DeSEQ_OM2_vs_OM1.csv",sep=",")


deseq_results<- filter(deseq_om2_om1,
                       !grepl("Incertae_Sedis",Genus))
unique_asvs<-deseq_results %>% pull(ASV) %>% unique()
unique_genus<-deseq_results %>% pull(Genus) %>% unique()

#### Part 2: Format Data ####

# Normalize microbiome data and select only the significant taxa

# CLR transformation
ps_binary <- subset_samples(phylo_soil_genus, LTSP.Treatment %in% c("OM1", "OM2")) ##Gemini used
ps_clr = ps_binary %>% microbiome::transform('clr') 

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
#predictors = df_final %>% select(-LTSP.Treatment,-CN.Ratio,-Total.Carbon) %>% as.data.frame()


outcome = df_final %>% pull(LTSP.Treatment) %>%
  factor(levels=c("OM1","OM2"))

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
source('3_randomforest_functions.R')

soil_model = run_rf(X = predictors, y = outcome, 
                    fold_list = folds,
                    hyper = tune_grid, 
                    rngseed = 421)
names(soil_model)

save(soil_model,file = "output/soil_model_OM1_OM2.Rdata")



####Figures ####
load("output/soil_model_OM1_OM2.Rdata")

roc_test =  roc(soil_model$test_labels$true_labels,
               soil_model$test_labels$predicted_probabilities)
roc_train = roc(soil_model$train_labels$true_labels,
                soil_model$train_labels$predicted_probabilities)

roc_plot<-ggplot() +
  # Training data: this is a type of control
  geom_line(aes(x = 1 - roc_train$specificities, 
                y = roc_train$sensitivities), 
            color = "red",size=1) +
  # Test data: tells us the strength of the prediction
  geom_line(aes(x = 1 - roc_test$specificities,
                y = roc_test$sensitivities), 
            color = "black",size=1) +
  geom_abline(slope = 1, intercept = 0, color = "gray", linetype = "dashed",size=1) +
  labs(x = "False Positive Rate", y = "True Positive Rate") +
  annotate("text", x = 0.7, y = 0.2, 
           label = sprintf("Train (red): %.2f (%.2f-%.2f)\nTest (black): %.2f (%.2f-%.2f)",
                           auc(roc_train), soil_model$auc_train_ci[1], soil_model$auc_train_ci[2],
                           auc(roc_test), soil_model$auc_test_ci[1], soil_model$auc_test_ci[2]), 
           size = 12/.pt) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size=10))
roc_plot

roc_data = data.frame(Dataset = 'RF Tutorial Data',
                      Training_AUC = round(soil_model$auc_train, 2),
                      Training_AUC_CI = paste0(round(soil_model$auc_train_ci[1], 2), "-", 
                                               round(soil_model$auc_train_ci[2], 2)),
                      Testing_AUC = round(soil_model$auc_test, 2),
                      Testing_AUC_CI = paste0(round(soil_model$auc_test_ci[1], 2), "-", 
                                              round(soil_model$auc_test_ci[2], 2)))

soil_model$importance

importance_plot_labels<-soil_model$importance


importance_plot<-soil_model$importance %>% 
  mutate(Feature = gsub("g__","",Feature)) %>%
# Data are automatically arranged by decreasing importance - turn it into a factor.
# Otherwise the features will show up alphabetically in the plot.
  mutate(Feature = factor(.$Feature,levels = .$Feature)) %>% 
  ggplot(aes(Feature,MeanDecreaseGini)) + #,fill=MeanDecreaseGini)) +
  geom_col() +
  theme_classic() +
  scale_y_continuous(expand=c(0,0)) +
  theme(axis.text.x = element_text(angle=50, vjust = 1, hjust=1),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12)) +
  ylab('Importance (Gini)') + xlab(NULL) +
  labs(fill = "Mean Decrease Gini")
importance_plot


#### Save plots ####
roc_plot
ggsave("figures/om1_om2_roc.png",
       units=c("in"),
       width = 4,
       height = 3)
importance_plot
ggsave("figures/om1_om2_importance.png",
       units=c("in"),
       width = 6,
       height = 3)
ggsave("figures/om1_om2_importance.svg",
       units=c("in"),
       width = 6,
       height = 3)

