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



#Load dataset and aggregate reads to Genus level
load("phylo_soil.RData")
ps <-  phylo_soil %>%
  tax_glom('Genus')

#### Part 2: Format Data ####

# Normalize microbiome data and select only the significant taxa
# CLR transformation
ps_clr = ps %>% microbiome::transform('clr') 

#Filter taxa
test<- tax_table(ps_clr) %>% rownames() %>% unique() %>% head(n=20)## REPLACE WITH ACTUAL TAXA
ps_filt = prune_taxa(test,ps_clr)

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
  pivot_wider(names_from = Genus, values_from = Abundance) %>%
  mutate(LTSP.Treatment = if_else(grepl("OM",LTSP.Treatment),"OM","REF"))

# Remove rows with NA values in the metadata
df_noNA = df_pivot %>% na.omit()

# Remove the sample ID column - otherwise the code will try to use it as an explanatory variable (just like the microbial genera).
# We do this after pivoting (try doing it before pivoting, see what happens!)
df_final = df_noNA %>% select(-Sample)


#### Part 3: Set up RF ####


# Predictors and outcome
predictors = df_final %>% select(-LTSP.Treatment)

outcome = df_final %>% pull(LTSP.Treatment) %>%
  factor(levels=c("REF","OM"))
  #factor(levels=c("REF","OM1","OM2"))

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

pd_model = run_rf(X = predictors, y = outcome, 
                  fold_list = folds,
                  hyper = tune_grid, 
                  rngseed = 421)

names(pd_model)

#### Part 5: Interpret Results ####

# Generate ROC curve

roc_test = roc(pd_model$test_labels$true_labels,
               pd_model$test_labels$predicted_probabilities)
roc_train = roc(pd_model$train_labels$true_labels,
                pd_model$train_labels$predicted_probabilities)

# True positive rate = sensitivity
# False positive rate = 1-specificity
ggplot() +
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
                           auc(roc_train), pd_model$auc_train_ci[1], pd_model$auc_train_ci[2],
                           auc(roc_test), pd_model$auc_test_ci[1], pd_model$auc_test_ci[2]), 
           size = 6) +
  theme_minimal(base_size=18)

#Compile results into table
roc_data = data.frame(Dataset = 'RF Tutorial Data',
                      Training_AUC = round(pd_model$auc_train, 2),
                      Training_AUC_CI = paste0(round(pd_model$auc_train_ci[1], 2), "-", 
                                               round(pd_model$auc_train_ci[2], 2)),
                      Testing_AUC = round(pd_model$auc_test, 2),
                      Testing_AUC_CI = paste0(round(pd_model$auc_test_ci[1], 2), "-", 
                                              round(pd_model$auc_test_ci[2], 2)))

#Plot importance values
pd_model$importance %>% 
  # Data are automatically arranged by decreasing importance - turn it into a factor.
  # Otherwise the features will show up alphabetically in the plot.
  mutate(Feature = factor(.$Feature,levels = .$Feature)) %>% 
  ggplot(aes(Feature,MeanDecreaseGini,fill=MeanDecreaseGini)) +
  geom_col() +
  theme_classic(base_size=18) +
  theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1)) +
  ylab('Importance (Gini)') + xlab(NULL)


