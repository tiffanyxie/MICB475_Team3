# Run random forest on each fold

# Average all the folds
average_rf = function(all_labels_train,all_labels_test,
                      feature_importance_values){
  
  
  # Combine the test data labels and predictions (each fold is currently in a separate df)
  # Each sample is only included in a test dataset once, so no need to average.
  test_labels = bind_rows(all_labels_test)
  
  # Combine the train data labels and predictions (each fold is currently in a separate df)
  # Each sample is included 9 times, so need to average.
  train_labels = bind_rows(all_labels_train) %>% 
    group_by(row,true_labels) %>% 
    summarize(REF = mean(REF),
              OM1 = mean(OM1),
              OM2 = mean(OM2)) %>% 
    ungroup()
  
  # Calculate the average importance values for each variable.
  # We will use Reduce to add together all values that are at the same row/column coordinate across all datasets in feature_importance_values, then divide by the number of datasets to get the mean.
  mean_feature_importance = Reduce("+", feature_importance_values) / length(feature_importance_values)
  
  # Convert to a data frame and add a Feature column
  importance_df = data.frame(Feature = rownames(mean_feature_importance), mean_feature_importance)
  # Remove row names
  rownames(importance_df) = NULL
  # Sort the results from most to least important. We will use MeanDecreaseGini.
  importance_df = importance_df %>% arrange(-MeanDecreaseGini)
  
  results = list(test_labels = test_labels,
                 train_labels = train_labels,
                 importance = importance_df)
  return(results)
}


run_rf = function(X, y, fold_list,
                  hyper, rngseed = 421,
                  kfold=T) {
  
  # X = predictors; y = outcome
  # fold_list = folds; hyper = tune_grid;
  # rngseed=421; kfold=T
  
  # Calculate number of total folds
  number_of_folds = length(fold_list)
  
  # Create series of empty vectors. We'll store model outputs here.
  # AUC scores of the ROC curves
  train_auc_scores = c()
  test_auc_scores = c()
  
  # These will contain data frames of the following:
  # For every sample, we will record a) its actual outcome (PD or Control) and 
  # b) the outcome predicted by the random forest model.
  # This will be used to calculate how accurate our model is.
  all_labels_train = list()
  all_labels_test = list()
  
  # Importance values
  feature_importance_values = list()
  
  for (fold in fold_list) {
    if (kfold == F){
      fold = fold_list[[1]]
    }
    
    # Create train and test datasets.
    X_train_fold = X[-fold, ]
    y_train_fold = y[-fold]
    X_test_fold = X[fold, ]
    y_test_fold = y[fold]
    
    # This will tell the RF command how to perform the RF.
    train_control = trainControl(method = "cv", # K-fold cross validation
                                 number = number_of_folds) # 10 folds

    
    # Use hyperparameter tuning to optimize each parameter.
    # Note that optimal settings are chosen based on ROC/AUC - prone to overfitting!
    set.seed(rngseed) # Reproducible randomness
    rf_model = suppressWarnings(caret::train(X_train_fold, y_train_fold, # training dataset
                                             method = "ranger",
                                             trControl = train_control, # Perform tuning
                                             tuneGrid = hyper,
                                             metric = "Accuracy" #Changed from ROC to train
    ))
    
    # Finally, run random forest using the optimal settings
    set.seed(rngseed) # Reproducible randomness
    final_model = randomForest(X_train_fold, y_train_fold, 
                               mtry = rf_model$bestTune$mtry,
                               splitrule = rf_model$bestTune$splitrule,
                               min.node.size = rf_model$bestTune$min.node.size,
                               importance = TRUE)
    
    
    
    # Calculate and save model statistics (TRAINING DATA)
    
    # ------- !!! change below here
    
    train_pred<-predict(final_model,X_train_fold,type = "prob")
    
  
    test_pred<-predict(final_model,X_test_fold,type = "prob")

    
    
    #train_pred_proba = predict(final_model, type = "prob")[, 2]
    #train_auc = auc(roc(y_train_fold, train_pred_proba))
    #train_auc_scores = c(train_auc_scores, train_auc)
    # Save predictions
    temp = tibble(row = c(1:nrow(X))[-fold], # Row from the original training dataset
                  true_labels = y_train_fold, # Actual outcomes
                  REF = train_pred[,1],
                  OM1 = train_pred[,2],
                  OM2 = train_pred[,3]) # Predicted outcomes
    all_labels_train[[length(all_labels_train) + 1]] = temp
    
    # Calculate and save model statistics (TESTING DATA)
    #test_pred = predict(final_model, X_test_fold, type = "prob")[, 2]
    #test_auc = auc(roc(y_test_fold, test_pred_proba))
    #test_auc_scores = c(test_auc_scores, test_auc)
    
    
    
    # Save predictions
    temp = tibble(row = c(1:nrow(X))[fold], # Row from the original testing dataset
                  true_labels = y_test_fold, # Actual outcomes
                  REF = test_pred[,1],
                  OM1 = test_pred[,2],
                  OM2 = test_pred[,3])
    all_labels_test[[length(all_labels_test) + 1]] = temp
    
    # Save feature importance values
    # Ctrl, PD: higher number = higher in that group
    # MeanDecrease: how does the model quality decrease if the variable is removed? Higher=more important
    feature_importance_values[[length(feature_importance_values) + 1]] = final_model$importance
    
    # ------- !!! change above here
    
  }
  
  # Average all the folds
  avg_result = average_rf(all_labels_train,all_labels_test,
                          feature_importance_values)
  
  return(avg_result)
}