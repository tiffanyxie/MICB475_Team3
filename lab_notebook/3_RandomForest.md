# Aim 3: Random Forest

## Aim
* To develop a random forest model using soil conditions and bacterial populations to predict level of OM Removal

## Code
[Main Script: Random_Forest.R](https://github.com/tiffanyxie/MICB475_Team3/blob/main/R_scripts/3_RandomForest.R)

[Random Forest Functions](https://github.com/tiffanyxie/MICB475_Team3/blob/main/R_scripts/3_randomforest_functions_modified.R)

Modifications to Random_Forest.R
* Applied Bessie's suggestions to allow multiple outcomes
1) y = factor(y, levels = c("Control", "PD")) -> y = factor(y)
2) Remove classProbs = TRUE and summaryFunction = twoClassSummary from trainControl()
3) Change metric = "ROC" to metric = "Accuracy" in train()

* Added confusion matrix to results interpretation and no longer generating ROC Curves

Modifications to random forest functions
* No longer calculating AUC in run_rf()
* No longer returning AUC values in average_rf()
* Modify run_rf() and average_rf() to return true values and probability of REF, OM1, and OM2
* No change to importance values, still returning importance values

## Results

* Using total 37 bacterial genuses determined to be differentially abundant between any two OM treatment levels via DESeq (p < 0.05, log2FC > 1)
* Soil conditions: pH, total carbon, total nitrogen, CN ratio, soil moisture content
Hyperparameters:
tune_grid = expand.grid(mtry = c(3,6,10), 
                        splitrule = c("gini","extratrees"),
                        min.node.size = c(2,3,4))

Importance Values
![Importance Value Plot](https://github.com/tiffanyxie/MICB475_Team3/blob/main/R_scripts/output/importance_plot_1.png)

**Training**

Confusion Matrix

|            | **True REF** | **True OM1** | **True OM2** |
|------------|---------|---------|---------|
| **Predicted REF**    | 3       | 2       | 6       |
| **Predicted OM1**    | 8       | 22      | 20      |
| **Predicted OM2**    | 4       | 21      | 18      |


| **Metric**               | **REF** | **OM1** | **OM2** |
|--------------------------|---------|---------|---------|
| Sensitivity              | 0.20000 | 0.4889  | 0.4091  |
| Specificity              | 0.91011 | 0.5254  | 0.5833  |
| Pos Pred Value           | 0.27273 | 0.4400  | 0.4186  |
| Neg Pred Value           | 0.87097 | 0.5741  | 0.5738  |
| Prevalence               | 0.14423 | 0.4327  | 0.4231  |
| Detection Rate           | 0.02885 | 0.2115  | 0.1731  |
| Detection Prevalence     | 0.10577 | 0.4808  | 0.4135  |
| Balanced Accuracy        | 0.55506 | 0.5072  | 0.4962  |

**Testing**

Confusion Matrix

|             | **True REF** | **True OM1** | **True OM2** |
|-------------|---------|---------|---------|
| **Predicted REF**     | 10      | 0       | 1       |
| **Predicted OM1**     | 3       | 31      | 16      |
| **Predicted OM2**     | 2       | 14      | 27      |


| **Metric**               | **REF** | **OM1** | **OM2** |
|--------------------------|---------|---------|---------|
| Sensitivity              | 0.66667 | 0.6889  | 0.6136  |
| Specificity              | 0.98876 | 0.6780  | 0.7333  |
| Pos Pred Value           | 0.90909 | 0.6200  | 0.6279  |
| Neg Pred Value           | 0.94624 | 0.7407  | 0.7213  |
| Prevalence               | 0.14423 | 0.4327  | 0.4231  |
| Detection Rate           | 0.09615 | 0.2981  | 0.2596  |
| Detection Prevalence     | 0.10577 | 0.4808  | 0.4135  |
| Balanced Accuracy        | 0.82772 | 0.6834  | 0.6735  |

