# March 19, 2026: Aim 2 ISA Update and Aim 3 Random Forest First Attempt

## Agenda 
1. Discuss ISA results
2. Discuss first attempt at generating RF model :D

## Updates on Project

Random Forest Updates!
* Modified random forest script to generate initial results with confusion matrix
* Want to check with Bessie and team to ensure our approach to modifying code is correct before diving into model results
* What metrics will we use to assess the accuracy of our model (ROC, confusion matrix)?

[Notes on Random Forest Code Modification and Initial Results](../lab_notebook/3_RandomForest.md) \
 [Random Forest Main Script](../R_scripts/3_RandomForest.R) \
[Random Forest Functions](../R_scripts/3_randomforest_functions_modified.R)

## Discussion Notes
* ISA results not particularly helpful
* Multi-class model is not very good at differentiating OM1 and OM2
* Total carbon and CN ratio are important predictive factors
* Can we do a binary model to differentiate OM1 vs REF, OM2 vs REF,and OM1 + OM2 vs REF?
* Generate ROC curves and calcluate AUC
