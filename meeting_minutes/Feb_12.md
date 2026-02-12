# Feb 12, 2026: Data processing and proposal update

## Agenda
* Provide update on data processing
* Decide on metadata filtering 
* Review proposal

## Update on data processing
* Finished QIIME2 pipeline up to taxonomic analysis (generated taxonomy.qzv and taxa barplots)
* Relevant outputs: [QIIME2 Processing Lab Notebook Entry](../lab_notebook/1_QIIME2_Processing.md)
* Filtered table to exclude mitochondria and chloroplast, but need to decide what metadata categories to exclude

## Proposal update
* Skeleton of introduction/background and basis for experimental aims completed, and started on proposed approach
* How do we approach the hypothesis/prediction as this is very exploratory?

Aims:
1. Diversity metric: To determine what variables impact alpha diversity by screening all columns in the soil dataset, then sorting by p-value to look for significance.
2. DESeq, indicator taxa, core microbiota: To determine what microbial taxa is associated with OM treatment levels and can be used as predictive features for our machine learning model.
3. Random forest model: To develop a machine learning model to predict OM level using microbial taxa and soil condition factors, and to evaluate its accuracy.



## Meeting Notes

