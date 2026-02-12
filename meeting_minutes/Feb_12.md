# Feb 12, 2026: Data processing and proposal update

## Agenda
* Provide update on data processing
* Decide on metadata filtering 
* Review proposal

## Update on data processing
* Finished QIIME2 pipeline up to taxonomic analysis (generated taxonomy.qzv and taxa barplots)
* Relevant outputs: [QIIME2 Processing Lab Notebook Entry](../lab_notebook/1_QIIME2_Processing.md)
* Filtered table to exclude mitochondria and chloroplast, but need to decide what metadata categories to exclude (likely only focus on organic layer?)

## Proposal update
* Skeleton of introduction/background and basis for experimental aims completed, and started on proposed approach
* How do we approach the hypothesis/prediction as this is very exploratory?

Aims:
1. Diversity metric: To determine what variables impact alpha diversity by screening all columns in the soil dataset, then sorting by p-value to look for significance.
2. DESeq, indicator taxa, core microbiota: To determine what microbial taxa is associated with OM treatment levels and can be used as predictive features for our machine learning model.
3. Random forest model: To develop a machine learning model to predict OM level using microbial taxa and soil condition factors, and to evaluate its accuracy.

Questions about random forest/methodology
* Confirmation: will we be able to build a model to compare > 2 levels of organic matter removal (Ref, OM1, OM2, OM3)?
  * If a pivot is necessary, would predicting soil depth be a worthwhile alternative?
* How much will we have to subset our table before feeding it into the random forest model?
* Confirmation: Remove herbicide use and A horizons? Remove compaction treatments?
* Samples for testing: exclude from DESeq, indicator species analysis, and core microbiome? or only exclude from random forest?

## Meeting Notes

