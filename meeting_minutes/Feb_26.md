# Feb 26, 2026: Data processing and proposal update

## Agenda
1. Show Bessie the initial diversity analyses for Aim 1  
2. Review proposal  
3. Discuss how to infer importance of metadata categories from alpha diversity metrics  
   - (We decided in the proposal that Aim 1 will shift from “choosing metadata categories” to “scanning how metadata influences diversity.”)
   - Confirm that the statistical analysis methods chosen for aim 1 are correct (kruskal wallis, spearman, permanova) and confirm our conclusion of using non parametric tests

## Update on the project
- Began preparing R scripts for diversity analyses
- Continued refining the approach for Aim 1 (alpha diversity + metadata influence)
- 
[Aim 1 Lab Notebook](https://github.com/tiffanyxie/MICB475_Team3/blob/main/lab_notebook/RStudio_diversity.md)

## Discussion notes
- Aim 1: not directly assessing impact of OM removal, but just checking how metadata variables impact alpha and beta diversity
- We will include all soil conditions in our machine learning model
- Subset ecozones: do Aim 1 to compare soil conditions within each ecozone
- Do non-parametric tests for Aim 1 within each ecozone 

## Relevant Analysis Scripts
- [Diversity Stats trials.R](../soil_export/Diversity%20Stats%20trials.R)
- [RStudio_processing.R](../soil_export/RStudio_processing.R)
