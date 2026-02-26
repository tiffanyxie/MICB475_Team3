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
- So far we found that ecozones exhibit different alpha diversity, in order to parse out soil conditions that impact diversity, we will do sub-ecozone analyses for Aim 1 (use non-parametric tests)
- We can also compare soil metadata conditions between OM removal to check if our results are the same
- Goal for next week: finish Aim 1 by next week, Spearman correlation between each soil condition and alpha diversity, beta diversity, compare soil metadata conditions between OM removal. 

## Relevant Analysis Scripts
- [Diversity Stats trials.R](../soil_export/Diversity%20Stats%20trials.R)
- [RStudio_processing.R](../soil_export/RStudio_processing.R)
