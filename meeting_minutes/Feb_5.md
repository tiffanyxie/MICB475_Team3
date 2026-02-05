# Feb. 5, 2026: Brainstorm ideas for project 2

## Agenda
* Review current project ideas and finalize idea for project 2

## Current Project Ideas with Soil Dataset
[Paper associated with dataset](https://pubmed.ncbi.nlm.nih.gov/28765786/)

[Code for visualizing metadata metrics associated with project ideas](https://github.com/tiffanyxie/MICB475_Team3/tree/main/project_brainstorm)

**1. Comparing the impact of LTSP treatmemt on the microbiome within each ecozone**
* Remove samples with C1 and C2 compaction treatments, A horizon sampling depth, and herbicide treatment

| Ecozone | IDFBC | SBSBC | PPCA | BSON | JPON | LPTX |
|--------|-------|-------|------|------|------|------|
| Samples per treatment | 9 | 9 | 9 | 8–9 | 1 (OM3), 7–9 | 9–12 |
| \# of Treatment types | 3 (REF, OM1, OM2) | 3 (REF, OM1, OM2) | 4 | 4 | 4 | 4 |
| \# of Sites | 3 | 3 | 3 | 3 | 3 | 3 |
| Elevation Range (m) | 1075–1180 | 780–1100 | 1135–1350 | 442–450 | 228–490 | 88 |
| \# of Soil types | 1 | 3 | 1 | 2 | 3 | 1 |
| \# of Tree Cover Types | 3 | 3 | 1 | 1 | 3 | 1 |
| \# of Climatic Zone Types | 1 | 1 | 1 | 1 | 1 | 1 |
| Moisture Range | 53.7–88 | 40–85 | 6.7–60.7 | 29.7–76.4 | 24–73.7 | 9–41.1 |
| CN ratio Range | 25.1–50.9 | 28.7–42.3 | 23.6–46.4 | 30.9–47.4 | 29.8–40.7 | 23.2–32.3 |

Visualization of metadata metrics (sample number, elevation, moisture content, CN ratio)
![Metadata_Horizon](../project_brainstorm/metadata_compare_ltsp_treatment.png)

**2. Comparing the impact of soil depth on the microbiome within each ecozone**
* The only independent variable is horizon (O vs A)
* O horizon is the top organic layer, A horizon is the 20 cm of mineral layer below that
* Samples Removed: C1 and C2 compaction treatments, OM3, herbicide treatment

| Ecozone | IDFBC | SBSBC | PPCA | BSON | JPON | LPTX |
|--------|-------|-------|------|------|------|------|
| Samples per horizon | 27 | 27 | 27 | 25 | 24 | 33 |
| \# of Sites | 3 | 3 | 3 | 3 | 3 | 3 |
| Elevation Range | 1075–1180 | 780–1100 | 1135–1350 | 442–450 | 228–490 | 88 |
| \# of Soil types | 1 | 3 | 1 | 2 | 3 | 1 |
| \# of Tree Cover Types | 3 | 3 | 1 | 1 | 3 | 1 |
| \# of Climatic Zone Types | 1 | 1 | 1 | 1 | 1 | 1 |
| \# of Treatment types | 3 | 3 | 3 | 3 | 3 | 3 |
| Moisture Range | 20.6–88 | 15–85 | 6.7–60.7 | 13.1–75.7 | 4–73.7 | 4–41.1 |
| CN ratio Range | 25.1–50.9 | 17.4–42.3 | 19.1–46.4 | 16.3–47.4 | 11.1–40.7 | 15.5–32.3 |

Visualization of metadata metrics (sample number, elevation, moisture content, CN ratio)
![Metadata_Horizon](../project_brainstorm/metadata_compare_horizon.png)

**3. Correlate Diversity with some continuous metric**
* Pool nothing
* Measure alpha diversity for each individual sample
* Put that on a scatter plot with diversity on the y axis and CN ratio, Moisture content, elevation, or annual precipitation on the x axis and see if a linear relationship arises

## Meeting Minutes
Ideas 1 and 2:
* Original authors assessed effect of removal and compaction levels and found significant impacts
* Horizon has been assessed and more diversity in organic later
* Eliminate ideas 1 and 2
  
Idea 3: Explore multiple metadata and see what has the biggest impact
* Run diversity and statistical analyses and outputs p value of significant diversity

Machine Learning Model: Predict organic matter removal
* Build machine learning of what level of organic matter removal has been applied to a region
* What taxa are important for predicting to organic matter removal?
* Use core analyses to select taxa of importance and use that to train machine learning model
* Which parameter is best for predicting organic matter removal?
* Paper established OM1, OM2, OM3 have different microbial compositions
* Create decision tree to determine if it is organic matter removal 1, 2, or 3
* What metadata columns and taxonomic groups go in to influence organic matter removal
* Suggest 70 total samples for everything
* Need to filter before running

Relevance:
* Rogue soil dataset -> tells you what level of organic material removal is present
* Can be applied to natural organic removal
* Tells you which variables have greatest predictive power
  
Pipeline:
1. Diversity metric
2. DESeq, Indicator taxa, core microbio - identify taxa that go into model
3. Random forest model - requires tuning and set-up

## Action Items
* Proposal: ensure we can justify each aim, how will these analyses help answer our question
* Have an overview of the proposal (e.g. how subsetting data)
* Process dataset (QIIME2) everytime go to decision point get Bessie
* Leave full week to discuss machine learning model before presentation


