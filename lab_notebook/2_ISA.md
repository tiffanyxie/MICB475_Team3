# Indicator Species Analysis (ISA)

## Aim:

- Identify microbial taxa that are significantly associated with LTSP treatment groups (REF, OM1, OM2) within each ecozone.
- Indicator Species Analysis (ISA) was used to detect taxa whose relative abundance patterns are strongly associated with specific treatments.

## Code:

[R Script](../R_scripts/2_ISA.R)

## Results:

Indicator species analysis was performed separately for each ecozone using genus-level relative abundance data.  
The ISA algorithm (`multipatt` from the **indicspecies** package) identifies taxa whose occurrence and abundance are significantly associated with particular treatment groups.

A total of **58 significant indicator taxa** (p ≤ 0.05) were detected across ecozones.

The indicator value (statistic) reflects the strength of association between a taxon and a specific treatment group.

### Indicator taxa summary

![ISA indicator taxa](../R_scripts/output/ISA_plot.png)

The figure shows the **top 10 taxa with the highest indicator values in each ecozone** for visualization purposes.  
Higher indicator values indicate a stronger association between a taxon and a specific LTSP treatment group.


still have to work on: 
- within each ecozone, also have to analyze the OM 1,2,3 removal
- change the name of top 10 species into genus name or smth
- change to use the unrarefied, compare om1 to rarefied, and om2 to rarefied, the value use 0.6
