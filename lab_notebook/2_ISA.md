# Indicator Species Analysis (ISA)

## Aim:

- Identify microbial taxa that are significantly associated with LTSP treatment groups (REF, OM1, OM2) or treatment combinations across BC as a whole after merging ecozones.
- Indicator Species Analysis (ISA) was used to detect taxa whose relative abundance patterns are significantly associated with specific treatment groups or treatment combinations.

## Code:

[R Script](../R_scripts/2_ISA.R)

## Results:

Indicator species analysis was performed on all BC samples combined, rather than separately by ecozone, using genus-level relative abundance data.  
Only samples from the REF, OM1, and OM2 treatment groups were retained for the analysis.

The ISA algorithm (`multipatt` from the **indicspecies** package) was used to identify genera whose occurrence and relative abundance were significantly associated with specific treatment groups or treatment combinations.

A total of **23 significant indicator taxa** (**p < 0.05**) were detected across BC.

The indicator value (`stat`) reflects the strength of association between a taxon and a given treatment group or treatment combination, but taxa were retained based on **statistical significance (p < 0.05)** rather than by selecting only the highest indicator values.

### Indicator taxa summary

![ISA indicator taxa](../R_scripts/output/2_ISA.png)

The figure shows **all significant indicator taxa across BC** to provide an overall overview of treatment associations.  
Unlike the previous version, the plot does not select only the top taxa by indicator value, and ties are retained rather than removed.  
Higher indicator values indicate a stronger association between a taxon and a specific LTSP treatment group or treatment combination.

adjustment has been made compared to version last week: 
1. merge the ecozones and set bc as a whole
2. keep p< 0.05, not the indicator score
3. generate graph with all data to have an overview
4. keep the tie situation instead of filtering them out
