library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)

#### Load data ####
#Use non-rarefied phyloseq

load("phylo_soil.RData")

#### "core" microbiome ####

## Convert to relative abundance 
soil_RA <- transform_sample_counts(phylo_soil, fun=function(x) x/sum(x))

## Filter dataset by LTSP Treatment  
soil_REF <- subset_samples(soil_RA, `LTSP.Treatment`=="OM1")
soil_OM1 <- subset_samples(soil_RA, `LTSP.Treatment`=="OM1")
soil_OM2 <- subset_samples(soil_RA, `LTSP.Treatment`=="OM2")


#For later reference: Abundance thresholds to test - 0 (presence/absence), 0.001 (filter out rare ASVs), 0.01 (abundant ASVs)

#Note that prevalence set to 0.2 did not give any ASVs. Prevalence was arbitrarily set to 0.15 ## 

## What ASVs are found in more than 15% of samples in each LTSP Treatment group?

soil_REF_ASVs <- core_members(soil_REF, detection=0, prevalence = 0.15)
soil_OM1_ASVs <- core_members(soil_OM1, detection=0, prevalence = 0.15)
soil_OM2_ASVs <- core_members(soil_OM2, detection=0, prevalence = 0.15)

## What are these ASVs? 
tax_table(prune_taxa(soil_REF_ASVs,phylo_soil))
tax_table(prune_taxa(soil_OM1_ASVs,phylo_soil))
tax_table(prune_taxa(soil_OM2_ASVs,phylo_soil))

## Plot those ASVs' relative abundance 

#In REF, it seems that the relative abundance is very low for the genus of interest?

prune_taxa(soil_REF_ASVs,soil_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`LTSP.Treatment`, scales ="free")

prune_taxa(soil_OM1_ASVs,soil_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`LTSP.Treatment`, scales ="free")

prune_taxa(soil_OM2_ASVs,soil_RA) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`LTSP.Treatment`, scales ="free")

## Create a Venn diagram using all the ASVs shared and unique to LTSP treatment groups 
soil_list_full <- list(REF = soil_REF_ASVs, OM1 = soil_OM1_ASVs, OM2 = soil_OM2_ASVs)

soil_venn <- ggVennDiagram(x = soil_list_full)

## Save venn diagram
ggsave("soil_venn.png", soil_venn)
