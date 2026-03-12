library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)

#### Load data ####
#Use non-rarefied phyloseq

load("phylo_soil.RData")

#### "core" microbiome ####

## Glom to Genus level
phylo_soil <- tax_glom(phylo_soil, taxrank = "Genus", NArm=FALSE)


## Convert to relative abundance 
soil_RA <- transform_sample_counts(phylo_soil, fun=function(x) x/sum(x))

## Filter dataset by LTSP Treatment  
soil_REF <- subset_samples(soil_RA, `LTSP.Treatment`=="REF")
soil_OM1 <- subset_samples(soil_RA, `LTSP.Treatment`=="OM1")
soil_OM2 <- subset_samples(soil_RA, `LTSP.Treatment`=="OM2")


#For reference: Abundance thresholds to test - 0 (presence/absence), 0.001 (filter out rare ASVs), 0.01 (abundant ASVs)


## What ASVs are found in more than 90% of samples in each LTSP Treatment group?

soil_REF_ASVs <- core_members(soil_REF, detection=0.001, prevalence = 0.9)
soil_OM1_ASVs <- core_members(soil_OM1, detection=0.001, prevalence = 0.9)
soil_OM2_ASVs <- core_members(soil_OM2, detection=0.001, prevalence = 0.9)


## What are these ASVs? 
tax_table(prune_taxa(soil_REF_ASVs,phylo_soil))
tax_table(prune_taxa(soil_OM1_ASVs,phylo_soil))
tax_table(prune_taxa(soil_OM2_ASVs,phylo_soil))



## Plot those ASVs' relative abundance 

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
#ggsave("soil_venn.png", soil_venn)

