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

## The below code is to identify which genera are present in each section of the Venn Diagram

#Make the soil list into a dataframe of Feature IDs
core_table <- data.frame(REF = as.vector(soil_list_full$REF), 
                         OM1 = as.vector(soil_list_full$OM1),
                         OM2 = c(as.vector(soil_list_full$OM2), NA)) #added a value to make the vector lengths equal

#Put Feature IDs into a single column to allow for joining to taxonomy table
core_long <- core_table %>%
  pivot_longer(cols = all_of(c("REF", "OM1", "OM2")),
               names_to = "Treatment",
               values_to = "Feature ID") %>%
  arrange(Treatment) %>%
  filter(`Feature ID` != is.na(`Feature ID`))

#Join dataframe to taxonomy table
tax_dat <- read_delim(file = "taxonomy.tsv", delim = "\t") %>%
  select(-Confidence) %>%
  separate(col=Taxon, sep = ";",
           into = c("Domain","Phylum","Class","Order","Family","Genus","Species"))

core_tax <- inner_join(core_long, tax_dat)

#Change data frame to report on the treatments each taxonomic category is present in
core_in_treatments <- core_tax%>%
  select(Treatment, Genus) %>% #can change Genus to whatever taxonomic level is of interest
  distinct() %>%                # removes duplicate treatment–genus pairs
  mutate(present = Treatment) %>%
  pivot_wider(
    names_from = Genus, #make sure to change the taxonomic level here as well
    values_from = present,
    values_fill = NA
  )
view(core_in_treatments)

#Same as above, but with raw ASVs instead of taxonomic levels
ASVs_in_treatments <- core_long%>%
  distinct() %>%                # removes duplicate treatment–genus pairs
  mutate(present = Treatment) %>%
  pivot_wider(
    names_from = `Feature ID`,
    values_from = present,
    values_fill = NA
  )
view(ASVs_in_treatments)
