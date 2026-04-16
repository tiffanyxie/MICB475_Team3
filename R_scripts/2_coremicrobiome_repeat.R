library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)
library(svglite)

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

## What FeatureIDs are found in more than 50% of samples in each LTSP Treatment group?

REF_core <- core_members(soil_REF, detection = 0.001, prevalence = 0.5)
OM1_core <- core_members(soil_OM1, detection = 0.001, prevalence = 0.5)
OM2_core <- core_members(soil_OM2, detection = 0.001, prevalence = 0.5)

soil_list_full <- list(REF = REF_core, OM1 = OM1_core, OM2 = OM2_core)

#used AI to trouble shoot to get genus level 
REF_genus <- tax_table(soil_REF)[REF_core, "Genus"]
OM1_genus <- tax_table(soil_OM1)[OM1_core, "Genus"]
OM2_genus <- tax_table(soil_OM2)[OM2_core, "Genus"]

soil_venn <- ggVennDiagram(list(
  REF = REF_genus,
  OM1 = OM1_genus,
  OM2 = OM2_genus))

#Used GenAi for troubleshooting formatting

soil_venn +
  scale_fill_gradient(
    low = "#97e5f1",
    high = "#0292a7",
    name = "Genus Count") + 
  theme(legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))

## Save venn diagram
ggsave("figures/core_microbiome.png",
       units=c("in"),
       width = 5.5,
       height = 5)
ggsave("figures/core_microbiome.svg",
       units=c("in"),
       width = 5.5,
       height = 5)

## The below code is to identify which genera are present in each section of the Venn Diagram


#Trouble shoot using AI
# Make soil list into long format of Feature IDs (FIXED: no data.frame padding needed)
core_long <- bind_rows(
  data.frame(FeatureID = soil_list_full$REF, Treatment = "REF"),
  data.frame(FeatureID = soil_list_full$OM1, Treatment = "OM1"),
  data.frame(FeatureID = soil_list_full$OM2, Treatment = "OM2")) %>%
  filter(!is.na(FeatureID))


# Join taxonomy table
tax_dat <- read_delim(file = "taxonomy.tsv", delim = "\t") %>%
  select(-Confidence) %>%
  rename(FeatureID = `Feature ID`) %>%
  separate(
    col = Taxon,
    sep = ";",
    into = c("Domain","Phylum","Class","Order","Family","Genus","Species"),
    fill = "right",
    remove = FALSE)

core_tax <- inner_join(core_long, tax_dat, by = "FeatureID")

#Change data frame to report on the treatments each taxonomic category is present in
core_in_treatments <- core_tax %>%
  select(Treatment, Genus) %>% #can change Genus to whatever taxonomic level is of interest
  filter(!is.na(Genus)) %>%
  distinct() %>% # removes duplicate treatment–genus pairs
  mutate(present = Treatment) %>%
  pivot_wider(
    names_from = Genus, #make sure to change the taxonomic level here as well
    values_from = present,
    values_fill = NA)

view(core_in_treatments)


#Same as above, but with raw ASVs instead of taxonomic levels
ASVs_in_treatments <- core_long %>%
  distinct(Treatment, FeatureID) %>% # removes duplicate treatment–genus pairs
  mutate(present = Treatment) %>%
  pivot_wider(
    names_from = FeatureID,
    values_from = present,
    values_fill = NA)

view(ASVs_in_treatments)



