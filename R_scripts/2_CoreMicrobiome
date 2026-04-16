#R script for Core Microbiome

                                                           
#### Load packages ####
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)
library(svglite)

#### Load data ####
load("phylo_soil_genus.RData")


#### Core Microbiome Analysis ####


# Convert to relative abundance 
soil_RA <- transform_sample_counts(phylo_soil_genus, 
                                   fun=function(x) x/sum(x))

# Subset samples to REF, OM1, and OM2
soil_REF <- subset_samples(soil_RA, `LTSP.Treatment`=="REF")
soil_OM1 <- subset_samples(soil_RA, `LTSP.Treatment`=="OM1")
soil_OM2 <- subset_samples(soil_RA, `LTSP.Treatment`=="OM2")


# Abundance thresholds to test 
#     - 0 (presence/absence)
#     - 0.001 (filter out rare ASVs)
#     - 0.01 (abundant ASVs)

# Core Microbiome Analysis: prevalence 50% and detection 0.1%
REF_core <- core_members(soil_REF, detection = 0.001, prevalence = 0.5)
OM1_core <- core_members(soil_OM1, detection = 0.001, prevalence = 0.5)
OM2_core <- core_members(soil_OM2, detection = 0.001, prevalence = 0.5)

# Combine core ASV results
soil_list_full <- list(REF = REF_core, OM1 = OM1_core, OM2 = OM2_core)

# Obtain genera associated with core ASVs
# Used AI to troubleshoot
REF_genus <- tax_table(soil_REF)[REF_core, "Genus"]
OM1_genus <- tax_table(soil_OM1)[OM1_core, "Genus"]
OM2_genus <- tax_table(soil_OM2)[OM2_core, "Genus"]

#### Results ####

# Generate venn diagram
soil_venn <- ggVennDiagram(list(
  REF = REF_genus,
  OM1 = OM1_genus,
  OM2 = OM2_genus))

#Used GenAi for troubleshooting formatting

# Plot venn diagram 
soil_venn +
  scale_fill_gradient(
    low = "#97e5f1",
    high = "#0292a7",
    name = "Genus Count") + 
  theme(legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))

# Save venn diagram
ggsave("figures/core_microbiome.png",
       units=c("in"),
       width = 5.5,
       height = 5)
ggsave("figures/core_microbiome.svg",
       units=c("in"),
       width = 5.5,
       height = 5)

# Obtain table of core genera associated with each section of venn diagram

ref_genera<-REF_genus %>% as.data.frame() %>% pull(Genus)
om1_genera<-OM1_genus %>% as.data.frame() %>% pull(Genus)
om2_genera<-OM2_genus %>% as.data.frame() %>% pull(Genus)

all_genus<- c(ref_genera,
              om1_genera,
              om2_genera) %>%
  unique()



core_genera_table<-data.frame(Genus = all_genus) %>%
  mutate(REF = Genus %in% ref_genera,
         OM1 = Genus %in% om1_genera,
         OM2 = Genus %in% om2_genera) %>%
  arrange(desc(REF),desc(OM1),desc(OM2))


REF_unique <-core_genera_table %>%
  filter(REF, !OM1, !OM2)

REF_OM1_unique <-core_genera_table %>%
  filter(REF, OM1, !OM2)

REF_OM2_unique <-core_genera_table %>%
  filter(REF, !OM1, OM2)

OM1_unique <-core_genera_table %>%
  filter(!REF, OM1, !OM2)

OM1_OM2_unique <-core_genera_table %>%
  filter(!REF, OM1, OM2)
