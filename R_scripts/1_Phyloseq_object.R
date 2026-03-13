# Create phyloseq object


#### Load libraries ####

library(phyloseq)
library(ape) # importing trees
library(tidyverse)
library(picante)

#### Load Data ####
otu <- read_delim(file = "feature-table.txt", delim = "\t", skip=1)
sampdat <- read_delim(file = "soil_metadata.txt", delim = "\t")
taxonomy <- read_delim(file = "taxonomy.tsv", delim = "\t")
phylotree <- read.tree("tree.nwk")


#### Preparing Data to Generate Phyloseq Object ####
#Format OTU table
otu_mat <- as.matrix(otu[,-1])
rownames(otu_mat) <- otu$`#OTU ID`
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)

#Format Sample Metadata
samp_df <- as.data.frame(sampdat[,-1])
rownames(samp_df) <- sampdat$`#SampleID`
#Replace pH 0 with NA
samp_df_clean <- samp_df %>% 
  mutate(pH = if_else(pH == 0, NA, pH))
SAMP <- sample_data(samp_df_clean)

#Format taxonomy table
tax_mat <- taxonomy %>%
  select(-Confidence) %>%
  separate(col=Taxon, sep = ";",
           into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix()
tax_mat <- tax_mat[,-1]
rownames(tax_mat) <- taxonomy$`Feature ID`
TAX <- tax_table(tax_mat)

#### Create phyloseq object####
phylo_soil <- phyloseq(OTU, SAMP, TAX, phylotree)
#Rarefy to a depth of 3310
phylo_soil_rare <- rarefy_even_depth(phylo_soil, rngseed = 67, sample.size = 3310)
# Glom to genus level
phylo_soil_genus<-phylo_soil %>% tax_glom('Genus')

#Check if pH successfully changed
get_variable(phylo_soil,"pH")
#Check other metadata variables
get_variable(phylo_soil,"Total.Carbon")
get_variable(phylo_soil,"Total.Nitrogen")
get_variable(phylo_soil,"CN.Ratio")
get_variable(phylo_soil,"Moisture.Content")


#### Save phyloseq objects ####
save(phylo_soil, file = "phylo_soil.RData")
save(phylo_soil_rare, file = "phylo_soil_rare.RData")
save(phylo_soil_genus, file = "phylo_soil_genus.RData")

load("phylo_soil.RData")
load("phylo_soil_rare.RData")
load("phylo_soil_genus.RData")
