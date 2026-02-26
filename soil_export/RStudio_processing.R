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
SAMP <- sample_data(samp_df)

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

save(phylo_soil, file = "phylo_soil.RData")
save(phylo_soil_rare, file = "phylo_soil_rare.RData")

load("phylo_soil.RData")
load("phylo_soil_rare.RData")

### Alpha diversity ####
#Prepare Metadata table to include Faith's PD
phylo_dist <- pd(t(otu_table(phylo_soil_rare)), phy_tree(phylo_soil_rare),
                 include.root=F) 
sample_data(phylo_soil_rare)$PD <- phylo_dist$PD


#Plot diversity by OM removal
plot_richness(phylo_soil_rare, x = "LTSP.Treatment", measures = c("Shannon","Chao1")) +
  xlab("OM Removal") +
  geom_boxplot()

ggplot(sample_data(phylo_soil_rare), aes(LTSP.Treatment, PD)) + 
  geom_boxplot() +
  xlab("OM Removal") +
  ylab("Phylogenetic Diversity")

#Plot diversity by Compaction
plot_richness(phylo_soil_rare, x = "Compaction.Treatment", measures = c("Shannon","Chao1")) +
  xlab("Compaction Treatment") +
  geom_boxplot()

ggplot(sample_data(phylo_soil_rare), aes(Compaction.Treatment, PD)) + 
  geom_boxplot() +
  xlab("Compaction") +
  ylab("Phylogenetic Diversity")

#Plot diversity by Site
plot_richness(phylo_soil_rare, x = "Site", measures = c("Shannon","Chao1")) +
  xlab("Site") +
  geom_boxplot()

ggplot(sample_data(phylo_soil_rare), aes(Site, PD)) + 
  geom_boxplot() +
  xlab("Site") +
  ylab("Phylogenetic Diversity")

#Plot diversity by Ecozone
plot_richness(phylo_soil_rare, x = "Ecozone", measures = c("Shannon","Chao1")) +
  xlab("Ecozone") +
  geom_boxplot()

ggplot(sample_data(phylo_soil_rare), aes(Ecozone, PD)) + 
  geom_boxplot() +
  xlab("Ecozone") +
  ylab("Phylogenetic Diversity")

#### Beta diversity #####

#Will need to revisit PERMANOVA usage for this tommorrow

#Obtain weighted unifrac distance metrics
wuni_dm <- distance(phylo_soil_rare, method = "wunifrac")
uni_dm <- distance(phylo_soil_rare, method = "unifrac")
bc_dm <- distance(phylo_soil_rare, method = "bray")

#Create coordinate tables
pcoa_wuni <- ordinate(phylo_soil_rare, method="PCoA", distance = wuni_dm)
pcoa_uni <- ordinate(phylo_soil_rare, method="PCoA", distance = uni_dm)
pcoa_bc <- ordinate(phylo_soil_rare, method="PCoA", distance = bc_dm)

pcoa_wuni_Comp <- ordinate(phylo_Comp, method="PCoA", distance = wuni_Comp)

#Generate figures
plot_ordination(phylo_soil_rare, pcoa_uni, 
                             color = "Ecozone", shape = "LTSP.Treatment",
                title = "Unweighted Unifrac Distance") +
  labs(pch="OM Removal", col = "Ecozone")

plot_ordination(phylo_soil_rare, pcoa_wuni, 
                color = "Ecozone", shape = "LTSP.Treatment",
                title = "Weighted Unifrac Distance") +
  labs(pch="OM Removal", col = "Ecozone")

plot_ordination(phylo_soil_rare, pcoa_bc, 
                color = "Ecozone", shape = "LTSP.Treatment",
                title = "Bray Curtis") +
  labs(pch="OM Removal", col = "Ecozone")

#### Taxonomy bar plots ####

# Convert to relative abundance
soil_RA <- transform_sample_counts(phylo_soil_rare, function(x) x/sum(x))

# To remove black bars, "glom" by phylum first
soil_phylum <- tax_glom(soil_RA, taxrank = "Phylum", NArm=FALSE)

plot_bar(soil_phylum, fill="Phylum") + 
  facet_wrap(.~LTSP.Treatment, scales = "free_x")

#ChatGPT-based code used below to identify abundant phyla
phylum_means <- psmelt(soil_phylum) %>%
  dplyr::group_by(Phylum) %>%
  dplyr::summarise(mean_abundance = mean(Abundance)) %>%
  dplyr::arrange(desc(mean_abundance))
head(phylum_means, 3)
