# Load packages needed for phyloseq data handling and indicator species analysis
library(phyloseq)
library(indicspecies)
library(dplyr)
library(readr)

# Set seed so results are reproducible
set.seed(67)

# Use the rarefied phyloseq object created earlier
ps <- phylo_soil_rare


# Function to run Indicator Species Analysis (ISA) for one ecozone
do_isa <- function(ps, ez){
  
  # Keep only samples that belong to the selected ecozone
  keep <- sample_names(ps)[ sample_data(ps)$Ecozone == ez ]
  ps_ez <- prune_samples(keep, ps)
  
  # Remove taxa with zero abundance after subsetting
  ps_ez <- prune_taxa(taxa_sums(ps_ez) > 0, ps_ez)
  
  # Aggregate taxa at the genus level
  ps_ez <- tax_glom(ps_ez, taxrank = "Genus")
  
  # Convert counts to relative abundance within each sample
  ps_ez <- transform_sample_counts(ps_ez, function(x) x / sum(x))
  
  # Extract OTU table as matrix for ISA
  otu <- as(otu_table(ps_ez), "matrix")
  if (taxa_are_rows(ps_ez)) otu <- t(otu)
  
  # Define treatment groups
  grp <- as.factor(sample_data(ps_ez)$LTSP.Treatment)
  
  # Run indicator species analysis using the indicspecies package
  isa <- multipatt(otu, grp, func = "r.g", control = how(nperm = 999))
  
  # Convert results to a dataframe
  out <- as.data.frame(isa$sign)
  out$taxon <- rownames(out)
  out$Ecozone <- ez
  
  # Keep only significant indicator taxa
  out <- dplyr::filter(out, p.value <= 0.05)
  
  out
}


# Get all ecozones present in the dataset
ecoz <- unique(as.character(sample_data(ps)$Ecozone))

# Run ISA separately for each ecozone and combine the results
res <- bind_rows(lapply(ecoz, function(ez) do_isa(ps, ez)))

# Display results
res

# Save indicator taxa results to a CSV file
write_csv(res, "ISA_indicators_by_Ecozone.csv")

# Plot indicator value for significant taxa

library(ggplot2)

ggplot(res, aes(x = reorder(taxon, stat), y = stat, fill = Ecozone)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(
    title = "Indicator taxa by ecozone",
    x = "Taxon",
    y = "Indicator value"
  ) +
  theme_minimal()

ggsave("~/Desktop/ISA_indicator_plot.png", width = 8, height = 6)