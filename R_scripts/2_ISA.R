# ISA by ecozone using unrarefied data
# only keep REF, OM1, OM2
# keep significant taxa (p <= 0.05)

# Load packages
library(phyloseq)
library(indicspecies)
library(dplyr)
library(readr)
library(ggplot2)
library(permute)

# Set seed for reproducibility
set.seed(67)

# Use the UNRAREFIED phyloseq object
ps <- phylo_soil

# Keep only samples from the selected ecozone
do_isa <- function(ps, ez) {
keep <- sample_names(ps)[sample_data(ps)$Ecozone == ez]
ps_ez <- prune_samples(keep, ps)
  
# Keep only REF, OM1, OM2
keep_trt <- sample_names(ps_ez)[sample_data(ps_ez)$LTSP.Treatment %in% c("REF", "OM1", "OM2")]
ps_ez <- prune_samples(keep_trt, ps_ez)
  
# Remove taxa with zero abundance after subsetting
ps_ez <- prune_taxa(taxa_sums(ps_ez) > 0, ps_ez)
  
# Aggregate taxa at the genus level
ps_ez <- tax_glom(ps_ez, taxrank = "Genus")
  
# Convert counts to relative abundance within each sample
ps_ez <- transform_sample_counts(ps_ez, function(x) {
  if (sum(x) == 0) return(x)
  x / sum(x)
})
  
# Extract OTU table as matrix for ISA
otu <- as(otu_table(ps_ez), "matrix")
if (taxa_are_rows(ps_ez)) otu <- t(otu)
  
# Define treatment groups to each sample
grp <- as.factor(as.character(sample_data(ps_ez)$LTSP.Treatment))
  
# Run indicator species analysis
isa <- multipatt(otu, grp, func = "IndVal.g", control = how(nperm = 999))
  
# Convert output to dataframe
out <- as.data.frame(isa$sign)
out$taxon <- rownames(out)
out$Ecozone <- ez
  
# Join taxonomy table so that Genus names can be used
tax_df <- as.data.frame(tax_table(ps_ez))
tax_df$taxon <- rownames(tax_df)
out <- left_join(out, tax_df, by = "taxon")
  
# Create genus label for plotting
out$GenusLabel <- ifelse(
  is.na(out$Genus) | out$Genus == "" | out$Genus == "NA",
  out$taxon,
  as.character(out$Genus)
)
  
  # Keep significant taxa only
  out <- out %>%
    filter(p.value <= 0.05)
  
  return(out)
}

# Run ISA separately for each ecozone, and combine into one data chart ready for plot
ecoz <- unique(as.character(sample_data(ps)$Ecozone))

res <- bind_rows(lapply(ecoz, function(ez) do_isa(ps, ez)))

# Save result table
write_csv(res, "ISA_indicators_by_Ecozone_REF_OM1_OM2_unrarefied_p0.05.csv")

# Create treatment label from ISA result columns
res2 <- res %>%
  mutate(
    AssociatedTreatment = case_when(
      s.REF == 1 & s.OM1 == 0 & s.OM2 == 0 ~ "REF",
      s.REF == 0 & s.OM1 == 1 & s.OM2 == 0 ~ "OM1",
      s.REF == 0 & s.OM1 == 0 & s.OM2 == 1 ~ "OM2",
      s.REF == 1 & s.OM1 == 1 & s.OM2 == 0 ~ "REF+OM1",
      s.REF == 1 & s.OM1 == 0 & s.OM2 == 1 ~ "REF+OM2",
      s.REF == 0 & s.OM1 == 1 & s.OM2 == 1 ~ "OM1+OM2",
      s.REF == 1 & s.OM1 == 1 & s.OM2 == 1 ~ "REF+OM1+OM2",
      TRUE ~ "Other"
    )
  )

# Keep only single-treatment indicators
plot_df <- res2 %>%
  filter(AssociatedTreatment %in% c("REF", "OM1", "OM2")) %>%
  group_by(Ecozone, AssociatedTreatment) %>%
  slice_max(order_by = stat, n = 5, with_ties = FALSE) %>%
  ungroup()

# Create plot
p <- ggplot(
  plot_df,
  aes(
    x = reorder(GenusLabel, stat),
    y = stat,
    fill = AssociatedTreatment
  )
) +
  geom_col() +
  coord_flip() +
  facet_grid(
    Ecozone ~ AssociatedTreatment,
    scales = "free_y",
    space = "free_y"
  ) +
  labs(
    title = "Top indicator genera by ecozone and treatment",
    x = "Genus",
    y = "Indicator value"
  ) +
  theme_classic() +
  scale_y_continuous(n.breaks = 4) +
  theme(
    text = element_text(size = 12),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 10),
    legend.position = "none"
  )

# Show plot
print(p)

# Save plot
ggsave(
  "ISA_indicator_plot_by_ecozone_and_treatment.png",
  plot = p,
  width = 12,
  height = 8
)
