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

# Run ISA on all BC samples together
do_isa <- function(ps) {
  
# Keep only REF, OM1, OM2
keep_trt <- sample_names(ps)[sample_data(ps)$LTSP.Treatment %in% c("REF", "OM1", "OM2")]
ps_sub <- prune_samples(keep_trt, ps)
  
# Remove taxa with zero abundance after subsetting
ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
  
# Aggregate taxa at the genus level
ps_sub <- tax_glom(ps_sub, taxrank = "Genus")
  
# Convert counts to relative abundance within each sample
ps_sub <- transform_sample_counts(ps_sub, function(x) {
  if (sum(x) == 0) return(x)
  x / sum(x)
})
  
# Extract OTU table as matrix for ISA
otu <- as(otu_table(ps_sub), "matrix")
if (taxa_are_rows(ps_sub)) otu <- t(otu)
  
# Define treatment groups for each sample
grp <- as.factor(as.character(sample_data(ps_sub)$LTSP.Treatment))
  
# Run indicator species analysis
isa <- multipatt(otu, grp, func = "IndVal.g", control = how(nperm = 999))
  
# Convert output to dataframe
out <- as.data.frame(isa$sign)
out$taxon <- rownames(out)
  
# Join taxonomy table so that Genus names can be used
tax_df <- as.data.frame(tax_table(ps_sub))
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
  filter(p.value < 0.05)
  
return(out)
}

# Run ISA on all BC data together
res <- do_isa(ps)

# Save result table
write_csv(res, "ISA_indicators_BC_REF_OM1_OM2_unrarefied_p_lt_0.05.csv")

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

# Keep all significant indicators for overview
plot_df <- res2 %>%
  filter(AssociatedTreatment != "Other")

# Optional: order treatment panels
plot_df$AssociatedTreatment <- factor(
  plot_df$AssociatedTreatment,
  levels = c("REF", "OM1", "OM2", "REF+OM1", "REF+OM2", "OM1+OM2", "REF+OM1+OM2")
)

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
    . ~ AssociatedTreatment,
    scales = "free_y",
    space = "free_y"
  ) +
  labs(
    title = "Significant indicator genera across BC by treatment combination",
    x = "Genus",
    y = "Indicator value"
  ) +
  theme_classic() +
  scale_y_continuous(
    breaks = c(0, 0.2, 0.4, 0.6),
    labels = c("0.0", "0.2", "0.4", "0.6")
  ) +
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
  "ISA_indicator_plot_BC_all_significant_by_treatment.png",
  plot = p,
  width = 16,
  height = 8
)
