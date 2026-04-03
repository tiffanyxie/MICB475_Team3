# Load packages
library(phyloseq)
library(indicspecies)
library(dplyr)
library(readr)
library(ggplot2)
library(permute)
library(patchwork)
library(svglite)

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
  plot_df %>% mutate(GenusLabel = gsub("g__","",GenusLabel)),
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


#### Format Figures ####

plot_df_unique<-plot_df %>% 
  mutate(GenusLabel = make.unique(gsub("g__","",GenusLabel))) %>%
  arrange(desc(stat)) %>%
  mutate(GenusLabel = factor(GenusLabel,levels=unique(GenusLabel)))

ref<-plot_df_unique %>%
  filter(AssociatedTreatment == "REF") %>%
  ggplot(aes(x = GenusLabel, y = stat)) +
  geom_bar(stat="identity") +
  ylab("Indicator Value") +
  theme_classic() +
  ggtitle("REF") +
  scale_y_continuous(expand=c(0,0),limits=c(0,0.8),breaks=seq(0,0.8,by=0.2)) +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = 12),
        axis.text = element_text(size=10),
        axis.title.x= element_blank(),
        axis.text.x =  element_text(angle = 45,hjust =1,size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11)) 
ref

om2<-plot_df_unique %>%
  filter(AssociatedTreatment == "OM2") %>%
  ggplot(aes(x = GenusLabel, y = stat)) +
  geom_bar(stat="identity") +
  ylab("Indicator Value") +
  theme_classic() +
  ggtitle("OM2") +
  scale_y_continuous(expand=c(0,0),limits=c(0,0.8),breaks=seq(0,0.8,by=0.2)) +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = 12),
        axis.text = element_text(size=10),
        axis.title.x= element_blank(),
        axis.text.x =  element_text(angle = 45,hjust =1,size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11))
om2

ref_om1<-plot_df_unique %>%
  filter(AssociatedTreatment == "REF+OM1") %>%
  ggplot(aes(x = GenusLabel, y = stat)) +
  geom_bar(stat="identity") +
  ylab("Indicator Value") +
  theme_classic() +
  ggtitle("REF + OM1") +
  scale_y_continuous(expand=c(0,0),limits=c(0,0.8),breaks=seq(0,0.8,by=0.2)) +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = 12),
        axis.text = element_text(size=10),
        axis.title.x= element_blank(),
        axis.text.x =  element_text(angle = 45,hjust =1,size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11))
ref_om1 

ref_om2<-plot_df_unique %>%
  filter(AssociatedTreatment == "REF+OM2") %>%
  ggplot(aes(x = GenusLabel, y = stat)) +
  geom_bar(stat="identity") +
  ylab("Indicator Value") +
  theme_classic() +
  ggtitle("REF + OM2") +
  scale_y_continuous(expand=c(0,0),limits=c(0,0.8),breaks=seq(0,0.8,by=0.2)) +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = 12),
        axis.text = element_text(size=10),
        axis.title.x= element_blank(),
        axis.text.x =  element_text(angle = 45,hjust =1,size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11))
ref_om2 

om1_om2<-plot_df_unique %>%
  filter(AssociatedTreatment == "OM1+OM2") %>%
  ggplot(aes(x = GenusLabel, y = stat)) +
  geom_bar(stat="identity") +
  ylab("Indicator Value") +
  theme_classic() +
  ggtitle("OM1 + OM2") +
  scale_y_continuous(expand=c(0,0),limits=c(0,0.8),breaks=seq(0,0.8,by=0.2)) +
  theme(plot.title = element_blank(),
        axis.title = element_text(size = 12),
        axis.title.x= element_blank(),
        axis.text = element_text(size=10),
        axis.text.x =  element_text(angle = 45,hjust =1,size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11))
om1_om2 


ref + plot_spacer() + om2 + plot_layout(ncol = 3, widths=c(4,0.1,2))
ggsave("figures/isa_individual_treatment.png",width=6,height=3,units=c('in'))

ref_om1 + plot_spacer() + ref_om2 + plot_spacer() + om1_om2  + 
  plot_layout(ncol = 5, widths=c(3,0.15,2,0.15,2),guides="collect")
ggsave("figures/isa_combined_treatment.png",width=6,height=3,units=c('in'))


all_om<-plot_df_unique %>%
  ggplot(aes(x = GenusLabel, y = stat,fill = AssociatedTreatment)) +
  geom_bar(stat="identity") +
  ylab("Indicator Value") +
  xlab("Genus") +
  labs(fill = "OM Treatment") +
  theme_classic() +
  scale_y_continuous(expand=c(0,0),limits=c(0,0.8),breaks=seq(0,0.8,by=0.2)) +
  theme(plot.title = element_text(size = 14,hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size=10),
        axis.text.x =  element_text(angle = 45,hjust =1,size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11))
all_om
ggsave("figures/Figure_2_ISA_combined.png",width=6,height=3,units=c('in'))


#### Original Plot ####
p <- ggplot(
  plot_df %>% mutate(GenusLabel = gsub("g__","",GenusLabel)),
  aes(
    x = reorder(GenusLabel, stat),
    y = stat,
    fill = AssociatedTreatment
  )
) +
  geom_col() +
  facet_grid(
    . ~ AssociatedTreatment,
    scales = "free_y",
    space = "free_y"
  ) +
  labs(
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
