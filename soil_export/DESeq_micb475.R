#!/usr/bin/env Rscript
library(tidyverse)
library(phyloseq)
library(DESeq2)
library(stringr)

#### Load data ####
load("phylo_soil.RData")

#### DESeq ####
soil_deseq_plus1 <- transform_sample_counts(phylo_soil, function(x) x+1)
soil_deseq <- phyloseq_to_deseq2(soil_deseq_plus1, ~`LTSP.Treatment`)
DESEQ_soil <- DESeq(soil_deseq)

# comparisons
res_OM1_vs_REF <- results(DESEQ_soil, tidy=TRUE, 
                          contrast = c("LTSP.Treatment", "OM1", "REF"))
res_OM2_vs_REF <- results(DESEQ_soil, tidy=TRUE, 
                          contrast = c("LTSP.Treatment", "OM2", "REF"))
res_OM2_vs_OM1 <- results(DESEQ_soil, tidy=TRUE, 
                          contrast = c("LTSP.Treatment", "OM2", "OM1"))

## volcano plot: effect size vs significance
res_OM1_vs_REF %>%
  filter(!is.na(padj)) %>%  # This line removes the rows causing the warning
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant)) +
  labs(title = "OM1 vs REF") +
  theme_minimal()

res_OM2_vs_REF %>%
  filter(!is.na(padj)) %>%  # This line removes the rows causing the warning
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant)) +
  labs(title = "OM2 vs REF") +
  theme_minimal()

res_OM2_vs_OM1 %>%
  filter(!is.na(padj)) %>%  # This line removes the rows causing the warning
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant)) +
  labs(title = "OM2 vs OM1") +
  theme_minimal()


# to get table of results
sigASVs_OM1_vs_REF <- res_OM1_vs_REF %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
sigASVs_OM2_vs_REF <- res_OM2_vs_REF %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)
sigASVs_OM2_vs_OM1 <- res_OM2_vs_OM1 %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)

# get only asv names
sigASVs_vec_OM1_vs_REF <- sigASVs_OM1_vs_REF %>%
  pull(ASV)
sigASVs_vec_OM2_vs_REF <- sigASVs_OM2_vs_REF %>%
  pull(ASV)
sigASVs_vec_OM2_vs_OM1 <- sigASVs_OM2_vs_OM1 %>%
  pull(ASV)

# prune phyloseq file
soil_deseq_OM1_vs_REF <- prune_taxa(sigASVs_vec_OM1_vs_REF,phylo_soil)
soil_deseq_OM2_vs_REF <- prune_taxa(sigASVs_vec_OM2_vs_REF,phylo_soil)
soil_deseq_OM2_vs_OM1 <- prune_taxa(sigASVs_vec_OM2_vs_OM1,phylo_soil)

sigASVs_final_OM1_vs_REF <- tax_table(soil_deseq_OM1_vs_REF) %>% 
  as.data.frame() %>%
  rownames_to_column(var = "ASV") %>%
  right_join(sigASVs_OM1_vs_REF, by = "ASV")
sigASVs_final_OM2_vs_REF <- tax_table(soil_deseq_OM2_vs_REF) %>% 
  as.data.frame() %>%
  rownames_to_column(var = "ASV") %>%
  right_join(sigASVs_OM2_vs_REF, by = "ASV")
sigASVs_final_OM2_vs_OM1 <- tax_table(soil_deseq_OM2_vs_OM1) %>% 
  as.data.frame() %>%
  rownames_to_column(var = "ASV") %>%
  right_join(sigASVs_OM2_vs_OM1, by = "ASV")

# clean, filter for Top 20, and set sorting
plot_data_OM1_vs_REF <- sigASVs_final_OM1_vs_REF %>%
  # remove NAs and the 'g__' prefix
  filter(!is.na(Genus) & Genus != "NA") %>% 
  mutate(Genus = str_replace(Genus, "g__", "")) %>%
  # handle duplicates (e.g. multiple Afipia strains)
  mutate(Genus = make.unique(Genus)) %>%
  # pick top 20 by significance
  filter(!is.na(padj)) %>%
  arrange(padj) %>% 
  head(20) %>%
  # re-sort
  arrange(log2FoldChange) %>% 
  mutate(Genus = factor(Genus, levels = unique(Genus)))

plot_data_OM2_vs_REF <- sigASVs_final_OM2_vs_REF %>%
  filter(!is.na(Genus) & Genus != "NA") %>% 
  mutate(Genus = str_replace(Genus, "g__", "")) %>%
  mutate(Genus = make.unique(Genus)) %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>% 
  head(20) %>%
  arrange(log2FoldChange) %>% 
  mutate(Genus = factor(Genus, levels = unique(Genus)))

plot_data_OM2_vs_OM1 <- sigASVs_final_OM2_vs_OM1 %>%
  # remove NAs and the 'g__' prefix
  filter(!is.na(Genus) & Genus != "NA") %>% 
  mutate(Genus = str_replace(Genus, "g__", "")) %>%
  # handle duplicates (e.g. multiple Afipia strains)
  mutate(Genus = make.unique(Genus)) %>%
  # pick top 20 by significance
  filter(!is.na(padj)) %>%
  arrange(padj) %>% 
  head(20) %>%
  # re-sort
  arrange(log2FoldChange) %>% 
  mutate(Genus = factor(Genus, levels = unique(Genus)))

# create bar plots
ggplot(plot_data_OM1_vs_REF, aes(x = log2FoldChange, y = Genus, fill = log2FoldChange > 0)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  geom_errorbar(aes(xmin = log2FoldChange - lfcSE, xmax = log2FoldChange + lfcSE), width = 0.2) +
  scale_fill_manual(values = c("TRUE" = "lightpink", "FALSE" = "lightblue"), 
                    labels = c("Lower in OM1", "Higher in OM1"),
                    name = "Abundance") +
  labs(title = "OM1 vs REF: Top 20 Significant Genera",
       subtitle = "Grouped by Abundance Direction") +
  theme_minimal()

ggplot(plot_data_OM2_vs_REF, aes(x = log2FoldChange, y = Genus, fill = log2FoldChange > 0)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  geom_errorbar(aes(xmin = log2FoldChange - lfcSE, xmax = log2FoldChange + lfcSE), width = 0.2) +
  scale_fill_manual(values = c("TRUE" = "lightpink", "FALSE" = "lightblue"), 
                    labels = c("Lower in OM2", "Higher in OM1"),
                    name = "Abundance") +
  labs(title = "OM2 vs REF: Top 20 Significant Genera",
       subtitle = "Grouped by Abundance Direction") +
  theme_minimal()

ggplot(plot_data_OM2_vs_OM1, aes(x = log2FoldChange, y = Genus, fill = log2FoldChange > 0)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  geom_errorbar(aes(xmin = log2FoldChange - lfcSE, xmax = log2FoldChange + lfcSE), width = 0.2) +
  scale_fill_manual(values = c("TRUE" = "lightpink", "FALSE" = "lightblue"), 
                    labels = c("Lower in OM2", "Higher in OM1"),
                    name = "Abundance") +
  labs(title = "OM2 vs OM1: Top 20 Significant Genera",
       subtitle = "Grouped by Abundance Direction") +
  theme_minimal()