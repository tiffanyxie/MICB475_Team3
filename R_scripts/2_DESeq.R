#!/usr/bin/env Rscript
library(tidyverse)
library(phyloseq)
library(DESeq2)

setwd("~/Desktop/soil_proj")

#### Load data ####
load("phylo_soil.RData")

soil_glom <- tax_glom(phylo_soil, taxrank = "Genus")

#### DESeq ####


## NOTE: If you get a zeros error, then you need to add '1' count to all reads
soil_deseq_plus1 <- transform_sample_counts(soil_glom, function(x) x+1)
soil_deseq <- phyloseq_to_deseq2(soil_deseq_plus1, ~`LTSP.Treatment`)
DESEQ_soil <- DESeq(soil_deseq)

# comparisons
res_OM1_vs_REF <- results(DESEQ_soil, tidy=TRUE, 
                          contrast = c("LTSP.Treatment", "OM1", "REF"))
res_OM2_vs_REF <- results(DESEQ_soil, tidy=TRUE, 
                          contrast = c("LTSP.Treatment", "OM2", "REF"))
res_OM2_vs_OM1 <- results(DESEQ_soil, tidy=TRUE, 
                          contrast = c("LTSP.Treatment", "OM2", "OM1"))

## Volcano plot: effect size VS significance
res_OM1_vs_REF %>%
  filter(!is.na(padj)) %>%  # This line removes the rows causing the warning
  mutate(significant = padj<0.05 & abs(log2FoldChange)>1) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant)) +
  labs(title = "OM1 vs REF") +
  theme_minimal()

res_OM2_vs_REF %>%
  filter(!is.na(padj)) %>%  # This line removes the rows causing the warning
  mutate(significant = padj<0.05 & abs(log2FoldChange)>1) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant)) +
  labs(title = "OM2 vs REF") +
  theme_minimal()

res_OM2_vs_OM1 %>%
  filter(!is.na(padj)) %>%  # This line removes the rows causing the warning
  mutate(significant = padj<0.05 & abs(log2FoldChange)>1) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant)) +
  labs(title = "OM2 vs OM1") +
  theme_minimal()

# To get table of results
sigASVs_OM1_vs_REF <- res_OM1_vs_REF %>% 
  filter(padj<0.05 & abs(log2FoldChange)>1) %>%
  dplyr::rename(ASV=row)
sigASVs_OM2_vs_REF <- res_OM2_vs_REF %>% 
  filter(padj<0.05 & abs(log2FoldChange)>1) %>%
  dplyr::rename(ASV=row)
sigASVs_OM2_vs_OM1 <- res_OM2_vs_OM1 %>% 
  filter(padj<0.05 & abs(log2FoldChange)>1) %>%
  dplyr::rename(ASV=row)

# see how much genus
nrow(sigASVs_OM1_vs_REF) # 12
nrow(sigASVs_OM2_vs_REF) # 28
nrow(sigASVs_OM2_vs_OM1) # 17

# Get only asv names
sigASVs_vec_OM1_vs_REF <- sigASVs_OM1_vs_REF %>%
  pull(ASV)
sigASVs_vec_OM2_vs_REF <- sigASVs_OM2_vs_REF %>%
  pull(ASV)
sigASVs_vec_OM2_vs_OM1 <- sigASVs_OM2_vs_OM1 %>%
  pull(ASV)

# Prune phyloseq file
soil_deseq_OM1_vs_REF <- prune_taxa(sigASVs_vec_OM1_vs_REF,soil_glom)
soil_deseq_OM2_vs_REF <- prune_taxa(sigASVs_vec_OM2_vs_REF,soil_glom)
soil_deseq_OM2_vs_OM1 <- prune_taxa(sigASVs_vec_OM2_vs_OM1,soil_glom)

sigASVs_vec_OM1_vs_REF <- tax_table(soil_deseq_OM1_vs_REF ) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_OM1_vs_REF) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))
sigASVs_vec_OM2_vs_REF <- tax_table(soil_deseq_OM2_vs_REF ) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_OM2_vs_REF) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))
sigASVs_vec_OM2_vs_OM1 <- tax_table(soil_deseq_OM2_vs_OM1 ) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_OM2_vs_OM1) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))


plot_data_OM1_vs_REF <- tax_table(soil_deseq_OM1_vs_REF) %>% as.data.frame() %>%
  rownames_to_column(var = "ASV") %>%
  right_join(sigASVs_OM1_vs_REF) %>% 
  filter(!is.na(Genus)) %>% 
  mutate(Genus = str_replace(Genus, "g__", "")) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))
plot_data_OM2_vs_REF <- tax_table(soil_deseq_OM2_vs_REF) %>% as.data.frame() %>%
  rownames_to_column(var = "ASV") %>%
  right_join(sigASVs_OM2_vs_REF) %>% 
  filter(!is.na(Genus)) %>% 
  mutate(Genus = str_replace(Genus, "g__", "")) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))
plot_data_OM2_vs_OM1 <- tax_table(soil_deseq_OM2_vs_OM1) %>% as.data.frame() %>%
  rownames_to_column(var = "ASV") %>%
  right_join(sigASVs_OM2_vs_OM1) %>% 
  filter(!is.na(Genus)) %>% 
  mutate(Genus = str_replace(Genus, "g__", "")) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

ggplot(plot_data_OM1_vs_REF, aes(x = log2FoldChange, y = Genus, fill = log2FoldChange > 0)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(xmin = log2FoldChange - lfcSE, xmax = log2FoldChange + lfcSE), width = 0.2) +
  scale_fill_manual(values = c("TRUE" = "lightpink", "FALSE" = "lightblue"), 
                    labels = c("Lower in OM1", "Higher in OM1"),
                    name = "Abundance") +
  labs(title = "OM1 vs REF") +
  theme_minimal()
ggplot(plot_data_OM2_vs_REF, aes(x = log2FoldChange, y = Genus, fill = log2FoldChange > 0)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(xmin = log2FoldChange - lfcSE, xmax = log2FoldChange + lfcSE), width = 0.2) +
  scale_fill_manual(values = c("TRUE" = "lightpink", "FALSE" = "lightblue"), 
                    labels = c("Lower in OM2", "Higher in OM2"),
                    name = "Abundance") +
  labs(title = "OM2 vs REF") +
  theme_minimal()
ggplot(plot_data_OM2_vs_OM1, aes(x = log2FoldChange, y = Genus, fill = log2FoldChange > 0)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(xmin = log2FoldChange - lfcSE, xmax = log2FoldChange + lfcSE), width = 0.2) +
  scale_fill_manual(values = c("TRUE" = "lightpink", "FALSE" = "lightblue"), 
                    labels = c("Lower in OM2", "Higher in OM2"),
                    name = "Abundance") +
  labs(title = "OM2 vs OM1") +
  theme_minimal()
