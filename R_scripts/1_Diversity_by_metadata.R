# R Script for Diversity Metrics

# Load libraries
library(phyloseq)
library(ape) # importing trees
library(tidyverse)
library(picante)
library(svglite)

load("phylo_soil.RData")
load("phylo_soil_rare.RData")

#Prepare Metadata table to include Faith's PD
phylo_dist <- pd(t(otu_table(phylo_soil_rare)), phy_tree(phylo_soil_rare),
                 include.root=F) 
sample_data(phylo_soil_rare)$PD <- phylo_dist$PD

#Create new Metadata table with other diversity metrics
alphadiv <- estimate_richness(phylo_soil_rare)
samp_dat <- sample_data(phylo_soil_rare)
samp_dat_wdiv <- data.frame(samp_dat, alphadiv)
#Pivot Diversity metrics to allow for faceting
dat_wdiv_long <- samp_dat_wdiv %>%
  select(-c(se.chao1, se.ACE, InvSimpson, Observed)) %>%
  pivot_longer(cols = all_of(c("PD", "Chao1", "ACE", "Shannon",
                               "Simpson", "Fisher")),
               names_to = "Metrics",
               values_to = "Metric_Values")


#Change OM treatment column to factor (used Claude AI for troubleshooting)
sample_data(phylo_soil_rare)$LTSP.Treatment<-factor(sample_data(phylo_soil_rare)$LTSP.Treatment,
                                                    levels=c("REF","OM1","OM2"))

#Subset data to specific ecozones
IDFBC_phylo_rare <- subset_samples(phylo_soil_rare, Ecozone == "IDFBC")
SBSBC_phylo_rare <- subset_samples(phylo_soil_rare, Ecozone == "SBSBC")

#Generate List of Diversity Metrics

metrics <- c("Shannon","Chao1", "ACE", "Simpson", "Fisher", "Observed","PD")



####Shannon's Diversity Figure####


kruskal.test(Shannon ~ `LTSP.Treatment`, data = subset(samp_dat_wdiv))

bc_shannon_plt<-plot_richness(phylo_soil_rare, x = "LTSP.Treatment", measures = c("Shannon")) +
  xlab("OM Removal") +
  ylab("Shannon's Diversity") +
  ggtitle("BC") +
  geom_boxplot() +
  theme_classic() +
  theme(plot.title = element_blank(), #element_text(size = 14,hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size=10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11),
        strip.text = element_blank(),
        strip.background = element_blank())
bc_shannon_plt
ggsave("figures/bc_diversity.png",units=c('in'),width = 4, height = 4)

#IDFBC
kruskal.test(Shannon ~ `LTSP.Treatment`, data = subset(samp_dat_wdiv, Ecozone == "IDFBC"))
idfbc_shannon_plt<-plot_richness(IDFBC_phylo_rare, x = "LTSP.Treatment", measures = c("Shannon")) +
  xlab("OM Removal") +
  ylab("Shannon's Diversity") +
  ggtitle("IDFBC") +
  geom_boxplot() +
  theme_classic() +
  theme(plot.title = element_blank(), #element_text(size = 14,hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size=10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11),
        strip.text = element_blank(), #Used claude to add code to remove facet labels
        strip.background = element_blank())

idfbc_shannon_plt
ggsave("figures/idfbc_diversity.png",units=c('in'),width = 4, height = 4)


#SBSBC
kruskal.test(Shannon ~ `LTSP.Treatment`, data = subset(samp_dat_wdiv, Ecozone == "SBSBC"))

sbsbc_shannon_plt<-plot_richness(SBSBC_phylo_rare, x = "LTSP.Treatment", measures = c("Shannon")) +
  xlab("OM Removal") +
  ylab("Shannon's Diversity") +
  ggtitle("SBSBC") +
  geom_boxplot() +
  theme_classic() +
  theme(plot.title = element_blank(), #element_text(size = 14,hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size=10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11),
        strip.text = element_blank(),
        strip.background = element_blank())
sbsbc_shannon_plt
ggsave("figures/sbsbc_diversity.png",units=c('in'),width = 4, height = 4)

#### Shannon Diversity Plot's Matching Axis ####

bc_shannon_plt<-plot_richness(phylo_soil_rare, x = "LTSP.Treatment", measures = c("Shannon")) +
  xlab("OM Removal") +
  ylab("Shannon's Diversity") +
  ggtitle("BC") +
  geom_boxplot() +
  theme_classic() +
  theme(plot.title = element_blank(), #element_text(size = 14,hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size=10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11),
        strip.text = element_blank(),
        strip.background = element_blank()) +
  scale_y_continuous(breaks=c(3.5,4,4.5,5),limits = c(3.5,5.1))

idfbc_shannon_plt<-plot_richness(IDFBC_phylo_rare, x = "LTSP.Treatment", measures = c("Shannon")) +
  xlab("OM Removal") +
  ylab("Shannon's Diversity") +
  ggtitle("IDFBC") +
  geom_boxplot() +
  theme_classic() +
  theme(plot.title = element_blank(), #element_text(size = 14,hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size=10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11),
        strip.text = element_blank(), #Used claude to add code to remove facet labels
        strip.background = element_blank()) +
  scale_y_continuous(breaks=c(3.5,4,4.5,5),limits = c(3.5,5.1))

sbsbc_shannon_plt<-plot_richness(SBSBC_phylo_rare, x = "LTSP.Treatment", measures = c("Shannon")) +
  xlab("OM Removal") +
  ylab("Shannon's Diversity") +
  ggtitle("SBSBC") +
  geom_boxplot() +
  theme_classic() +
  theme(plot.title = element_blank(), #element_text(size = 14,hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size=10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11),
        strip.text = element_blank(),
        strip.background = element_blank()) +
  scale_y_continuous(breaks=c(3.5,4,4.5,5),limits = c(3.5,5.1))

bc_shannon_plt + idfbc_shannon_plt + sbsbc_shannon_plt 
ggsave("figures/combined_diversity.png",width = 8, height = 4, units=c('in'))

#### Plot alpha diversity by OM removal####
#IDFBC Ecozone
plot_richness(IDFBC_phylo_rare, x = "LTSP.Treatment", nrow = 2, measures = c("Shannon","Chao1", "ACE", "Simpson", "Fisher", "Observed"),
              title = "IDFBC Alpha by OMR") +
  xlab("OM Removal") +
  geom_boxplot()

ggplot(sample_data(IDFBC_phylo_rare), aes(LTSP.Treatment, PD)) + 
  geom_boxplot() +
  labs(title = "IDFBC Faith's by OMR", x = "OM Removal", y = "Phylogenetic Diversity")

kruskal.test(Observed ~ `LTSP.Treatment`, data = subset(samp_dat_wdiv, Ecozone == "IDFBC"))
kruskal.test(Chao1 ~ `LTSP.Treatment`, data = subset(samp_dat_wdiv, Ecozone == "IDFBC"))
kruskal.test(ACE ~ `LTSP.Treatment`, data = subset(samp_dat_wdiv, Ecozone == "IDFBC"))
kruskal.test(Shannon ~ `LTSP.Treatment`, data = subset(samp_dat_wdiv, Ecozone == "IDFBC"))
kruskal.test(Simpson ~ `LTSP.Treatment`, data = subset(samp_dat_wdiv, Ecozone == "IDFBC"))
kruskal.test(Fisher ~ `LTSP.Treatment`, data = subset(samp_dat_wdiv, Ecozone == "IDFBC"))
kruskal.test(PD ~ `LTSP.Treatment`, data = subset(samp_dat_wdiv, Ecozone == "IDFBC"))

#SBSBC Ecozone
plot_richness(SBSBC_phylo_rare, x = "LTSP.Treatment", nrow = 2, measures = c("Shannon","Chao1", "ACE", "Simpson", "Fisher", "Observed"),
              title = "SBSBC Alpha by OMR") +
  xlab("OM Removal") +
  geom_boxplot()

ggplot(sample_data(SBSBC_phylo_rare), aes(LTSP.Treatment, PD)) + 
  geom_boxplot() +
  labs(title = "SBSBC Faith's by OMR", x = "OM Removal", y = "Phylogenetic Diversity")

kruskal.test(Observed ~ `LTSP.Treatment`, data = subset(samp_dat_wdiv, Ecozone == "SBSBC"))
kruskal.test(Chao1 ~ `LTSP.Treatment`, data = subset(samp_dat_wdiv, Ecozone == "SBSBC"))
kruskal.test(ACE ~ `LTSP.Treatment`, data = subset(samp_dat_wdiv, Ecozone == "SBSBC"))
kruskal.test(Shannon ~ `LTSP.Treatment`, data = subset(samp_dat_wdiv, Ecozone == "SBSBC"))
kruskal.test(Simpson ~ `LTSP.Treatment`, data = subset(samp_dat_wdiv, Ecozone == "SBSBC"))
kruskal.test(Fisher ~ `LTSP.Treatment`, data = subset(samp_dat_wdiv, Ecozone == "SBSBC"))
kruskal.test(PD ~ `LTSP.Treatment`, data = subset(samp_dat_wdiv, Ecozone == "SBSBC"))

#All of BC
plot_richness(phylo_soil_rare, x = "LTSP.Treatment", nrow = 2, measures = c("Shannon","Chao1", "ACE", "Simpson", "Fisher", "Observed"),
              title = "All of BC's Alpha by OMR") +
  xlab("OM Removal") +
  geom_boxplot()

ggplot(sample_data(phylo_soil_rare), aes(LTSP.Treatment, PD)) + 
  geom_boxplot() +
  labs(title = "All BC Faith's by OMR", x = "OM Removal", y = "Phylogenetic Diversity")

kruskal.test(Observed ~ `LTSP.Treatment`, data = samp_dat_wdiv)
kruskal.test(Chao1 ~ `LTSP.Treatment`, data = samp_dat_wdiv)
kruskal.test(ACE ~ `LTSP.Treatment`, data = samp_dat_wdiv)
kruskal.test(Shannon ~ `LTSP.Treatment`, data = samp_dat_wdiv)
kruskal.test(Simpson ~ `LTSP.Treatment`, data = samp_dat_wdiv)
kruskal.test(Fisher ~ `LTSP.Treatment`, data = samp_dat_wdiv)
kruskal.test(PD ~ `LTSP.Treatment`, data = samp_dat_wdiv)

#### Alpha diversity by Compaction####
#IDFBC Ecozone
plot_richness(IDFBC_phylo_rare, x = "Compaction.Treatment", nrow = 2, measures = c("Shannon","Chao1", "ACE", "Simpson", "Fisher", "Observed"),
              title = "IDFBC Alpha by Compaction") +
  xlab("Compaction Treatment") +
  geom_boxplot()

ggplot(sample_data(IDFBC_phylo_rare), aes(Compaction.Treatment, PD)) + 
  geom_boxplot() +
  labs(title = "IDFBC Faith's by Compaction", x = "Compaction Treatment", y = "Phylogenetic Diversity")

kruskal.test(Observed ~ `Compaction.Treatment`, data = subset(samp_dat_wdiv, Ecozone == "IDFBC"))
kruskal.test(Chao1 ~ `Compaction.Treatment`, data = subset(samp_dat_wdiv, Ecozone == "IDFBC"))
 anova_ob_vs_comp_log <- aov(lm(log(Chao1) ~ `Compaction.Treatment`, data=subset(samp_dat_wdiv, Ecozone == "IDFBC")))
 summary(anova_ob_vs_comp_log)
 TukeyHSD(anova_ob_vs_comp_log)
kruskal.test(ACE ~ `Compaction.Treatment`, data = subset(samp_dat_wdiv, Ecozone == "IDFBC"))
 anova_ob_vs_comp_log <- aov(lm(log(ACE) ~ `Compaction.Treatment`, data=subset(samp_dat_wdiv, Ecozone == "IDFBC")))
 summary(anova_ob_vs_comp_log)
 TukeyHSD(anova_ob_vs_comp_log)
kruskal.test(Shannon ~ `Compaction.Treatment`, data = subset(samp_dat_wdiv, Ecozone == "IDFBC"))
kruskal.test(Simpson ~ `Compaction.Treatment`, data = subset(samp_dat_wdiv, Ecozone == "IDFBC"))
kruskal.test(Fisher ~ `Compaction.Treatment`, data = subset(samp_dat_wdiv, Ecozone == "IDFBC"))
kruskal.test(PD ~ `Compaction.Treatment`, data = subset(samp_dat_wdiv, Ecozone == "IDFBC"))

#SBSBC Ecozone
plot_richness(SBSBC_phylo_rare, x = "Compaction.Treatment", nrow = 2, measures = c("Shannon","Chao1", "ACE", "Simpson", "Fisher", "Observed"),
              title = "SBSBC Alpha by Compaction") +
  xlab("Compaction Treatment") +
  geom_boxplot()

ggplot(sample_data(SBSBC_phylo_rare), aes(Compaction.Treatment, PD)) + 
  geom_boxplot() +
  labs(title = "SBSBC Faith's by Compaction", x = "Compaction Treatment", y = "Phylogenetic Diversity")

kruskal.test(Observed ~ `Compaction.Treatment`, data = subset(samp_dat_wdiv, Ecozone == "SBSBC"))
kruskal.test(Chao1 ~ `Compaction.Treatment`, data = subset(samp_dat_wdiv, Ecozone == "SBSBC"))
kruskal.test(ACE ~ `Compaction.Treatment`, data = subset(samp_dat_wdiv, Ecozone == "SBSBC"))
kruskal.test(Shannon ~ `Compaction.Treatment`, data = subset(samp_dat_wdiv, Ecozone == "SBSBC"))
kruskal.test(Simpson ~ `Compaction.Treatment`, data = subset(samp_dat_wdiv, Ecozone == "SBSBC"))
kruskal.test(Fisher ~ `Compaction.Treatment`, data = subset(samp_dat_wdiv, Ecozone == "SBSBC"))
kruskal.test(PD ~ `Compaction.Treatment`, data = subset(samp_dat_wdiv, Ecozone == "SBSBC"))

#All of BC
plot_richness(phylo_soil_rare, x = "Compaction.Treatment", nrow = 2, measures = c("Shannon","Chao1", "ACE", "Simpson", "Fisher", "Observed"),
              title = "All of BC's Alpha by Compaction") +
  xlab("Compaction Treatment") +
  geom_boxplot()

ggplot(sample_data(phylo_soil_rare), aes(Compaction.Treatment, PD)) + 
  geom_boxplot() +
  labs(title = "All BC Faith's by Compaction", x = "Compaction Treatment", y = "Phylogenetic Diversity")

kruskal.test(Observed ~ `Compaction.Treatment`, data = samp_dat_wdiv)
kruskal.test(Chao1 ~ `Compaction.Treatment`, data = samp_dat_wdiv)
kruskal.test(ACE ~ `Compaction.Treatment`, data = samp_dat_wdiv)
kruskal.test(Shannon ~ `Compaction.Treatment`, data = samp_dat_wdiv)
kruskal.test(Simpson ~ `Compaction.Treatment`, data = samp_dat_wdiv)
kruskal.test(Fisher ~ `Compaction.Treatment`, data = samp_dat_wdiv)
kruskal.test(PD ~ `Compaction.Treatment`, data = samp_dat_wdiv)


#### Alpha diversity by Site####
#IDFBC Ecozone
plot_richness(IDFBC_phylo_rare, x = "Site", nrow = 2, measures = c("Shannon","Chao1", "ACE", "Simpson", "Fisher", "Observed"),
              title = "IDFBC Alpha by Site") +
  xlab("Site") +
  geom_boxplot()

ggplot(sample_data(IDFBC_phylo_rare), aes(Site, PD)) + 
  geom_boxplot() +
  labs(title = "IDFBC Faith's by Site", x = "Site", y = "Phylogenetic Diversity")

kruskal.test(Observed ~ `Site`, data = subset(samp_dat_wdiv, Ecozone == "IDFBC"))
kruskal.test(Chao1 ~ `Site`, data = subset(samp_dat_wdiv, Ecozone == "IDFBC"))
kruskal.test(ACE ~ `Site`, data = subset(samp_dat_wdiv, Ecozone == "IDFBC"))
kruskal.test(Shannon ~ `Site`, data = subset(samp_dat_wdiv, Ecozone == "IDFBC"))
kruskal.test(Simpson ~ `Site`, data = subset(samp_dat_wdiv, Ecozone == "IDFBC"))
kruskal.test(Fisher ~ `Site`, data = subset(samp_dat_wdiv, Ecozone == "IDFBC"))
kruskal.test(PD ~ `Site`, data = subset(samp_dat_wdiv, Ecozone == "IDFBC"))
  anova_ob_vs_site_log <- aov(lm(log(PD) ~ `Site`, data=subset(samp_dat_wdiv, Ecozone == "IDFBC")))
  summary(anova_ob_vs_site_log)
  TukeyHSD(anova_ob_vs_site_log)


#SBSBC Ecozone
plot_richness(SBSBC_phylo_rare, x = "Site", nrow = 2, measures = c("Shannon","Chao1", "ACE", "Simpson", "Fisher", "Observed"),
              title = "SBSBC Alpha by Site") +
  xlab("Site") +
  geom_boxplot()

ggplot(sample_data(SBSBC_phylo_rare), aes(Site, PD)) + 
  geom_boxplot() +
  labs(title = "SBSBC Faith's by Site", x = "Site", y = "Phylogenetic Diversity")

kruskal.test(Observed ~ `Site`, data = subset(samp_dat_wdiv, Ecozone == "SBSBC"))
kruskal.test(Chao1 ~ `Site`, data = subset(samp_dat_wdiv, Ecozone == "SBSBC"))
kruskal.test(ACE ~ `Site`, data = subset(samp_dat_wdiv, Ecozone == "SBSBC"))
kruskal.test(Shannon ~ `Site`, data = subset(samp_dat_wdiv, Ecozone == "SBSBC"))
kruskal.test(Simpson ~ `Site`, data = subset(samp_dat_wdiv, Ecozone == "SBSBC"))
kruskal.test(Fisher ~ `Site`, data = subset(samp_dat_wdiv, Ecozone == "SBSBC"))
kruskal.test(PD ~ `Site`, data = subset(samp_dat_wdiv, Ecozone == "SBSBC"))

#All BC
plot_richness(phylo_soil_rare, x = "Site", nrow = 2, measures = c("Shannon","Chao1", "ACE", "Simpson", "Fisher", "Observed"),
              title = "All of BC's Alpha by Site") +
  xlab("Site") +
  geom_boxplot()

ggplot(sample_data(phylo_soil_rare), aes(Site, PD)) + 
  geom_boxplot() +
  labs(title = "All BC Faith's by Site", x = "Site", y = "Phylogenetic Diversity")

kruskal.test(Observed ~ `Site`, data = samp_dat_wdiv)
kruskal.test(Chao1 ~ `Site`, data = samp_dat_wdiv)
kruskal.test(ACE ~ `Site`, data = samp_dat_wdiv)
kruskal.test(Shannon ~ `Site`, data = samp_dat_wdiv)
kruskal.test(Simpson ~ `Site`, data = samp_dat_wdiv)
kruskal.test(Fisher ~ `Site`, data = samp_dat_wdiv)
kruskal.test(PD ~ `Site`, data = samp_dat_wdiv)

#### Alpha diversity by Moisture.Content####
#IDFBC Ecozone
ggplot(subset(dat_wdiv_long, Ecozone == "IDFBC"), aes(x=Moisture.Content, y=Metric_Values)) +
  geom_point() +
  labs(title = "IDFBC Alpha Diversity by Moisture", y = "") +
  geom_smooth(method = "lm")+ 
  facet_wrap(~Metrics, scales = "free_y")

cor.test(subset(samp_dat_wdiv, Ecozone == "IDFBC")$Moisture.Content, 
         subset(samp_dat_wdiv, Ecozone == "IDFBC")$Shannon, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, Ecozone == "IDFBC")$Moisture.Content, 
         subset(samp_dat_wdiv, Ecozone == "IDFBC")$Chao1, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, Ecozone == "IDFBC")$Moisture.Content, 
         subset(samp_dat_wdiv, Ecozone == "IDFBC")$ACE, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, Ecozone == "IDFBC")$Moisture.Content, 
         subset(samp_dat_wdiv, Ecozone == "IDFBC")$Simpson, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, Ecozone == "IDFBC")$Moisture.Content, 
         subset(samp_dat_wdiv, Ecozone == "IDFBC")$Fisher, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, Ecozone == "IDFBC")$Moisture.Content, 
         subset(samp_dat_wdiv, Ecozone == "IDFBC")$PD, method = "spearman", exact = FALSE)

#SBSBC Ecozone
ggplot(subset(dat_wdiv_long, Ecozone == "SBSBC"), aes(x=Moisture.Content, y=Metric_Values)) +
  geom_point() +
  labs(title = "SBSBC Alpha Diversity by Moisture", y = "") +
  geom_smooth(method = "lm")+ 
  facet_wrap(~Metrics, scales = "free_y")

cor.test(subset(samp_dat_wdiv, Ecozone == "SBSBC")$Moisture.Content, 
         subset(samp_dat_wdiv, Ecozone == "SBSBC")$Shannon, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, Ecozone == "SBSBC")$Moisture.Content, 
         subset(samp_dat_wdiv, Ecozone == "SBSBC")$Chao1, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, Ecozone == "SBSBC")$Moisture.Content, 
         subset(samp_dat_wdiv, Ecozone == "SBSBC")$ACE, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, Ecozone == "SBSBC")$Moisture.Content, 
         subset(samp_dat_wdiv, Ecozone == "SBSBC")$Simpson, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, Ecozone == "SBSBC")$Moisture.Content, 
         subset(samp_dat_wdiv, Ecozone == "SBSBC")$Fisher, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, Ecozone == "SBSBC")$Moisture.Content, 
         subset(samp_dat_wdiv, Ecozone == "SBSBC")$PD, method = "spearman", exact = FALSE)

#All BC
ggplot(subset(dat_wdiv_long), aes(x=Moisture.Content, y=Metric_Values)) +
  geom_point() +
  labs(title = "BC Alpha Diversity by Moisture", y = "") +
  geom_smooth(method = "lm")+ 
  facet_wrap(~Metrics, scales = "free_y")

cor.test(samp_dat_wdiv$Moisture.Content, 
         samp_dat_wdiv$Shannon, method = "spearman", exact = FALSE)

cor.test(samp_dat_wdiv$Moisture.Content, 
         samp_dat_wdiv$Chao1, method = "spearman", exact = FALSE)

cor.test(samp_dat_wdiv$Moisture.Content, 
         samp_dat_wdiv$ACE, method = "spearman", exact = FALSE)

cor.test(samp_dat_wdiv$Moisture.Content, 
         samp_dat_wdiv$Simpson, method = "spearman", exact = FALSE)

cor.test(samp_dat_wdiv$Moisture.Content, 
         samp_dat_wdiv$Fisher, method = "spearman", exact = FALSE)

cor.test(samp_dat_wdiv$Moisture.Content, 
         samp_dat_wdiv$PD, method = "spearman", exact = FALSE)

#### Alpha diversity by Total.Carbon####
#IDFBC Ecozone
ggplot(subset(dat_wdiv_long, Ecozone == "IDFBC"), aes(x=Total.Carbon, y=Metric_Values)) +
  geom_point() +
  labs(title = "IDFBC Alpha Diversity by Total.Carbon", y = "") +
  geom_smooth(method = "lm")+ 
  facet_wrap(~Metrics, scales = "free_y")

cor.test(subset(samp_dat_wdiv, Ecozone == "IDFBC")$Total.Carbon, 
         subset(samp_dat_wdiv, Ecozone == "IDFBC")$Shannon, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, Ecozone == "IDFBC")$Total.Carbon, 
         subset(samp_dat_wdiv, Ecozone == "IDFBC")$Chao1, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, Ecozone == "IDFBC")$Total.Carbon, 
         subset(samp_dat_wdiv, Ecozone == "IDFBC")$ACE, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, Ecozone == "IDFBC")$Total.Carbon, 
         subset(samp_dat_wdiv, Ecozone == "IDFBC")$Simpson, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, Ecozone == "IDFBC")$Total.Carbon, 
         subset(samp_dat_wdiv, Ecozone == "IDFBC")$Fisher, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, Ecozone == "IDFBC")$Total.Carbon, 
         subset(samp_dat_wdiv, Ecozone == "IDFBC")$PD, method = "spearman", exact = FALSE)

#SBSBC Ecozone
ggplot(subset(dat_wdiv_long, Ecozone == "SBSBC"), aes(x=Total.Carbon, y=Metric_Values)) +
  geom_point() +
  labs(title = "SBSBC Alpha Diversity by Total.Carbon", y = "") +
  geom_smooth(method = "lm")+ 
  facet_wrap(~Metrics, scales = "free_y")

cor.test(subset(samp_dat_wdiv, Ecozone == "SBSBC")$Total.Carbon, 
         subset(samp_dat_wdiv, Ecozone == "SBSBC")$Shannon, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, Ecozone == "SBSBC")$Total.Carbon, 
         subset(samp_dat_wdiv, Ecozone == "SBSBC")$Chao1, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, Ecozone == "SBSBC")$Total.Carbon, 
         subset(samp_dat_wdiv, Ecozone == "SBSBC")$ACE, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, Ecozone == "SBSBC")$Total.Carbon, 
         subset(samp_dat_wdiv, Ecozone == "SBSBC")$Simpson, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, Ecozone == "SBSBC")$Total.Carbon, 
         subset(samp_dat_wdiv, Ecozone == "SBSBC")$Fisher, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, Ecozone == "SBSBC")$Total.Carbon, 
         subset(samp_dat_wdiv, Ecozone == "SBSBC")$PD, method = "spearman", exact = FALSE)

#All BC
ggplot(subset(dat_wdiv_long), aes(x=Total.Carbon, y=Metric_Values)) +
  geom_point() +
  labs(title = "BC Alpha Diversity by Total.Carbon", y = "") +
  geom_smooth(method = "lm")+ 
  facet_wrap(~Metrics, scales = "free_y")

cor.test(samp_dat_wdiv$Total.Carbon, 
         samp_dat_wdiv$Shannon, method = "spearman", exact = FALSE)

cor.test(samp_dat_wdiv$Total.Carbon, 
         samp_dat_wdiv$Chao1, method = "spearman", exact = FALSE)

cor.test(samp_dat_wdiv$Total.Carbon, 
         samp_dat_wdiv$ACE, method = "spearman", exact = FALSE)

cor.test(samp_dat_wdiv$Total.Carbon, 
         samp_dat_wdiv$Simpson, method = "spearman", exact = FALSE)

cor.test(samp_dat_wdiv$Total.Carbon, 
         samp_dat_wdiv$Fisher, method = "spearman", exact = FALSE)

cor.test(samp_dat_wdiv$Total.Carbon, 
         samp_dat_wdiv$PD, method = "spearman", exact = FALSE)

#### Alpha diversity by Total.Nitrogen####
#IDFBC Ecozone
ggplot(subset(dat_wdiv_long, Ecozone == "IDFBC"), aes(x=Total.Nitrogen, y=Metric_Values)) +
  geom_point() +
  labs(title = "IDFBC Alpha Diversity by Total.Nitrogen", y = "") +
  geom_smooth(method = "lm")+ 
  facet_wrap(~Metrics, scales = "free_y")

cor.test(subset(samp_dat_wdiv, Ecozone == "IDFBC")$Total.Nitrogen, 
         subset(samp_dat_wdiv, Ecozone == "IDFBC")$Shannon, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, Ecozone == "IDFBC")$Total.Nitrogen, 
         subset(samp_dat_wdiv, Ecozone == "IDFBC")$Chao1, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, Ecozone == "IDFBC")$Total.Nitrogen, 
         subset(samp_dat_wdiv, Ecozone == "IDFBC")$ACE, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, Ecozone == "IDFBC")$Total.Nitrogen, 
         subset(samp_dat_wdiv, Ecozone == "IDFBC")$Simpson, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, Ecozone == "IDFBC")$Total.Nitrogen, 
         subset(samp_dat_wdiv, Ecozone == "IDFBC")$Fisher, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, Ecozone == "IDFBC")$Total.Nitrogen, 
         subset(samp_dat_wdiv, Ecozone == "IDFBC")$PD, method = "spearman", exact = FALSE)

#SBSBC Ecozone
ggplot(subset(dat_wdiv_long, Ecozone == "SBSBC"), aes(x=Total.Nitrogen, y=Metric_Values)) +
  geom_point() +
  labs(title = "SBSBC Alpha Diversity by Total.Nitrogen", y = "") +
  geom_smooth(method = "lm")+ 
  facet_wrap(~Metrics, scales = "free_y")

cor.test(subset(samp_dat_wdiv, Ecozone == "SBSBC")$Total.Nitrogen, 
         subset(samp_dat_wdiv, Ecozone == "SBSBC")$Shannon, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, Ecozone == "SBSBC")$Total.Nitrogen, 
         subset(samp_dat_wdiv, Ecozone == "SBSBC")$Chao1, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, Ecozone == "SBSBC")$Total.Nitrogen, 
         subset(samp_dat_wdiv, Ecozone == "SBSBC")$ACE, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, Ecozone == "SBSBC")$Total.Nitrogen, 
         subset(samp_dat_wdiv, Ecozone == "SBSBC")$Simpson, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, Ecozone == "SBSBC")$Total.Nitrogen, 
         subset(samp_dat_wdiv, Ecozone == "SBSBC")$Fisher, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, Ecozone == "SBSBC")$Total.Nitrogen, 
         subset(samp_dat_wdiv, Ecozone == "SBSBC")$PD, method = "spearman", exact = FALSE)

#All BC
ggplot(subset(dat_wdiv_long), aes(x=Total.Nitrogen, y=Metric_Values)) +
  geom_point() +
  labs(title = "BC Alpha Diversity by Total.Nitrogen", y = "") +
  geom_smooth(method = "lm")+ 
  facet_wrap(~Metrics, scales = "free_y")

cor.test(samp_dat_wdiv$Total.Nitrogen, 
         samp_dat_wdiv$Shannon, method = "spearman", exact = FALSE)

cor.test(samp_dat_wdiv$Total.Nitrogen, 
         samp_dat_wdiv$Chao1, method = "spearman", exact = FALSE)

cor.test(samp_dat_wdiv$Total.Nitrogen, 
         samp_dat_wdiv$ACE, method = "spearman", exact = FALSE)

cor.test(samp_dat_wdiv$Total.Nitrogen, 
         samp_dat_wdiv$Simpson, method = "spearman", exact = FALSE)

cor.test(samp_dat_wdiv$Total.Nitrogen, 
         samp_dat_wdiv$Fisher, method = "spearman", exact = FALSE)

cor.test(samp_dat_wdiv$Total.Nitrogen, 
         samp_dat_wdiv$PD, method = "spearman", exact = FALSE)

#### Alpha diversity by CN.Ratio####
#IDFBC Ecozone
ggplot(subset(dat_wdiv_long, Ecozone == "IDFBC"), aes(x=CN.Ratio, y=Metric_Values)) +
  geom_point() +
  labs(title = "IDFBC Alpha Diversity by CN.Ratio", y = "") +
  geom_smooth(method = "lm")+ 
  facet_wrap(~Metrics, scales = "free_y")

cor.test(subset(samp_dat_wdiv, Ecozone == "IDFBC")$CN.Ratio, 
         subset(samp_dat_wdiv, Ecozone == "IDFBC")$Shannon, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, Ecozone == "IDFBC")$CN.Ratio, 
         subset(samp_dat_wdiv, Ecozone == "IDFBC")$Chao1, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, Ecozone == "IDFBC")$CN.Ratio, 
         subset(samp_dat_wdiv, Ecozone == "IDFBC")$ACE, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, Ecozone == "IDFBC")$CN.Ratio, 
         subset(samp_dat_wdiv, Ecozone == "IDFBC")$Simpson, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, Ecozone == "IDFBC")$CN.Ratio, 
         subset(samp_dat_wdiv, Ecozone == "IDFBC")$Fisher, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, Ecozone == "IDFBC")$CN.Ratio, 
         subset(samp_dat_wdiv, Ecozone == "IDFBC")$PD, method = "spearman", exact = FALSE)

#SBSBC Ecozone
ggplot(subset(dat_wdiv_long, Ecozone == "SBSBC"), aes(x=CN.Ratio, y=Metric_Values)) +
  geom_point() +
  labs(title = "SBSBC Alpha Diversity by CN.Ratio", y = "") +
  geom_smooth(method = "lm")+ 
  facet_wrap(~Metrics, scales = "free_y")

cor.test(subset(samp_dat_wdiv, Ecozone == "SBSBC")$CN.Ratio, 
         subset(samp_dat_wdiv, Ecozone == "SBSBC")$Shannon, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, Ecozone == "SBSBC")$CN.Ratio, 
         subset(samp_dat_wdiv, Ecozone == "SBSBC")$Chao1, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, Ecozone == "SBSBC")$CN.Ratio, 
         subset(samp_dat_wdiv, Ecozone == "SBSBC")$ACE, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, Ecozone == "SBSBC")$CN.Ratio, 
         subset(samp_dat_wdiv, Ecozone == "SBSBC")$Simpson, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, Ecozone == "SBSBC")$CN.Ratio, 
         subset(samp_dat_wdiv, Ecozone == "SBSBC")$Fisher, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, Ecozone == "SBSBC")$CN.Ratio, 
         subset(samp_dat_wdiv, Ecozone == "SBSBC")$PD, method = "spearman", exact = FALSE)

#All BC
ggplot(subset(dat_wdiv_long), aes(x=CN.Ratio, y=Metric_Values)) +
  geom_point() +
  labs(title = "BC Alpha Diversity by CN.Ratio", y = "") +
  geom_smooth(method = "lm")+ 
  facet_wrap(~Metrics, scales = "free_y")

cor.test(samp_dat_wdiv$CN.Ratio, 
         samp_dat_wdiv$Shannon, method = "spearman", exact = FALSE)

cor.test(samp_dat_wdiv$CN.Ratio, 
         samp_dat_wdiv$Chao1, method = "spearman", exact = FALSE)

cor.test(samp_dat_wdiv$CN.Ratio, 
         samp_dat_wdiv$ACE, method = "spearman", exact = FALSE)

cor.test(samp_dat_wdiv$CN.Ratio, 
         samp_dat_wdiv$Simpson, method = "spearman", exact = FALSE)

cor.test(samp_dat_wdiv$CN.Ratio, 
         samp_dat_wdiv$Fisher, method = "spearman", exact = FALSE)

cor.test(samp_dat_wdiv$CN.Ratio, 
         samp_dat_wdiv$PD, method = "spearman", exact = FALSE)

#### Alpha diversity by pH####
#IDFBC Ecozone
subset(dat_wdiv_long, Ecozone == "IDFBC")%>%
  filter(pH > 0) %>% #Unlikely the pH was actually zero, most likely was not measured
  ggplot(aes(x=pH, y=Metric_Values)) +
  geom_point() +
  labs(title = "IDFBC Alpha Diversity by pH", y = "") +
  geom_smooth(method = "lm")+ 
  facet_wrap(~Metrics, scales = "free_y")

cor.test(subset(samp_dat_wdiv, Ecozone == "IDFBC" & pH > 0)$pH, 
         subset(samp_dat_wdiv, Ecozone == "IDFBC" & pH > 0)$Shannon, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, Ecozone == "IDFBC" & pH > 0)$pH, 
         subset(samp_dat_wdiv, Ecozone == "IDFBC" & pH > 0)$Chao1, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, Ecozone == "IDFBC" & pH > 0)$pH, 
         subset(samp_dat_wdiv, Ecozone == "IDFBC" & pH > 0)$ACE, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, Ecozone == "IDFBC" & pH > 0)$pH, 
         subset(samp_dat_wdiv, Ecozone == "IDFBC" & pH > 0)$Simpson, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, Ecozone == "IDFBC" & pH > 0)$pH, 
         subset(samp_dat_wdiv, Ecozone == "IDFBC" & pH > 0)$Fisher, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, Ecozone == "IDFBC" & pH > 0)$pH, 
         subset(samp_dat_wdiv, Ecozone == "IDFBC" & pH > 0)$PD, method = "spearman", exact = FALSE)

#SBSBC Ecozone
subset(dat_wdiv_long, Ecozone == "SBSBC")%>%
  filter(pH > 0) %>%
  ggplot(aes(x=pH, y=Metric_Values)) +
  geom_point() +
  labs(title = "SBSBC Alpha Diversity by pH", y = "") +
  geom_smooth(method = "lm")+ 
  facet_wrap(~Metrics, scales = "free_y")

cor.test(subset(samp_dat_wdiv, Ecozone == "SBSBC" & pH > 0)$pH, 
         subset(samp_dat_wdiv, Ecozone == "SBSBC" & pH > 0)$Shannon, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, Ecozone == "SBSBC" & pH > 0)$pH, 
         subset(samp_dat_wdiv, Ecozone == "SBSBC" & pH > 0)$Chao1, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, Ecozone == "SBSBC" & pH > 0)$pH, 
         subset(samp_dat_wdiv, Ecozone == "SBSBC" & pH > 0)$ACE, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, Ecozone == "SBSBC" & pH > 0)$pH, 
         subset(samp_dat_wdiv, Ecozone == "SBSBC" & pH > 0)$Simpson, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, Ecozone == "SBSBC" & pH > 0)$pH, 
         subset(samp_dat_wdiv, Ecozone == "SBSBC" & pH > 0)$Fisher, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, Ecozone == "SBSBC" & pH > 0)$pH, 
         subset(samp_dat_wdiv, Ecozone == "SBSBC" & pH > 0)$PD, method = "spearman", exact = FALSE)

#All BC
filter(dat_wdiv_long, pH > 0) %>%
  ggplot(aes(x=pH, y=Metric_Values)) +
  geom_point() +
  labs(title = "BC Alpha Diversity by pH", y = "") +
  geom_smooth(method = "lm")+ 
  facet_wrap(~Metrics, scales = "free_y")
  
filter(dat_wdiv_long, pH > 0) %>%
  ggplot(aes(y=pH, x=Ecozone)) +
  geom_boxplot() +
  labs(title = "pH by Ecozone", y = "")

cor.test(subset(samp_dat_wdiv, pH > 0)$pH, 
         subset(samp_dat_wdiv, pH > 0)$Shannon, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, pH > 0)$pH, 
         subset(samp_dat_wdiv, pH > 0)$Chao1, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, pH > 0)$pH, 
         subset(samp_dat_wdiv, pH > 0)$ACE, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, pH > 0)$pH, 
         subset(samp_dat_wdiv, pH > 0)$Simpson, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, pH > 0)$pH, 
         subset(samp_dat_wdiv, pH > 0)$Fisher, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, pH > 0)$pH, 
         subset(samp_dat_wdiv, pH > 0)$PD, method = "spearman", exact = FALSE)

#### Alpha diversity by Soil.Bulk.Density####
#IDFBC Ecozone
ggplot(subset(dat_wdiv_long, Ecozone == "IDFBC"), aes(x=Soil.Bulk.Density, y=Metric_Values)) +
  geom_point() +
  labs(title = "IDFBC Alpha Diversity by Soil.Bulk.Density", y = "") +
  geom_smooth(method = "lm")+ 
  facet_wrap(~Metrics, scales = "free_y")
  
cor.test(subset(samp_dat_wdiv, Ecozone == "IDFBC")$Soil.Bulk.Density, 
         subset(samp_dat_wdiv, Ecozone == "IDFBC")$Shannon, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, Ecozone == "IDFBC")$Soil.Bulk.Density, 
         subset(samp_dat_wdiv, Ecozone == "IDFBC")$Chao1, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, Ecozone == "IDFBC")$Soil.Bulk.Density, 
         subset(samp_dat_wdiv, Ecozone == "IDFBC")$ACE, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, Ecozone == "IDFBC")$Soil.Bulk.Density, 
         subset(samp_dat_wdiv, Ecozone == "IDFBC")$Simpson, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, Ecozone == "IDFBC")$Soil.Bulk.Density, 
         subset(samp_dat_wdiv, Ecozone == "IDFBC")$Fisher, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, Ecozone == "IDFBC")$Soil.Bulk.Density, 
         subset(samp_dat_wdiv, Ecozone == "IDFBC")$PD, method = "spearman", exact = FALSE)

#SBSBC Ecozone
ggplot(subset(dat_wdiv_long, Ecozone == "SBSBC"), aes(x=Soil.Bulk.Density, y=Metric_Values)) +
  geom_point() +
  labs(title = "SBSBC Alpha Diversity by Soil.Bulk.Density", y = "") +
  geom_smooth(method = "lm")+ 
  facet_wrap(~Metrics, scales = "free_y")
  
cor.test(subset(samp_dat_wdiv, Ecozone == "SBSBC")$Soil.Bulk.Density, 
         subset(samp_dat_wdiv, Ecozone == "SBSBC")$Shannon, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, Ecozone == "SBSBC")$Soil.Bulk.Density, 
         subset(samp_dat_wdiv, Ecozone == "SBSBC")$Chao1, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, Ecozone == "SBSBC")$Soil.Bulk.Density, 
         subset(samp_dat_wdiv, Ecozone == "SBSBC")$ACE, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, Ecozone == "SBSBC")$Soil.Bulk.Density, 
         subset(samp_dat_wdiv, Ecozone == "SBSBC")$Simpson, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, Ecozone == "SBSBC")$Soil.Bulk.Density, 
         subset(samp_dat_wdiv, Ecozone == "SBSBC")$Fisher, method = "spearman", exact = FALSE)

cor.test(subset(samp_dat_wdiv, Ecozone == "SBSBC")$Soil.Bulk.Density, 
         subset(samp_dat_wdiv, Ecozone == "SBSBC")$PD, method = "spearman", exact = FALSE)

#All BC
ggplot(subset(dat_wdiv_long), aes(x=Soil.Bulk.Density, y=Metric_Values)) +
  geom_point() +
  labs(title = "BC Alpha Diversity by Soil.Bulk.Density", y = "") +
  geom_smooth(method = "lm") + 
  facet_wrap(~Metrics, scales = "free_y")

cor.test(samp_dat_wdiv$Soil.Bulk.Density, 
         samp_dat_wdiv$Shannon, method = "spearman", exact = FALSE)

cor.test(samp_dat_wdiv$Soil.Bulk.Density, 
         samp_dat_wdiv$Chao1, method = "spearman", exact = FALSE)

cor.test(samp_dat_wdiv$Soil.Bulk.Density, 
         samp_dat_wdiv$ACE, method = "spearman", exact = FALSE)

cor.test(samp_dat_wdiv$Soil.Bulk.Density, 
         samp_dat_wdiv$Simpson, method = "spearman", exact = FALSE)

cor.test(samp_dat_wdiv$Soil.Bulk.Density, 
         samp_dat_wdiv$Fisher, method = "spearman", exact = FALSE)

cor.test(samp_dat_wdiv$Soil.Bulk.Density, 
         samp_dat_wdiv$PD, method = "spearman", exact = FALSE)

#### Generate Coordinate Tables for each subset to use in beta diversity plots####
#IDFBC Ecozone
IDFBC_jac <- ordinate(IDFBC_phylo_rare, method="PCoA", 
                      distance = distance(IDFBC_phylo_rare, method = "jaccard"))
IDFBC_bc <- ordinate(IDFBC_phylo_rare, method="PCoA", 
                     distance = distance(IDFBC_phylo_rare, method = "bray"))
IDFBC_uni <- ordinate(IDFBC_phylo_rare, method="PCoA", 
                      distance = distance(IDFBC_phylo_rare, method = "unifrac"))
IDFBC_wuni <- ordinate(IDFBC_phylo_rare, method="PCoA", 
                       distance = distance(IDFBC_phylo_rare, method = "wunifrac"))

#SBSBC Ecozone
SBSBC_jac <- ordinate(SBSBC_phylo_rare, method="PCoA", 
                      distance = distance(SBSBC_phylo_rare, method = "jaccard"))
SBSBC_bc <- ordinate(SBSBC_phylo_rare, method="PCoA", 
                     distance = distance(SBSBC_phylo_rare, method = "bray"))
SBSBC_uni <- ordinate(SBSBC_phylo_rare, method="PCoA", 
                      distance = distance(SBSBC_phylo_rare, method = "unifrac"))
SBSBC_wuni <- ordinate(SBSBC_phylo_rare, method="PCoA", 
                       distance = distance(SBSBC_phylo_rare, method = "wunifrac"))

#All BC
BC_jac <- ordinate(phylo_soil_rare, method="PCoA", 
                   distance = distance(phylo_soil_rare, method = "jaccard"))
BC_bc <- ordinate(phylo_soil_rare, method="PCoA", 
                  distance = distance(phylo_soil_rare, method = "bray"))
BC_uni <- ordinate(phylo_soil_rare, method="PCoA", 
                   distance = distance(phylo_soil_rare, method = "unifrac"))
BC_wuni <- ordinate(phylo_soil_rare, method="PCoA", 
                    distance = distance(phylo_soil_rare, method = "wunifrac"))

#### Plot beta diversity ####

#IDFBC
plot_ordination(IDFBC_phylo_rare, IDFBC_jac, 
                color = "Site", shape = "LTSP.Treatment",
                title = "Jaccard Index") +
  labs(pch="OM Removal", col = "Site")

plot_ordination(IDFBC_phylo_rare, IDFBC_bc, 
                color = "Site", shape = "LTSP.Treatment",
                title = "Bray Curtis Dissimilarity") +
  labs(pch="OM Removal", col = "Site")

plot_ordination(IDFBC_phylo_rare, IDFBC_uni, 
                color = "Site", shape = "LTSP.Treatment",
                title = "Unweighted Unifrac Distance") +
  labs(pch="OM Removal", col = "Site")

plot_ordination(IDFBC_phylo_rare, IDFBC_wuni, 
                color = "Site", shape = "LTSP.Treatment",
                title = "Weighted Unifrac Distance") +
  labs(pch="OM Removal", col = "Site")

#SBSBC
plot_ordination(SBSBC_phylo_rare, SBSBC_jac, 
                color = "Site", shape = "LTSP.Treatment",
                title = "Jaccard Index") +
  labs(pch="OM Removal", col = "Site")

plot_ordination(SBSBC_phylo_rare, SBSBC_bc, 
                color = "Site", shape = "LTSP.Treatment",
                title = "Bray Curtis Dissimilarity") +
  labs(pch="OM Removal", col = "Site")


plot_ordination(SBSBC_phylo_rare, SBSBC_uni, 
                color = "Site", shape = "LTSP.Treatment",
                title = "Unweighted Unifrac Distance") +
  labs(pch="OM Removal", col = "Site")

plot_ordination(SBSBC_phylo_rare, SBSBC_wuni, 
                color = "Site", shape = "LTSP.Treatment",
                title = "Weighted Unifrac Distance") +
  labs(pch="OM Removal", col = "Site")

#All BC
plot_ordination(phylo_soil_rare, BC_jac, 
                color = "Site", shape = "LTSP.Treatment",
                title = "Jaccard Index") +
  labs(pch="OM Removal", col = "Site")

plot_ordination(phylo_soil_rare, BC_bc, 
                color = "Site", shape = "LTSP.Treatment",
                title = "Bray Curtis Dissimilarity") +
  labs(pch="OM Removal", col = "Site")


plot_ordination(phylo_soil_rare, BC_uni, 
                color = "Site", shape = "LTSP.Treatment",
                title = "Unweighted Unifrac Distance") +
  labs(pch="OM Removal", col = "Site")

plot_ordination(phylo_soil_rare, BC_wuni, 
                color = "Site", shape = "LTSP.Treatment",
                title = "Weighted Unifrac Distance") +
  labs(pch="OM Removal", col = "Site")



#### Generate dsitance matrices for each subset to use with PERMANOVA ####
#IDFBC
IDFBC_dm_wuni <- distance(IDFBC_phylo_rare, method = "wunifrac")
IDFBC_dm_uni <- distance(IDFBC_phylo_rare, method = "unifrac")
IDFBC_dm_bray <- distance(IDFBC_phylo_rare, method = "bray")
IDFBC_dm_jac <- distance(IDFBC_phylo_rare, method = "jaccard")

#SBSBC
SBSBC_dm_wuni <- distance(SBSBC_phylo_rare, method = "wunifrac")
SBSBC_dm_uni <- distance(SBSBC_phylo_rare, method = "unifrac")
SBSBC_dm_bray <- distance(SBSBC_phylo_rare, method = "bray")
SBSBC_dm_jac <- distance(SBSBC_phylo_rare, method = "jaccard")

#ALl BC
BC_dm_wuni <- distance(phylo_soil_rare, method = "wunifrac")
BC_dm_uni <- distance(phylo_soil_rare, method = "unifrac")
BC_dm_bray <- distance(phylo_soil_rare, method = "bray")
BC_dm_jac <- distance(phylo_soil_rare, method = "jaccard")

#### Run PERMANOVA ####
#IDFBC
adonis2(IDFBC_dm_wuni ~ Site*LTSP.Treatment*Compaction.Treatment,
        data = subset(samp_dat_wdiv, Ecozone == "IDFBC"),
        by = "terms")

adonis2(IDFBC_dm_uni ~ Site*LTSP.Treatment*Compaction.Treatment,
        data = subset(samp_dat_wdiv, Ecozone == "IDFBC"),
        by = "terms")

adonis2(IDFBC_dm_bray ~ Site*LTSP.Treatment*Compaction.Treatment,
        data = subset(samp_dat_wdiv, Ecozone == "IDFBC"),
        by = "terms")

adonis2(IDFBC_dm_jac ~ Site*LTSP.Treatment*Compaction.Treatment,
        data = subset(samp_dat_wdiv, Ecozone == "IDFBC"),
        by = "terms")


#SBSBC
adonis2(SBSBC_dm_wuni ~ Site*LTSP.Treatment*Compaction.Treatment,
        data = subset(samp_dat_wdiv, Ecozone == "SBSBC"),
        by = "terms")

adonis2(SBSBC_dm_uni ~ Site*LTSP.Treatment*Compaction.Treatment,
        data = subset(samp_dat_wdiv, Ecozone == "SBSBC"),
        by = "terms")

adonis2(SBSBC_dm_bray ~ Site*LTSP.Treatment*Compaction.Treatment,
        data = subset(samp_dat_wdiv, Ecozone == "SBSBC"),
        by = "terms")

adonis2(SBSBC_dm_jac ~ Site*LTSP.Treatment*Compaction.Treatment,
        data = subset(samp_dat_wdiv, Ecozone == "SBSBC"),
        by = "terms")


#All BC
adonis2(BC_dm_wuni ~ Site*LTSP.Treatment*Compaction.Treatment,
        data = samp_dat_wdiv,
        by = "terms")

adonis2(BC_dm_uni ~ Site*LTSP.Treatment*Compaction.Treatment,
        data = samp_dat_wdiv,
        by = "terms")

adonis2(BC_dm_bray ~ Site*LTSP.Treatment*Compaction.Treatment,
        data = samp_dat_wdiv,
        by = "terms")

adonis2(BC_dm_jac ~ Site*LTSP.Treatment*Compaction.Treatment,
        data = samp_dat_wdiv,
        by = "terms")
