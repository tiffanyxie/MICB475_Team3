#!/usr/bin/env Rscript
library(tidyverse)
library(phyloseq)
library(ggsignif)

# Load data
load("phylo_soil_rare.RData")

######## Comparison of two means with t-test (Parametric) ##########
# Let's do very simple plot with t-test
plot_richness(phylo_soil_rare, x = "Ecozone", measures="Observed") + geom_boxplot()

# Need to extract information
alphadiv <- estimate_richness(phylo_soil_rare)
samp_dat <- sample_data(phylo_soil_rare)
samp_dat_wdiv <- data.frame(samp_dat, alphadiv)

View(samp_dat_wdiv)

# t.test()
t.test(Shannon ~ Ecozone, data=samp_dat_wdiv) #I like this one more

# Note: you can set variances to be equal for a "classic" t-test
t.test(Shannon ~ subject, data=samp_dat_wdiv, var.equal=TRUE) #But we won't use this


######## Comparison of two means with Wilcoxon (Non-parametric) ##########
# Let's say we want to know whether richness (Observed) differs between subjects
# rich <- estimate_richness(mpt, measures = c("Shannon", "Chao1", "Observed"))
# samp_dat <- data.frame(samp_dat, rich)

# Microbial count data is generally NON-NORMAL
# In fact, it is even more complex because microbial data is usually in RELATIVE ABUNDANCE
allCounts <- as.vector(otu_table(phylo_soil_rare))
allCounts <- allCounts[allCounts>0]
hist(allCounts)
hist(log(allCounts))

# Visualize the relationship between richness and subject
samp_dat_wdiv %>%
  filter(!is.na(Observed)) %>%
  ggplot(aes(x=Ecozone, y=Observed))+
  geom_point() 

# Attempt a t-test assuming normal distribution
t.test(Observed ~ Ecozone, data=samp_dat_wdiv)

# Non-parametric distribution means it does not follow a distribution-- it may
# be irregular, skewed, or both
Colour <- c("#8338ec", "#000000", "#FFEB3B", "#1f78b4",
            "#b2df8a", "#e31a1c")

ggplot(samp_dat_wdiv, aes(fill = LTSP.Treatment, colour = Ecozone)) +
         geom_histogram(aes(x=Observed), bins=20, linewidth = 1) +
  scale_color_manual(values = c(Colour[4], Colour[2]))+
  labs(x = "Features Observed in Sample", y = "Frequency")

## There are two ways of dealing with non-parametric data distributions:

# (1) Transform your data (usually with a log function)
ggplot(samp_dat_wdiv, aes(fill = Ecozone)) +
  geom_histogram(aes(x=log(Observed)), bins=20)

# Log transform then perform a t-test
t.test(log(Observed) ~ Ecozone, data=samp_dat_wdiv)


# Let's see what transformed data looks like:
yesno_sampdat %>%
  filter(!is.na(Observed)) %>%
  ggplot(aes(x=Ecozone, y=log(Observed)))+
  geom_boxplot() +
  geom_jitter() #spread out scatterplot

#Using Wilcoxon on LTSP Treatment by converting to OM Removal Yes or No

yesno_sampdat <- mutate(samp_dat_wdiv, 
                        OMR = if_else(LTSP.Treatment == "REF", "No", "Yes"))

# (2) Or you can use a non-parametric test, which usually uses "ranks" or random sampling
# instead of pre-defined distributions
wilcox.test(Observed ~ Ecozone, data=yesno_sampdat, exact = FALSE)
# Notice that transformation does not affect significance
wilcox.test(log(Observed) ~ Ecozone, data=samp_dat_wdiv, exact = FALSE)

#### Comparison of >2 means with ANOVA (Parametric) ####
# Let's create a Shannon plot for the body sites
ggplot(samp_dat_wdiv) +
  geom_boxplot(aes(x=`body.site`, y=Shannon)) + geom_point(aes(x=`body.site`, y=Shannon))

# Set up our linear model 
lm_shannon_vs_site <- lm(Shannon ~ `body.site`, dat=samp_dat_wdiv)
# Calculate AOV
anova_shannon_vs_site <- aov(lm_shannon_vs_site)
# Summarize to determine if there are significant differences
summary(anova_shannon_vs_site)
# Determine which groups are significant
tukey_sum <- TukeyHSD(anova_shannon_vs_site)

# Mapping the significance to a ggplot
shan_bodysite <- ggplot(samp_dat_wdiv, aes(x=`body.site`, y=Shannon)) +
  geom_boxplot() +
  geom_signif(comparisons = list(c("tongue","right palm"), c("tongue", "left palm")),
              y_position = c(4.2, 4.5),
              annotations = c("***","****"))

#### Comparison of >2 means with Kruskall-Wallis (Non-parametric) ####
ggplot(samp_dat_wdiv, aes(x=`Site`, y=Observed)) +
  geom_boxplot() +
  geom_point()

# Or you run a kruskal-wallis test to show that you get greater significance
kruskal_obs <- kruskal.test(Observed ~ `Site`, data = samp_dat_wdiv)

# log the data and run ANOVA again to see that you also get better significance
# than just the ANOVA without the transformation
lm_ob_vs_site_log <- lm(log(Observed) ~ `Compaction.Treatment`, data=samp_dat_wdiv)
anova_ob_vs_site_log <- aov(lm_ob_vs_site_log)
summary(anova_ob_vs_site_log)
TukeyHSD(anova_ob_vs_site_log)

# let's graph the observed features to the body site and annotate significance
# based on the logged ANOVA above
# Mapping the significance to a ggplot
ggplot(samp_dat_wdiv, aes(x=`body.site`, y=Observed)) +
  geom_boxplot() +
  geom_signif(comparisons = list(c("tongue","right palm"), c("tongue", "left palm"), c("tongue","gut")),
              y_position = c(148, 160, 172),
              annotations = c("0.001","0.0001","0.002"))

#### Correlations between two variables with Pearson (Parametric) ####
# prepare the data by removing NAs and zeros
samp_dat_wdiv <- samp_dat_wdiv %>%
  filter(!is.na(PD), !is.na(`CN.Ratio`), !is.na(Observed), PD!=0)

# Let's look at Shannon versus days since the experiment started
ggplot(samp_dat_wdiv,aes(x=CN.Ratio, y=PD)) +
  geom_point() +
  geom_smooth(method = "lm")

# Default option for cor.test is a parametric test called "Pearson"
# cor.test( Xvar, Yvar )
cor.test(samp_dat_wdiv$days.since.experiment.start, samp_dat_wdiv$Shannon)


#### Correlations between two variables with Spearman (Non-parametric) ####
# Let's look at how days of the month affects alpha diversity
ggplot(samp_dat_wdiv,aes(x=Moisture.Content, y=Observed)) +
  geom_point()+
  geom_smooth(method="lm")

# Again, we have two options: transform, or use a non-parametric test
ggplot(samp_dat_wdiv, aes(x=log(Moisture.Content), y=Observed)) +
  geom_point() +
  geom_smooth(method = "lm")
cor.test(log(samp_dat_wdiv$`Moisture.Content`), samp_dat_wdiv$Observed)

# Or, we can simply use a non-parametric test:
# cor.test( Xvar, Yvar, method = 'spearman')
cor.test(samp_dat_wdiv$Moisture.Content, 
         samp_dat_wdiv$Observed, method = "spearman", exact = FALSE)


### PERMANOVA (Permutational ANOVA) ####
# non-parametric version of ANOVA
# Takes a distance matrix, which can be calculated with any kind of metric you want
# e.g. Bray, Jaccard, Unifrac
# Need the package, "vegan"
# Use phyloseq to calculate weighted Unifrac distance matrix
?UniFrac
dm_unifrac <- UniFrac(phylo_soil_rare, weighted=TRUE)

dm_bray <- distance(phylo_soil_rare, method = "bray")

# plot the above as an ordination to a PCoA plot
ord.unifrac <- ordinate(phylo_soil_rare, method="PCoA", distance="wunifrac")

ord.bc <- ordinate(phylo_soil_rare, method="PCoA", distance="bray")
# Also use other metrics: for example, the vegan package includes bray and jaccard

plot_ordination(phylo_soil_rare, ord.unifrac, color="Ecozone", shape = "Site")

# run the permanova on the above matrix for weighted unifrac
?adonis2
adonis2(dm_unifrac ~ Ecozone + Site + LTSP.Treatment + Ecozone:LTSP.Treatment + Site:LTSP.Treatment,
        data = samp_dat_wdiv,
        by = "terms")

adonis2(dm_bray ~ Ecozone + Site + LTSP.Treatment + Ecozone:LTSP.Treatment + Site:LTSP.Treatment,
        data = samp_dat_wdiv,
        by = "terms")

# re-plot the above PCoA with ellipses to show a significant difference 
# between body sites using ggplot2
plot_ordination(phylo_soil_rare, ord.unifrac, color = "Ecozone", shape = "Site") +
  stat_ellipse(type = "t")

# can also use the ggforce package's geom_mark_ellipse function
plot_ordination(phylo_soil_rare, ord.unifrac, color = "Ecozone", shape="Site") +
  geom_mark_ellipse(aes())



### 2-way ANOVA (not being used) ####

# let's plot the subject and reported antibiotic usage against Shannon
# used facet_grid to customize the order of antibiotic usage so it's not alphabetical
ggplot(subset(samp_dat_wdiv, Ecozone == "SBSBC")) + geom_boxplot(aes(x=Site, y=Shannon)) +
  facet_grid(~factor(`LTSP.Treatment`))

ggplot(subset(samp_dat_wdiv, Ecozone == "IDFBC")) + geom_boxplot(aes(x=Site, y=Shannon)) +
  facet_grid(~factor(`LTSP.Treatment`))

# run the 2-way ANOVA
ml_anti_sub <- lm(Shannon ~ `Site`*`LTSP.Treatment`, data=samp_dat_wdiv)
summary(aov(ml_anti_sub))
TukeyHSD(aov(ml_anti_sub))