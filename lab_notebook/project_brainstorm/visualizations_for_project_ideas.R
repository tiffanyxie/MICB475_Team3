# Visualizing metadata metrics to brainstorm ideas for project with soil data
# Tiffany Xie
# Last edited: Feb. 4, 2026

#### Load libraries ####
library(readxl)
library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyverse)

#### Load soil metadata and format ####
soil<-read_excel("soil_metadata.xlsx")
soil$`CN Ratio`<-as.numeric(soil$`CN Ratio`)
soil$`Moisture Content`<-as.numeric(soil$`Moisture Content`)
soil$`Elevation`<-as.numeric(soil$`Elevation`)
soil <- soil %>%
  mutate(Ecozone = factor(Ecozone,levels=c("IDFBC","SBSBC",
                                           "PPCA","BSON","JPON","LPTX")))

#### Idea 1: Comparing LTSP Treatment ####

#Filter out C1 and C2, A horizon sampling depth, herbicide treatment
soil_filt_treatment<-soil %>%
  filter(!`Compaction Treatment` %in% c("C1","C2")) %>%
  filter(Horizon != "A horizon") %>%
  filter(`Herbicide Use`==0)

# Plot # of samples per treatment and ecozone
sample_plot_treatment<-soil_filt_treatment %>%
  ggplot(aes(x=`LTSP Treatment`)) +
  geom_histogram(stat="count") +
  facet_wrap(~`Ecozone`,ncol =3) +
  scale_y_continuous(lim=c(0,13),breaks=seq(0,14,by=2),expand = 0) +
  theme_classic() +
  ylab("# of Samples")
sample_plot_treatment

#Plot elevation of samples in treatment and ecozone
elevation_plot_treatment<-soil_filt_treatment %>%
  ggplot(aes(x=`LTSP Treatment`,y=Elevation)) +
  geom_boxplot() +
  geom_point() +
  facet_wrap(~`Ecozone`,ncol =3) +
  scale_y_continuous(breaks=seq(0,1200,by=400),expand = 0) +
  theme_classic() +
  ylab("Elevation (m)")
elevation_plot_treatment

#Plot moisture content of samples in treatment and ecozone
moisture_plot_treatment<-soil_filt_treatment %>%
  ggplot(aes(x=`LTSP Treatment`,y=`Moisture Content`)) +
  geom_boxplot() +
  geom_point() +
  facet_wrap(~`Ecozone`,ncol =3) +
  theme_classic() +
  scale_y_continuous(expand = 0) +
  ylab("Moisture Content (%)")
moisture_plot_treatment

#Plot CN ratio of samples in treatment and ecozone
cn_plot_treatment<-soil_filt_treatment %>%
  mutate(Ecozone = factor(Ecozone,levels=c("IDFBC","SBSBC",
                                           "PPCA","BSON","JPON","LPTX"))) %>%
  filter(!is.na(`CN Ratio`)) %>%
  ggplot(aes(x=`LTSP Treatment`,y=`CN Ratio`)) +
  geom_boxplot() +
  geom_point() +
  facet_wrap(~`Ecozone`,ncol =3) +
  theme_classic() +
  scale_y_continuous(lim=c(0,52),expand = 0) +
  ylab("CN Ratio")
cn_plot_treatment

# Combine all metadata plots for project idea 1
combined_treatment <- sample_plot_treatment + elevation_plot_treatment + moisture_plot_treatment + cn_plot_treatment
combined_treatment
ggsave("metadata_compare_ltsp_treatment.png",width = 10, height = 8)

#### Idea 2: Comparing Soil Depth ####

#Filter out C1 and C2, OM3 treatment, herbicide treatment
soil_filt_depth<-soil %>%
  filter(!`Compaction Treatment` %in% c("C1","C2")) %>%
  filter(`LTSP Treatment` !="OM3") %>%
  filter(`Herbicide Use`==0)


# Plot # of samples per horizon and ecozone
sample_plot_depth<-soil_filt_depth %>%
  ggplot(aes(x=`Horizon`)) +
  geom_histogram(stat="count") +
  facet_wrap(~`Ecozone`,ncol =3) +
  scale_y_continuous(expand = 0) +
  theme_classic() +
  ylab("# of Samples")
sample_plot_depth

#Plot elevation of samples in horizon and ecozone
elevation_plot_depth<-soil_filt_depth %>%
  ggplot(aes(x=`Horizon`,y=Elevation)) +
  geom_boxplot() +
  geom_point() +
  facet_wrap(~`Ecozone`,ncol =3) +
  scale_y_continuous(expand = 0,lim=c(0,1300)) +
  theme_classic() +
  ylab("Elevation (m)")
elevation_plot_depth

#Plot moisture of samples in horizon and ecozone
moisture_plot_depth<-soil_filt_depth %>%
  ggplot(aes(x=`Horizon`,y=`Moisture Content`)) +
  geom_boxplot() +
  geom_point() +
  facet_wrap(~`Ecozone`,ncol =3) +
  scale_y_continuous(expand=0,lim=c(0, 90)) +
  theme_classic() +
  ylab("Moisture Content (%)")
moisture_plot_depth

#Plot CN ratio of samples in horizon and ecozone
cn_plot_depth<-soil_filt_depth %>%
  filter(!is.na(`CN Ratio`)) %>%
  ggplot(aes(x=`Horizon`,y=`CN Ratio`)) +
  geom_boxplot() +
  geom_point() +
  facet_wrap(~`Ecozone`,ncol =3) +
  scale_y_continuous(expand = 0,lim=c(0,55)) +
  theme_classic() +
  ylab("CN Ratio")

cn_plot_depth

# Combine all metadata plots for project idea 2
combined_depth <- sample_plot_depth + elevation_plot_depth + 
  moisture_plot_depth + cn_plot_depth
combined_depth
ggsave("metadata_compare_horizon.png",width = 10, height = 8)