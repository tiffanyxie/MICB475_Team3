
# R Script to Compare metadata categories across OM treatments
# Last modified: March 11, 2026

#### Load Packages ####
library(tidyverse)
library(readxl)
library(ggplot2)
library(ggtext)
library(FSA)
library(patchwork)

#### Load metadata ####
soil_meta_raw<-read_delim("soil_metadata.txt")
#Replace spaces in column names
colnames(soil_meta_raw)<-str_replace_all(colnames(soil_meta_raw)," ", "_") 
colnames(soil_meta_raw)<-str_replace_all(colnames(soil_meta_raw),"\\(|\\)","")

#Clean up column names and convert horizion, ecozone, and treatments into factors
soil_meta<-soil_meta_raw %>%
  rename(Sample = `#SampleID`) %>%
  mutate(Horizon = factor(Horizon, levels =c("A horizon","O horizon"))) %>%
  mutate(Ecozone = factor(Ecozone,levels=c("IDFBC","SBSBC",
                                           "PPCA","BSON","JPON","LPTX"))) %>%
  mutate(LTSP_Treatment = factor(LTSP_Treatment,levels=c("REF","OM1","OM2","OM3"))) %>%
  mutate(Compaction_Treatment = factor(Compaction_Treatment, levels = c("REF","C0","C1","C2")))

#### Filter to BC, O horizon, no herbicide samples ####
bc_soil<-soil_meta %>%
  filter(Ecozone %in% c("IDFBC","SBSBC")) %>%
  filter(Horizon == "O horizon") %>%
  filter(Herbicide_Use == 0)

idfbc<-bc_soil %>% filter(Ecozone == "IDFBC")
sbsbc<-bc_soil %>% filter(Ecozone == "SBSBC")

#### Boxplots comparing moisture between OM treatment ####


plot_moisture<-function(df,plt_title){
  
  #Filter to samples with moisture content level
  df_filt <- df %>%
    filter(!is.na(Moisture_Content))
  
  #Get number of samples in each TLSP treatment
  sample_num<-df_filt %>%
    group_by(LTSP_Treatment) %>%
    summarise(Samples = n()) %>%
    mutate(Label = paste0(LTSP_Treatment,"\n (n = ",Samples,")"))
  
  #Create boxplot
  box_plot<- df_filt %>%
  ggplot(aes(x=LTSP_Treatment,y=Moisture_Content)) +
    scale_x_discrete(labels=sample_num$Label) +
    xlab("LTSP Treatment") +
    ylab("Moisture Content (%)") +
    ggtitle(plt_title) +
    theme_classic() +
    geom_boxplot()
  
  return(box_plot)
}

bc_moisture<-plot_moisture(bc_soil,"BC")
idfbc_moisture<-plot_moisture(idfbc,"IDFBC")
sbs_moisture<-plot_moisture(sbsbc,"SBSBC")
bc_moisture + idfbc_moisture + sbs_moisture

#### Statistical tests comparing moisture between OM treatments ####

#Overall BC soil
kruskal.test(Moisture_Content ~ LTSP_Treatment,bc_soil) 
dunnTest(Moisture_Content ~ LTSP_Treatment,bc_soil) 


#IDFBC 
kruskal.test(Moisture_Content ~ LTSP_Treatment,idfbc) 
dunnTest(Moisture_Content ~ LTSP_Treatment,idfbc) 

#SBSBC 
kruskal.test(Moisture_Content ~ LTSP_Treatment,sbsbc) 
dunnTest(Moisture_Content ~ LTSP_Treatment,sbsbc) 




#### Boxplots comparing total carbon between OM treatment ####
plot_carbon<-function(df,plt_title){
  
  #Filter to samples with moisture content level
  df_filt <-df %>%
    filter(!is.na(Total_Carbon))
  
  #Get number of samples in each TLSP treatment
  sample_num<-df_filt %>%
    group_by(LTSP_Treatment) %>%
    summarise(Samples = n()) %>%
    mutate(Label = paste0(LTSP_Treatment,"\n (n = ",Samples,")"))
  
  #Create boxplot
  box_plot<- df_filt %>%
    ggplot(aes(x=LTSP_Treatment,y=Total_Carbon)) +
    scale_x_discrete(labels=sample_num$Label) +
    xlab("LTSP Treatment") +
    ylab("Total Carbon") +
    ggtitle(plt_title) +
    theme_classic() +
    geom_boxplot()
  
  return(box_plot)
}

bc_carbon<-plot_carbon(bc_soil,"BC")
idfbc_carbon<-plot_carbon(idfbc,"IDFBC")
sbsbc_carbon<-plot_carbon(sbsbc,"SBSBC")

bc_carbon + idfbc_carbon + sbsbc_carbon

#### Statistical tests comparing carbon between OM treatments ####

#Overall BC soil
kruskal.test(Total_Carbon ~ LTSP_Treatment,bc_soil)
dunnTest(Total_Carbon ~ LTSP_Treatment,bc_soil)


#IDFBC 
kruskal.test(Total_Carbon ~ LTSP_Treatment,idfbc) 
dunnTest(Total_Carbon ~ LTSP_Treatment,idfbc)

#SBSBC 
kruskal.test(Total_Carbon ~ LTSP_Treatment,sbsbc) 
dunnTest(Total_Carbon ~ LTSP_Treatment,sbsbc) 


#### Boxplots comparing total nitrogen between OM treatment ####
plot_nitrogen<-function(df,plt_title){
  
  #Filter to samples with moisture content level
  df_filt <-df %>%
    filter(!is.na(Total_Nitrogen))
  
  #Get number of samples in each TLSP treatment
  sample_num<-df_filt %>%
    group_by(LTSP_Treatment) %>%
    summarise(Samples = n()) %>%
    mutate(Label = paste0(LTSP_Treatment,"\n (n = ",Samples,")"))
  
  #Create boxplot
  box_plot<- df_filt %>%
    ggplot(aes(x=LTSP_Treatment,y=Total_Nitrogen)) +
    scale_x_discrete(labels=sample_num$Label) +
    xlab("LTSP Treatment") +
    ylab("Total Nitrogen") +
    ggtitle(plt_title) +
    theme_classic() +
    geom_boxplot()
  
  return(box_plot)
}

bc_nitrogen<-plot_nitrogen(bc_soil,"BC")
idfbc_nitrogen<-plot_nitrogen(idfbc,"IDFBC")
sbsbc_nitrogen<-plot_nitrogen(sbsbc,"SBSBC")

bc_nitrogen + idfbc_nitrogen + sbsbc_nitrogen

#### Statistical tests comparing nitrogen between OM treatments ####

#Overall BC soil
kruskal.test(Total_Nitrogen ~ LTSP_Treatment,bc_soil) 
dunnTest(Total_Nitrogen ~ LTSP_Treatment,bc_soil) 


#IDFBC 
kruskal.test(Total_Nitrogen ~ LTSP_Treatment,idfbc) 
dunnTest(Total_Nitrogen ~ LTSP_Treatment,idfbc) 

#SBSBC 
kruskal.test(Total_Nitrogen ~ LTSP_Treatment,sbsbc) 
dunnTest(Total_Nitrogen ~ LTSP_Treatment,sbsbc) 



#### Boxplots comparing CN Ratio between OM treatment ####
plot_cn_ratio<-function(df,plt_title){
  
  #Filter to samples with moisture content level
  df_filt <-df %>%
    filter(!is.na(CN_Ratio))
  
  #Get number of samples in each TLSP treatment
  sample_num<-df_filt %>%
    group_by(LTSP_Treatment) %>%
    summarise(Samples = n()) %>%
    mutate(Label = paste0(LTSP_Treatment,"\n (n = ",Samples,")"))
  
  #Create boxplot
  box_plot<- df_filt %>%
    ggplot(aes(x=LTSP_Treatment,y=CN_Ratio)) +
    scale_x_discrete(labels=sample_num$Label) +
    xlab("LTSP Treatment") +
    ylab("C/N Ratio") +
    ggtitle(plt_title) +
    theme_classic() +
    geom_boxplot()
  
  return(box_plot)
}

bc_cn_ratio<-plot_cn_ratio(bc_soil,"BC")
idfbc_cn_ratio<-plot_cn_ratio(idfbc,"IDFBC")
sbsbc_cn_ratio<-plot_cn_ratio(sbsbc,"SBSBC")

bc_cn_ratio + idfbc_cn_ratio + sbsbc_cn_ratio

#### Statistical tests comparing CN Ratio between OM treatments ####

#Overall BC soil
kruskal.test(CN_Ratio ~ LTSP_Treatment,bc_soil) 
dunnTest(CN_Ratio ~ LTSP_Treatment,bc_soil) 


#IDFBC 
kruskal.test(CN_Ratio ~ LTSP_Treatment,idfbc) 
dunnTest(CN_Ratio ~ LTSP_Treatment,idfbc) 

#SBSBC 
kruskal.test(CN_Ratio ~ LTSP_Treatment,sbsbc) 
dunnTest(CN_Ratio ~ LTSP_Treatment,sbsbc) 







#### Boxplots comparing pH  between OM treatment ####
plot_ph<-function(df,plt_title){
  
  #Filter to samples with moisture content level
  df_filt <-df %>%
    filter(!is.na(pH)) %>%
    filter(pH != 0)
  
  #Get number of samples in each TLSP treatment
  sample_num<-df_filt %>%
    group_by(LTSP_Treatment) %>%
    summarise(Samples = n()) %>%
    mutate(Label = paste0(LTSP_Treatment,"\n (n = ",Samples,")"))
  
  #Create boxplot
  box_plot<- df_filt %>%
    ggplot(aes(x=LTSP_Treatment,y=pH)) +
    scale_x_discrete(labels=sample_num$Label) +
    xlab("LTSP Treatment") +
    ylab("pH") +
    ggtitle(plt_title) +
    theme_classic() +
    geom_boxplot()
  
  return(box_plot)
}

bc_ph<-plot_ph(bc_soil,"BC")
idfbc_ph<-plot_ph(idfbc,"IDFBC")
sbsbc_ph<-plot_ph(sbsbc,"SBSBC")

bc_ph + idfbc_ph + sbsbc_ph

#### Statistical tests comparing CN Ratio between OM treatments ####

#Overall BC soil
kruskal.test(pH ~ LTSP_Treatment,bc_soil %>% filter(pH != 0)) 
dunnTest(pH ~ LTSP_Treatment,bc_soil %>% filter(pH != 0)) 


#IDFBC 
kruskal.test(pH ~ LTSP_Treatment,idfbc %>% filter(pH != 0)) 
dunnTest(pH ~ LTSP_Treatment,idfbc %>% filter(pH != 0)) 

#SBSBC 
kruskal.test(pH ~ LTSP_Treatment,sbsbc %>% filter(pH != 0)) 
dunnTest(pH ~ LTSP_Treatment,sbsbc %>% filter(pH != 0)) 


