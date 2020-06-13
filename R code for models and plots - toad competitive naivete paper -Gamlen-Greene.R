################## R code for by Roseanna Gamlen-Greene 05.06.2020 
#roseanna.gamlen.greene@gmail.com

#for UNPUBLISHED paper by Gamlen-Greene and Richardson (CURRENTLY UNDER REVIEW IN ECOLOGY): "Co-occurrence history matters: island toads are competitively naïve to introduced frogs despite mainland co-existence"


#R version 3.6.2 (2019-12-12) -- "Dark and Stormy Night"
#Copyright (C) 2019 The R Foundation for Statistical Computing
#Platform: x86_64-w64-mingw32/x64 (64-bit)

setwd()
rm(list=ls())

######################
######### PACKAGES #########
#######################

library(ggplot2)
library(ggthemes)
library(reshape2)
library(gridExtra)
library(plyr)
library(dplyr)
library(visreg)
library(lme4)
library(nlme)
library(MuMIn)
library(lmerTest)
library(sjPlot)
library(sjmisc)
library(effects)
library(car)
library(optimx)
library(arules)
library(rstanarm)
library(tidyr)
library(dataRetrieval)
library(cowplot)
library(boot)
library(rstanarm)
library(utf8)
library(emmeans)
library(tidyr)
library(naniar)
library(phia)
library(multcomp)
library(data.table)
library(tidyverse)
library(RCurl)
library(afex)
library(tibble)
library(ggstatsplot)

#################################
############### 1. TOAD GROWTH RATE ##########
#######################################
rm(list=ls())

#import data
Rawdata_02_03_20<-read.csv(file="Rawdata_02_03_20.csv",header=T,row.names=NULL,sep=",")


#order levels of factor
Rawdata_02_03_20$Competition_level<- factor(Rawdata_02_03_20$Competition_level, levels = c("Low_control","FrogML","FrogHG","High_control"))
Rawdata_02_03_20<- droplevels(Rawdata_02_03_20)

#days as factor
Rawdata_02_03_20$Daysfactor<-as.factor(datafile$Days)

############ calculate toad growth rate
Rawdata_02_03_20$GR_finalmin_intial_degreedays_mgDD<-((Rawdata_02_03_20$Av_juv_biomass-Rawdata_02_03_20$X25.5.18_experiment_start_tapole_av_individual_weight)/Rawdata_02_03_20$SumDD)*1000

#### percetnage tadpole additions due to diving beetle predation
Rawdata_02_03_20$toadTadpoles_added_as_top.up_june7thand8th_percentagetotal<-(Rawdata_02_03_20$Tadpoles_added_as_top.up_june7thand8th/Rawdata_02_03_20$Total_number_individuals_start)*100

#########################
############### SCALING ##########
######################

#starting weight
Xtemp <- Rawdata_02_03_20$X25.5.18_experiment_start_tapole_av_individual_weight 
Xtemp

Xscaled <- (Xtemp - mean(Xtemp))/sd(Xtemp)
Rawdata_02_03_20$scaled_X25.5.18_experiment_start_tapole_av_individual_weight  <- Xscaled

#######################################
############# Linear Mixed Effects Model - TOAD GROWTH RATE ############
#########################################

model_GR_toads_FinalmInitial_degreeDays.Frank.3<-lmer(
  GR_finalmin_intial_degreedays_mgDD~
    Region*Competition_level+ 
    toadTadpoles_added_as_top.up_june7thand8th_percentagetotal+ 
    scaled_X25.5.18_experiment_start_tapole_av_individual_weight+
    Days+
    (1|Population)+ 
    (1+Days|Tank), #random effects
  data=datafile,na.action=na.exclude,REML=TRUE,control=lmerControl(optimizer="optimx", optCtrl=list(method='nlminb'),check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))

anova(model_GR_toads_FinalmInitial_degreeDays.Frank.3)
AIC(model_GR_toads_FinalmInitial_degreeDays.Frank.3) #-1887.154

#model_GR_toads_FinalmInitial_degreeDays.Frank.3
#capture.output(anova(model_GR_toads_FinalmInitial_degreeDays.Frank.3),file="anova_model_GR_toads_FinalmInitial_degreeDays.Frank.3_REMLFALSE.txt")
#capture.output(model_GR_toads_FinalmInitial_degreeDays.Frank.3,file="model_GR_toads_FinalmInitial_degreeDays.Frank.3_REMLisFALSE.txt")

#use REML=TRUE for the custom contrasts and for reporting F stats
#capture.output(anova(model_GR_toads_FinalmInitial_degreeDays.Frank.3),file="anova_model_GR_toads_FinalmInitial_degreeDays.Frank.3_REMLisTRUE.txt")
#capture.output(model_GR_toads_FinalmInitial_degreeDays.Frank.3,file="model_GR_toads_FinalmInitial_degreeDays.Frank.3_REMLisTRUE.txt")


#######################################################################
############################################################################
################## CUSTOM CONTRASTS TOAD GROWTH RATE ####################
############################################################################
############################################################################


#REML=TRUE
#model_GR_toads_FinalmInitial_degreeDays.Frank.3
emm3_model_RGRDD_finalmin_intial_sumdegreedays=emmeans(model_GR_toads_FinalmInitial_degreeDays.Frank.3,specs=~ Competition_level*Region)


#turn to txt
#capture.output(emm3_model_RGRDD_finalmin_intial_sumdegreedays,file="CIs_eemeans_model_GR_toads_FinalmInitial_degreeDays.Frank.3.txt")


#make into dataframe
emm3_model_RGRDD_finalmin_intial_sumdegreedaysdf<-as.data.frame(emm3_model_RGRDD_finalmin_intial_sumdegreedays)

#make new column
emm3_model_RGRDD_finalmin_intial_sumdegreedaysdf<-cbind(Model = "model_GR_toads_FinalmInitial_degreeDays.Frank.3", emm3_model_RGRDD_finalmin_intial_sumdegreedaysdf)

#emmean values
#region HG, comp level low control - single out estimate for each level of region/competition to insect into contrast table later to calculate percent
HGLow_controlHG_model_RGRDD<-emm3_model_RGRDD_finalmin_intial_sumdegreedaysdf[1,4]
HGFrogMLHG_model_RGRDD<-emm3_model_RGRDD_finalmin_intial_sumdegreedaysdf[2,4]
HGFrogHGHG_model_RGRDD<-emm3_model_RGRDD_finalmin_intial_sumdegreedaysdf[3,4]
HGHigh_controlHG_model_RGRDD<-emm3_model_RGRDD_finalmin_intial_sumdegreedaysdf[4,4]
MLLow_controlML_model_RGRDD<-emm3_model_RGRDD_finalmin_intial_sumdegreedaysdf[5,4] 
MLFrogMLML_model_RGRDD<-emm3_model_RGRDD_finalmin_intial_sumdegreedaysdf[6,4]
MLFrogHGML_model_RGRDD<-emm3_model_RGRDD_finalmin_intial_sumdegreedaysdf[7,4]
MLHigh_controlML_model_RGRDD<-emm3_model_RGRDD_finalmin_intial_sumdegreedaysdf[8,4]
HGFrog_overallHG_model_RGRDD<-(emm3_model_RGRDD_finalmin_intial_sumdegreedaysdf[2,4]+emm3_model_RGRDD_finalmin_intial_sumdegreedaysdf[3,4])/2
MLFrog_overallML_model_RGRDD<-(emm3_model_RGRDD_finalmin_intial_sumdegreedaysdf[6,4]+emm3_model_RGRDD_finalmin_intial_sumdegreedaysdf[7,4])/2


#create custom contrasts
Low_controlHG<-rep(c(1,0), times = c(1,7))
FrogMLHG<-rep(c(0,1,0), times = c(1,1,6))
FrogHGHG<-rep(c(0,1,0), times = c(2,1,5))
High_controlHG<-rep(c(0,1,0), times = c(3,1,4))
Low_controlML<-rep(c(0,1,0), times = c(4,1,3))
FrogMLML<-rep(c(0,1,0), times = c(5,1,2))
FrogHGML<-rep(c(0,1,0), times = c(6,1,1))
High_controlML<-rep(c(0,1), times = c(7,1))

#create joined frog response variable
Frog_overallHG=(FrogMLHG+FrogHGHG)/2
Frog_overallML=(FrogMLML+FrogHGML)/2

#bonferroni
(contrast_model_RGRDD_finalmin_intial_sumdegreedays_bonferroni2<-contrast(emm3_model_RGRDD_finalmin_intial_sumdegreedays,adjust="bonferroni",method=list(
  "FrogHGHG - FrogMLHG" = FrogHGHG - FrogMLHG, #not different
  "FrogHGML - FrogMLML" = FrogHGML - FrogMLML, #not different
  "Low_controlHG - Frog_overallHG" = Low_controlHG - Frog_overallHG,
  "Low_controlML - Frog_overallML"=Low_controlML - Frog_overallML, 
  "Frog_overallHG - High_controlHG"=Frog_overallHG - High_controlHG,   
  "Frog_overallML - High_controlML"=Frog_overallML - High_controlML,
  "Low_controlHG - High_controlHG" = Low_controlHG - High_controlHG,
  "Low_controlML - High_controlML"=Low_controlML - High_controlML
)))

#fdr
(contrast_model_RGRDD_finalmin_intial_sumdegreedays_fdr<-contrast(emm3_model_RGRDD_finalmin_intial_sumdegreedays,adjust="fdr",method=list(
  "FrogHGHG - FrogMLHG" = FrogHGHG - FrogMLHG, #not different
  "FrogHGML - FrogMLML" = FrogHGML - FrogMLML, #not different
  "Low_controlHG - Frog_overallHG" = Low_controlHG - Frog_overallHG,
  "Low_controlML - Frog_overallML"=Low_controlML - Frog_overallML, 
  "Frog_overallHG - High_controlHG"=Frog_overallHG - High_controlHG,   
  "Frog_overallML - High_controlML"=Frog_overallML - High_controlML,
  "Low_controlHG - High_controlHG" = Low_controlHG - High_controlHG,
  "Low_controlML - High_controlML"=Low_controlML - High_controlML
)))

#mvt
(contrast_model_RGRDD_finalmin_intial_sumdegreedays_mvt<-contrast(emm3_model_RGRDD_finalmin_intial_sumdegreedays,adjust="mvt",method=list(
  "FrogHGHG - FrogMLHG" = FrogHGHG - FrogMLHG, #not different
  "FrogHGML - FrogMLML" = FrogHGML - FrogMLML, #not different
  "Low_controlHG - Frog_overallHG" = Low_controlHG - Frog_overallHG,
  "Low_controlML - Frog_overallML"=Low_controlML - Frog_overallML, 
  "Frog_overallHG - High_controlHG"=Frog_overallHG - High_controlHG,   
  "Frog_overallML - High_controlML"=Frog_overallML - High_controlML,
  "Low_controlHG - High_controlHG" = Low_controlHG - High_controlHG,
  "Low_controlML - High_controlML"=Low_controlML - High_controlML
)))

#capture
#capture.output(contrast_model_RGRDD_finalmin_intial_sumdegreedays_mvt,file="contrast_model_RGRDD_finalmin_intial_sumdegreedays_mvt_model_GR_toads_FinalmInitial_degreeDays.Frank.3.txt")


#######################
################## PERCENTS ##################
####################

#make into dataframe
(contrast_model_RGRDD_finalmin_intial_sumdegreedaysdf<-as.data.frame(contrast_model_RGRDD_finalmin_intial_sumdegreedays_mvt))


#make new column
contrast_model_RGRDD_finalmin_intial_sumdegreedaysdf<-cbind(Model = "model_model_GR_toads_FinalmInitial_degreeDays.Frank.3", contrast_model_RGRDD_finalmin_intial_sumdegreedaysdf)

#HGLow_controlHG
(contrast_model_RGRDD_finalmin_intial_sumdegreedaysdf$emmean_HGLow_controlHG<-rep(HGLow_controlHG_model_RGRDD,length(contrast_model_RGRDD_finalmin_intial_sumdegreedaysdf$p.value)))
#percent
(contrast_model_RGRDD_finalmin_intial_sumdegreedaysdf$percent_HGLow_controlHG<-(contrast_model_RGRDD_finalmin_intial_sumdegreedaysdf$estimate/contrast_model_RGRDD_finalmin_intial_sumdegreedaysdf$emmean_HGLow_controlHG)*100)
#value
((HGLow_controlHG_model_RGRDD_finalmin_intial_sumdegreedaysdf_value<-contrast_model_RGRDD_finalmin_intial_sumdegreedaysdf[c(3,7),c(1:9)]))



#MLLow_controlML
(contrast_model_RGRDD_finalmin_intial_sumdegreedaysdf$emmean_MLLow_controlML<-rep(MLLow_controlML_model_RGRDD,length(contrast_model_RGRDD_finalmin_intial_sumdegreedaysdf$p.value)))
#percent
(contrast_model_RGRDD_finalmin_intial_sumdegreedaysdf$percent_MLLow_controlML<-(contrast_model_RGRDD_finalmin_intial_sumdegreedaysdf$estimate/contrast_model_RGRDD_finalmin_intial_sumdegreedaysdf$emmean_MLLow_controlML)*100)
#value
(MLLow_controlML_model_RGRDD_finalmin_intial_sumdegreedaysdf_value<-contrast_model_RGRDD_finalmin_intial_sumdegreedaysdf[c(4,8),c(1:7,10:11)])

#HGFrogoverallHG
(contrast_model_RGRDD_finalmin_intial_sumdegreedaysdf$emmean_HGFrog_overallHG<-rep(HGFrog_overallHG_model_RGRDD,length(contrast_model_RGRDD_finalmin_intial_sumdegreedaysdf$p.value)))
#percent
(contrast_model_RGRDD_finalmin_intial_sumdegreedaysdf$percent_HGFrog_overallHG<-(contrast_model_RGRDD_finalmin_intial_sumdegreedaysdf$estimate/contrast_model_RGRDD_finalmin_intial_sumdegreedaysdf$emmean_HGFrog_overallHG)*100)
#value
((HGFrog_overallHG_model_RGRDD_finalmin_intial_sumdegreedaysdf_value<-contrast_model_RGRDD_finalmin_intial_sumdegreedaysdf[c(5),c(1:7,12:13)]))

#MLFrogoverallML
contrast_model_RGRDD_finalmin_intial_sumdegreedaysdf$emmean_MLFrog_overallML<-rep(MLFrog_overallML_model_RGRDD,length(contrast_model_RGRDD_finalmin_intial_sumdegreedaysdf$p.value))
#percent
contrast_model_RGRDD_finalmin_intial_sumdegreedaysdf$percent_MLFrog_overallML<-(contrast_model_RGRDD_finalmin_intial_sumdegreedaysdf$estimate/contrast_model_RGRDD_finalmin_intial_sumdegreedaysdf$emmean_MLFrog_overallML)*100
#value
(MLFrog_overallML_model_RGRDD_finalmin_intial_sumdegreedaysdf_value<-contrast_model_RGRDD_finalmin_intial_sumdegreedaysdf[c(6),c(1:7,14:15)])

#write.csv(contrast_model_RGRDD_finalmin_intial_sumdegreedaysdf,"contrast_model_GR_toads_FinalmInitial_degreeDays.Frank.csv")

#bind together
HGLow_controlHG_model_RGRDD_finalmin_intial_sumdegreedaysdf_value
MLLow_controlML_model_RGRDD_finalmin_intial_sumdegreedaysdf_value
HGFrog_overallHG_model_RGRDD_finalmin_intial_sumdegreedaysdf_value
MLFrog_overallML_model_RGRDD_finalmin_intial_sumdegreedaysdf_value

(contrast_model_RGRDD_finalmin_intial_sumdegreedaysdf_selected<-bind_rows(
  HGLow_controlHG_model_RGRDD_finalmin_intial_sumdegreedaysdf_value,
  MLLow_controlML_model_RGRDD_finalmin_intial_sumdegreedaysdf_value,
  HGFrog_overallHG_model_RGRDD_finalmin_intial_sumdegreedaysdf_value,
  MLFrog_overallML_model_RGRDD_finalmin_intial_sumdegreedaysdf_value))



#   model_GR_toads_FinalmInitial_degreeDays.Frank.3
#write.csv(contrast_model_RGRDD_finalmin_intial_sumdegreedaysdf_selected,"contrast_selected_mvt_model_GR_toads_FinalmInitial_degreeDays.Frank.3_REMLTRUE.csv")



#############
#### CONFIDENCE INTERVALS ####
##############
(contrast_emm3_model_RGRDD_finalmin_intial_sumdegreedays_CI<-contrast(emm3_model_RGRDD_finalmin_intial_sumdegreedays,method=list(
  "Low_controlHG - Frog_overallHG" = Low_controlHG - Frog_overallHG,
  "Low_controlML - Frog_overallML"=Low_controlML - Frog_overallML,  
  "Frog_overallHG - High_controlHG"=Frog_overallHG - High_controlHG,   
  "Frog_overallML - High_controlML"=Frog_overallML - High_controlML,
  "Low_controlHG - High_controlHG" = Low_controlHG - High_controlHG,
  "Low_controlML - High_controlML"=Low_controlML - High_controlML),
  adjust = "bonferroni")%>%
   confint())

#capture.output(contrast_emm3_model_RGRDD_finalmin_intial_sumdegreedays_CI,file="contrast_model_GR_toads_FinalmInitial_degreeDays.4_CI.txt")

#################################
#################### PLOTS ######
################################


labels <- c(ML = "Mainland Toads", HG = "Haida Gwaii Toads")


#GR_finalmin_intial_degreedays_mgDD
ylim1 = boxplot.stats(Rawdata_02_03_20$GR_finalmin_intial_degreedays_mgDD)$stats[c(1, 5)]



(boxplot_GR_finalmin_intial_degreedays_mgDD<- Rawdata_02_03_20%>%
    mutate(Treatment_three = fct_relevel(Treatment_three, 
                                         "Low_control", "Frog","High_control")) %>%
    ggplot(aes(x=Treatment_three, y=GR_finalmin_intial_degreedays_mgDD)) +
    geom_boxplot(lwd=0.6) +
    scale_fill_manual("grey")+
    coord_cartesian(ylim = ylim1*1.12)+
    facet_wrap(~Region,labeller = labeller(Region = labels))+
    labs(x="Competition Treatment")+
    labs(y="Toad growth rate (mg/degree day)")+
    scale_x_discrete(breaks=c("Low_control",  "Frog","High_control"),
                     labels=c("Low intraspecific\n80 toads","Interspecific\n16 frogs\n80 toads", "High intraspecific\n160 toads" ))+
    scale_colour_discrete(name="Toad Region",breaks=c("HG", "ML"),
                          labels=c("Haida Gwaii", "Mainland"))+
    theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()))



##############################################################################
##############################################################################
########## 2. Toad weight at metamorphosis ###################################
##############################################################################
##############################################################################

rm(list=ls())

Rawdata_02_03_20<-read.csv(file="Rawdata_02_03_20.csv",header=T,row.names=NULL,sep=",")



Rawdata_02_03_20$toadTadpoles_added_as_top.up_june7thand8th_percentagetotal<-(Rawdata_02_03_20$Tadpoles_added_as_top.up_june7thand8th/Rawdata_02_03_20$Total_number_individuals_start)*100


Rawdata_02_03_20$Av_juv_biomass_mg<-Rawdata_02_03_20$Av_juv_biomass*1000

#days as factor
Rawdata_02_03_20$Daysfactor<-as.factor(Rawdata_02_03_20$Days)





#order levels of factor
Rawdata_02_03_20$Competition_level<- factor(Rawdata_02_03_20$Competition_level, levels = c("Low_control","FrogML","FrogHG","High_control"))
Rawdata_02_03_20<- droplevels(Rawdata_02_03_20)



#starting weight
Xtemp <- Rawdata_02_03_20$X25.5.18_experiment_start_tapole_av_individual_weight 
Xtemp

Xscaled <- (Xtemp - mean(Xtemp))/sd(Xtemp)
Rawdata_02_03_20$scaled_X25.5.18_experiment_start_tapole_av_individual_weight  <- Xscaled


##############################################################
############# LINEAR MIXED EFFECTS MODEL - TOAD WEIGHT AT METAMORPHOSIS ##########################################################
##########################################

model_avjuvbiomass_toads_Frank.2<-lmer(
  Av_juv_biomass~
    Region*Competition_level+ 
    toadTadpoles_added_as_top.up_june7thand8th_percentagetotal+ 
    scaled_X25.5.18_experiment_start_tapole_av_individual_weight+
    Mean_averageTemp+ #fixed effects
    Days+
    Mean_averageTemp+
    (1|Population)+ 
    (1+Days|Tank), #random effects
  data=Rawdata_02_03_20_rmbigouts,na.action=na.exclude,REML=TRUE,control=lmerControl(optimizer="optimx", optCtrl=list(method='nlminb'),check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))#

anova(model_avjuvbiomass_toads_Frank.2)
AIC(model_avjuvbiomass_toads_Frank.2) #

#capture.output(anova(model_avjuvbiomass_toads_Frank.2),file="anova_model_avjuvbiomass_toads_Frank.2_REMLFALSE.txt")
#capture.output(model_avjuvbiomass_toads_Frank.2,file="model_avjuvbiomass_toads_Frank.2_REMLFALSE.txt")

#capture.output(anova(model_avjuvbiomass_toads_Frank.2),file="anova_model_avjuvbiomass_toads_Frank.2_REMLTRUE.txt")
#capture.output(model_avjuvbiomass_toads_Frank.2,file="model_avjuvbiomass_toads_Frank.2_REMLTRUE.txt")


############################################################################
################################## CUSTOM CONTRASTS ########################
############################################################################


#model_avjuvbiomass_toads_Frank.2 #REML = TRUE
emm3_model_avjuvweight=emmeans(model_avjuvbiomass_toads_Frank.2,specs=~ Competition_level*Region)

#capture
#capture.output(emm3_model_avjuvweight,file="CIS_emmeans_model_avjuvbiomass_toads_Frank.2.txt")

#make into dataframe
emm3_model_avjuvweightdf<-as.data.frame(emm3_model_avjuvweight)

#make new column
emm3_model_avjuvweightdf<-cbind(Model = "model_avjuvbiomass_toads_Frank.2", emm3_model_avjuvweightdf)


#region HG, comp level low control - single out estimate for each level of region/competition to insect into contrast table later to calculate percent
HGLow_controlHG_model_avjuvweightdf<-emm3_model_avjuvweightdf[1,4]
HGFrogMLHG_model_avjuvweightdf<-emm3_model_avjuvweightdf[2,4]
HGFrogHGHG_model_avjuvweightdf<-emm3_model_avjuvweightdf[3,4]
HGHigh_controlHG_model_avjuvweightdf<-emm3_model_avjuvweightdf[4,4]
MLLow_controlML_model_avjuvweightdf<-emm3_model_avjuvweightdf[5,4] 
MLFrogMLML_model_avjuvweightdf<-emm3_model_avjuvweightdf[6,4]
MLFrogHGML_model_avjuvweightdf<-emm3_model_avjuvweightdf[7,4]
MLHigh_controlML_model_avjuvweightdf<-emm3_model_avjuvweightdf[8,4]
HGFrog_overallHG_model_avjuvweightdf<-(emm3_model_avjuvweightdf[2,4]+emm3_model_avjuvweightdf[3,4])/2
MLFrog_overallML_model_avjuvweightdf<-(emm3_model_avjuvweightdf[6,4]+emm3_model_avjuvweightdf[7,4])/2


#create custom contrasts
Low_controlHG<-rep(c(1,0), times = c(1,7))
FrogMLHG<-rep(c(0,1,0), times = c(1,1,6))
FrogHGHG<-rep(c(0,1,0), times = c(2,1,5))
High_controlHG<-rep(c(0,1,0), times = c(3,1,4))
Low_controlML<-rep(c(0,1,0), times = c(4,1,3))
FrogMLML<-rep(c(0,1,0), times = c(5,1,2))
FrogHGML<-rep(c(0,1,0), times = c(6,1,1))
High_controlML<-rep(c(0,1), times = c(7,1))

#create joined frog response variable
Frog_overallHG=(FrogMLHG+FrogHGHG)/2
Frog_overallML=(FrogMLML+FrogHGML)/2


#bonferroni
(contrast_model_avjuvweight_bonferroni2<-contrast(emm3_model_avjuvweight,adjust="bonferroni",method=list(
  "FrogHGHG - FrogMLHG" = FrogHGHG - FrogMLHG, #not different
  "FrogHGML - FrogMLML" = FrogHGML - FrogMLML, #not different
  "Low_controlHG - Frog_overallHG" = Low_controlHG - Frog_overallHG,
  "Low_controlML - Frog_overallML"=Low_controlML - Frog_overallML, 
  "Frog_overallHG - High_controlHG"=Frog_overallHG - High_controlHG,   
  "Frog_overallML - High_controlML"=Frog_overallML - High_controlML,
  "Low_controlHG - High_controlHG" = Low_controlHG - High_controlHG,
  "Low_controlML - High_controlML"=Low_controlML - High_controlML
)))

#fdr
(contrast_model_avjuvweight_fdr<-contrast(emm3_model_avjuvweight,adjust="fdr",method=list(
  "FrogHGHG - FrogMLHG" = FrogHGHG - FrogMLHG, #not different
  "FrogHGML - FrogMLML" = FrogHGML - FrogMLML, #not different
  "Low_controlHG - Frog_overallHG" = Low_controlHG - Frog_overallHG,
  "Low_controlML - Frog_overallML"=Low_controlML - Frog_overallML, 
  "Frog_overallHG - High_controlHG"=Frog_overallHG - High_controlHG,   
  "Frog_overallML - High_controlML"=Frog_overallML - High_controlML,
  "Low_controlHG - High_controlHG" = Low_controlHG - High_controlHG,
  "Low_controlML - High_controlML"=Low_controlML - High_controlML
)))

#mvt
(contrast_model_avjuvweight_mvt<-contrast(emm3_model_avjuvweight,adjust="mvt",method=list(
  "FrogHGHG - FrogMLHG" = FrogHGHG - FrogMLHG, #not different
  "FrogHGML - FrogMLML" = FrogHGML - FrogMLML, #not different
  "Low_controlHG - Frog_overallHG" = Low_controlHG - Frog_overallHG,
  "Low_controlML - Frog_overallML"=Low_controlML - Frog_overallML, 
  "Frog_overallHG - High_controlHG"=Frog_overallHG - High_controlHG,   
  "Frog_overallML - High_controlML"=Frog_overallML - High_controlML,
  "Low_controlHG - High_controlHG" = Low_controlHG - High_controlHG,
  "Low_controlML - High_controlML"=Low_controlML - High_controlML
)))

#capture
#capture.output(contrast_model_avjuvweight_mvt,file="contrast_model_avjuvweight_mvt_model_avjuvbiomass_toads_Frank.2.txt")



#######################################
################## PERCENTS ##################
#######################################

#make into dataframe
(contrast_model_avjuvweightdf<-as.data.frame(contrast_model_avjuvweight_mvt))

#make new column
contrast_model_avjuvweightdf<-cbind(Model = "model_avjuvbiomass_toads_Frank.2", contrast_model_avjuvweightdf)

#HGLow_controlHG
(contrast_model_avjuvweightdf$emmean_HGLow_controlHG<-rep(HGLow_controlHG_model_avjuvweightdf,length(contrast_model_avjuvweightdf$p.value)))
#percent
(contrast_model_avjuvweightdf$percent_HGLow_controlHG<-(contrast_model_avjuvweightdf$estimate/contrast_model_avjuvweightdf$emmean_HGLow_controlHG)*100)
#value
((HGLow_controlHG_model_avjuvweightdf_value<-contrast_model_avjuvweightdf[c(3,7),c(1:9)]))



#MLLow_controlML
(contrast_model_avjuvweightdf$emmean_MLLow_controlML<-rep(MLLow_controlML_model_avjuvweightdf,length(contrast_model_avjuvweightdf$p.value)))
#percent
(contrast_model_avjuvweightdf$percent_MLLow_controlML<-(contrast_model_avjuvweightdf$estimate/contrast_model_avjuvweightdf$emmean_MLLow_controlML)*100)
#value
(MLLow_controlML_model_avjuvweightdf_value<-contrast_model_avjuvweightdf[c(4,8),c(1:7,10:11)])

#HGFrogoverallHG
(contrast_model_avjuvweightdf$emmean_HGFrog_overallHG<-rep(HGFrog_overallHG_model_avjuvweightdf,length(contrast_model_avjuvweightdf$p.value)))
#percent
(contrast_model_avjuvweightdf$percent_HGFrog_overallHG<-(contrast_model_avjuvweightdf$estimate/contrast_model_avjuvweightdf$emmean_HGFrog_overallHG)*100)
#value
((HGFrog_overallHG_model_avjuvweightdf_value<-contrast_model_avjuvweightdf[c(5),c(1:7,12:13)]))

#MLFrogoverallML
contrast_model_avjuvweightdf$emmean_MLFrog_overallML<-rep(MLFrog_overallML_model_avjuvweightdf,length(contrast_model_avjuvweightdf$p.value))
#percent
contrast_model_avjuvweightdf$percent_MLFrog_overallML<-(contrast_model_avjuvweightdf$estimate/contrast_model_avjuvweightdf$emmean_MLFrog_overallML)*100
#value
(MLFrog_overallML_model_avjuvweightdf_value<-contrast_model_avjuvweightdf[c(6),c(1:7,14:15)])

#write.csv(contrast_model_avjuvweightdf,"contrast_model_GR_toads_FinalmInitial_degreeDays.Frank.csv")

#bind together
HGLow_controlHG_model_avjuvweightdf_value
MLLow_controlML_model_avjuvweightdf_value
HGFrog_overallHG_model_avjuvweightdf_value
MLFrog_overallML_model_avjuvweightdf_value

(contrast_model_avjuvweightdf_selected<-bind_rows(
  HGLow_controlHG_model_avjuvweightdf_value,
  MLLow_controlML_model_avjuvweightdf_value,
  HGFrog_overallHG_model_avjuvweightdf_value,
  MLFrog_overallML_model_avjuvweightdf_value))



#   model_avjuvbiomass_toads_Frank.2
#write.csv(contrast_model_avjuvweightdf_selected,"contrast_selected_mvt_model_avjuvbiomass_toads_Frank.2_june8th_REMTRUE.csv")


#################################
#################### PLOTS ######
################################
labels <- c(ML = "Mainland Toads", HG = "Haida Gwaii Toads")


ylim1 = boxplot.stats(Rawdata_02_03_20$Av_juv_biomass_mg)$stats[c(1, 5)]

(boxplot_Competition_level_treat3_JUVWEIGHT_mg_30_03_20_blank <- Rawdata_02_03_20%>%
    mutate(Treatment_three = fct_relevel(Treatment_three, 
                                         "Low_control", "Frog","High_control")) %>%
    ggplot(aes(x=Treatment_three, y=Av_juv_biomass_mg)) +
    geom_boxplot(lwd=0.6) +
    scale_fill_manual("grey")+
    coord_cartesian(ylim = ylim1*1.1)+
    facet_wrap(~Region,labeller = labeller(Region = labels))+
    labs(x="Competition Treatment")+
    labs(y="Average toad weight at metamorphosis (mg)")+
    scale_x_discrete(breaks=c("Low_control",  "Frog","High_control"),
                     labels=c("Low intraspecific\n80 toads","Interspecific\n16 frogs\n80 toads", "High intraspecific\n160 toads" ))+
    scale_colour_discrete(name="Toad Region",breaks=c("HG", "ML"),
                          labels=c("Haida Gwaii", "Mainland"))+
    theme_bw()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()))




###################################################################
############################# 3. MEDIAN TIME TO METAMORPHOSIS ######
#############################################################

rm(list=ls())

Rawdata_02_03_20<-read.csv(file="Rawdata_02_03_20.csv",header=T,row.names=NULL,sep=",")

Rawdata_02_03_20$toadTadpoles_added_as_top.up_june7thand8th_percentagetotal<-(Rawdata_02_03_20$Tadpoles_added_as_top.up_june7thand8th/Rawdata_02_03_20$Total_number_individuals_start)*100



#order levels of factor
Rawdata_02_03_20$Competition_level<- factor(Rawdata_02_03_20, levels = c("Low_control","FrogML","FrogHG","High_control"))
Rawdata_02_03_20<- droplevels(Rawdata_02_03_20)


########## how to calc MedianDays ################################
#make each toad have its own row and then run the code below
Rawdata_02_03_20_indivrows2 <- as.data.frame(lapply(Rawdata_02_03_20, rep, Rawdata_02_03_20$Count_juv_emerged))


#cal median
Rawdata_02_03_20_indivrows2_median2<-setDT(Rawdata_02_03_20_indivrows2)[,list(MeanDays=mean(Days), MaxDays=max(Days), MinDays=min(Days), MedianDays=as.numeric(median(Days)), StdDays=sd(Days)), by=Tank]

#join to large dataset
Rawdata_02_03_20_indivrows2_median2_large<-left_join(subset_dataRGRdd_renamedsomecolumns_2_10_20,Rawdata_02_03_20_indivrows2_median2,by="Tank")


#Write csv
#write.csv(Rawdata_02_03_20_indivrows2_median2_large,"Rawdata_02_03_20_indivrows2_median2_large.csv")



###########################
############### SCALING #######
#################################

#toadTadpoles_added_as_top.up_june7thand8th_percentagetotal

Xtemp <- Rawdata_02_03_20_indivrows2_median2_large$toadTadpoles_added_as_top.up_june7thand8th_percentagetotal
Xtemp

Xscaled <- (Xtemp - mean(Xtemp))/sd(Xtemp)
Rawdata_02_03_20_indivrows2_median2_large$scaled_toadTadpoles_added_as_top.up_june7thand8th_percentagetotal <- Xscaled




######################################
############# LMMER ############
#########################################

timetotoadmeta_toad_emergencemodel_median_june5_nocaptive<-lmer(MedianDays~
           Region*Competition_level +
      scaled_toadTadpoles_added_as_top.up_june7thand8th_percentagetotal+
                 X25.5.18_experiment_start_tapole_av_individual_weight+
                    Mean_averageTemp+ #fixed effects
                        (1|Population),
     data=Rawdata_02_03_20_indivrows2_median2_large,na.action=na.exclude,REML=TRUE,control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4))) 

anova(timetotoadmeta_toad_emergencemodel_median_june5_nocaptive)
AIC(timetotoadmeta_toad_emergencemodel_median_june5_nocaptive)  #441.7858


#interaction plot
timetotoadmeta_toad_emergencemodel_median_june5_nocaptive_interaction<- 
  effect('Region*Competition_level', 
              timetotoadmeta_toad_emergencemodel_median_june5_nocaptive,
                   se=TRUE)
(timetotoadmeta_toad_emergencemodel_median_june5_nocaptive_interactionplot<-
plot(timetotoadmeta_toad_emergencemodel_median_june5_nocaptive_interaction, multiline = TRUE))

#capture.output(anova(timetotoadmeta_toad_emergencemodel_median_june5_nocaptive),file="anova_timetotoadmeta_toad_emergencemodel_median_june5_nocaptive_REMLTRUE.txt")
#capture.output(timetotoadmeta_toad_emergencemodel_median_june5_nocaptive,file="timetotoadmeta_toad_emergencemodel_median_june5_nocaptive_REMLTRUE.txt")


#capture.output(anova(timetotoadmeta_toad_emergencemodel_median_june5_nocaptive),file="anova_timetotoadmeta_toad_emergencemodel_median_june5_nocaptive_REMLFALSE.txt")
#capture.output(timetotoadmeta_toad_emergencemodel_median_june5_nocaptive,file="timetotoadmeta_toad_emergencemodel_median_june5_nocaptive_REMLFALSE.txt")



############################################################################
################################## CUSTOM CONTRASTS ########################
############################################################################

#model_timetotoadmeta_toad_emergence_toads_Frank
emm3_model_timetotoadmeta_toad_emergence_mediandays=emmeans(timetotoadmeta_toad_emergencemodel_median_june5_nocaptive,specs=~ Competition_level*Region)


#capture.output(emm3_model_timetotoadmeta_toad_emergence_mediandays,file="CI_emmeans_emm3_model_timetotoadmeta_toad_emergence_mediandays_timetotoadmeta_toad_emergencemodel_median_june5_nocaptive.txt")

#make into dataframe
emm3_model_timetotoadmeta_toad_emergence_mediandaysdf<-as.data.frame(emm3_model_timetotoadmeta_toad_emergence_mediandays)

#make new column
emm3_model_timetotoadmeta_toad_emergence_mediandaysdf<-cbind(Model = "timetotoadmeta_toad_emergencemodel_median_june5_nocaptive", emm3_model_timetotoadmeta_toad_emergence_mediandaysdf)


#region HG, comp level low control - single out estimate for each level of region/competition to insect into contrast table later to calculate percent
HGLow_controlHG_model_timetotoadmeta_toad_emergencedf<-emm3_model_timetotoadmeta_toad_emergence_mediandaysdf[1,4]
HGFrogMLHG_model_timetotoadmeta_toad_emergencedf<-emm3_model_timetotoadmeta_toad_emergence_mediandaysdf[2,4]
HGFrogHGHG_model_timetotoadmeta_toad_emergencedf<-emm3_model_timetotoadmeta_toad_emergence_mediandaysdf[3,4]
HGHigh_controlHG_model_timetotoadmeta_toad_emergencedf<-emm3_model_timetotoadmeta_toad_emergence_mediandaysdf[4,4]
MLLow_controlML_model_timetotoadmeta_toad_emergencedf<-emm3_model_timetotoadmeta_toad_emergence_mediandaysdf[5,4]
MLFrogMLML_model_timetotoadmeta_toad_emergencedf<-emm3_model_timetotoadmeta_toad_emergence_mediandaysdf[6,4]
MLFrogHGML_model_timetotoadmeta_toad_emergencedf<-emm3_model_timetotoadmeta_toad_emergence_mediandaysdf[7,4]
MLHigh_controlML_model_timetotoadmeta_toad_emergencedf<-emm3_model_timetotoadmeta_toad_emergence_mediandaysdf[8,4]
HGFrog_overallHG_model_timetotoadmeta_toad_emergencedf<-(emm3_model_timetotoadmeta_toad_emergence_mediandaysdf[2,4]+emm3_model_timetotoadmeta_toad_emergence_mediandaysdf[3,4])/2
MLFrog_overallML_model_timetotoadmeta_toad_emergencedf<-(emm3_model_timetotoadmeta_toad_emergence_mediandaysdf[6,4]+emm3_model_timetotoadmeta_toad_emergence_mediandaysdf[7,4])/2


#create custom contrasts
Low_controlHG<-rep(c(1,0), times = c(1,7))
FrogMLHG<-rep(c(0,1,0), times = c(1,1,6))
FrogHGHG<-rep(c(0,1,0), times = c(2,1,5))
High_controlHG<-rep(c(0,1,0), times = c(3,1,4))
Low_controlML<-rep(c(0,1,0), times = c(4,1,3))
FrogMLML<-rep(c(0,1,0), times = c(5,1,2))
FrogHGML<-rep(c(0,1,0), times = c(6,1,1))
High_controlML<-rep(c(0,1), times = c(7,1))

#create joined frog response variable
Frog_overallHG=(FrogMLHG+FrogHGHG)/2
Frog_overallML=(FrogMLML+FrogHGML)/2


#bonferroni
(contrast_model_timetotoadmeta_toad_emergence_bonferroni2<-contrast(emm3_model_timetotoadmeta_toad_emergence_mediandays,adjust="bonferroni",method=list(
  "FrogHGHG - FrogMLHG" = FrogHGHG - FrogMLHG, #not different
  "FrogHGML - FrogMLML" = FrogHGML - FrogMLML, #not different
  "Low_controlHG - Frog_overallHG" = Low_controlHG - Frog_overallHG,
  "Low_controlML - Frog_overallML"=Low_controlML - Frog_overallML, 
  "Frog_overallHG - High_controlHG"=Frog_overallHG - High_controlHG,   
  "Frog_overallML - High_controlML"=Frog_overallML - High_controlML,
  "Low_controlHG - High_controlHG" = Low_controlHG - High_controlHG,
  "Low_controlML - High_controlML"=Low_controlML - High_controlML
  
)))

#mvt
(contrast_model_timetotoadmeta_toad_emergence_mvt<-contrast(emm3_model_timetotoadmeta_toad_emergence_mediandays,adjust="mvt",method=list(
  "FrogHGHG - FrogMLHG" = FrogHGHG - FrogMLHG, #not different
  "FrogHGML - FrogMLML" = FrogHGML - FrogMLML, #not different
  "Low_controlHG - Frog_overallHG" = Low_controlHG - Frog_overallHG,
  "Low_controlML - Frog_overallML"=Low_controlML - Frog_overallML, 
  "Frog_overallHG - High_controlHG"=Frog_overallHG - High_controlHG,   
  "Frog_overallML - High_controlML"=Frog_overallML - High_controlML,
  "Low_controlHG - High_controlHG" = Low_controlHG - High_controlHG,
  "Low_controlML - High_controlML"=Low_controlML - High_controlML
  
)))

#fdr
(contrast_model_timetotoadmeta_toad_emergence_fdr<-contrast(emm3_model_timetotoadmeta_toad_emergence_mediandays,adjust="fdr",method=list(
  "FrogHGHG - FrogMLHG" = FrogHGHG - FrogMLHG, #not different
  "FrogHGML - FrogMLML" = FrogHGML - FrogMLML, #not different
  "Low_controlHG - Frog_overallHG" = Low_controlHG - Frog_overallHG,
  "Low_controlML - Frog_overallML"=Low_controlML - Frog_overallML, 
  "Frog_overallHG - High_controlHG"=Frog_overallHG - High_controlHG,   
  "Frog_overallML - High_controlML"=Frog_overallML - High_controlML,
  "Low_controlHG - High_controlHG" = Low_controlHG - High_controlHG,
  "Low_controlML - High_controlML"=Low_controlML - High_controlML
)))

#capture
#capture.output(contrast_model_timetotoadmeta_toad_emergence_mvt,file="contrast_model_timetotoadmeta_toad_emergence_mvt_june5_nocaptive.txt")



#make into dataframe
(contrast_model_timetotoadmeta_toad_emergencedf<-as.data.frame(contrast_model_timetotoadmeta_toad_emergence_mvt))


#make new column
contrast_model_timetotoadmeta_toad_emergencedf<-cbind(Model = "timetotoadmeta_toad_emergencemodel_median_june5_nocaptive", contrast_model_timetotoadmeta_toad_emergencedf)


head(contrast_model_timetotoadmeta_toad_emergencedf)



#CIS confidence intervals

#plot emmean model estimates
(plot_eemeans_timetotoadmeta_toad_emergencemodel_median_june5_nocaptive_CI<-emmip(timetotoadmeta_toad_emergencemodel_median_june5_nocaptive, ~Competition_level| Region*Competition_level))

# CIs
(eemeans_timetotoadmeta_toad_emergencemodel_median_june5_nocaptive_CI<-emmeans(timetotoadmeta_toad_emergencemodel_median_june5_nocaptive, ~ Competition_level| Region*Competition_level))
#apture.output(eemeans_timetotoadmeta_toad_emergencemodel_median_june5_nocaptive_CI, file="eemeans_timetotoadmeta_toad_emergencemodel_median_june5_nocaptive_CI.txt")

#plot the CIs
(eemeansplot_timetotoadmeta_toad_emergencemodel_median_june5_nocaptive_CI<-plot(emmeans(timetotoadmeta_toad_emergencemodel_median_june5_nocaptive, ~ Competition_level| Region*Competition_level)))


#######################
################## PERCENTS ##################
####################


#HGLow_controlHG
(contrast_model_timetotoadmeta_toad_emergencedf$emmean_HGLow_controlHG<-rep(HGLow_controlHG_model_timetotoadmeta_toad_emergencedf,length(contrast_model_timetotoadmeta_toad_emergencedf$p.value)))
#percent
(contrast_model_timetotoadmeta_toad_emergencedf$percent_HGLow_controlHG<-(contrast_model_timetotoadmeta_toad_emergencedf$estimate/contrast_model_timetotoadmeta_toad_emergencedf$emmean_HGLow_controlHG)*100)
#value
((HGLow_controlHG_model_timetotoadmeta_toad_emergencedf_value<-contrast_model_timetotoadmeta_toad_emergencedf[c(3,7),c(1:9)]))



#MLLow_controlML
(contrast_model_timetotoadmeta_toad_emergencedf$emmean_MLLow_controlML<-rep(MLLow_controlML_model_timetotoadmeta_toad_emergencedf,length(contrast_model_timetotoadmeta_toad_emergencedf$p.value)))
#percent
(contrast_model_timetotoadmeta_toad_emergencedf$percent_MLLow_controlML<-(contrast_model_timetotoadmeta_toad_emergencedf$estimate/contrast_model_timetotoadmeta_toad_emergencedf$emmean_MLLow_controlML)*100)
#value
(MLLow_controlML_model_timetotoadmeta_toad_emergencedf_value<-contrast_model_timetotoadmeta_toad_emergencedf[c(4,8),c(1:7,10:11)])

#HGFrogoverallHG
(contrast_model_timetotoadmeta_toad_emergencedf$emmean_HGFrog_overallHG<-rep(HGFrog_overallHG_model_timetotoadmeta_toad_emergencedf,length(contrast_model_timetotoadmeta_toad_emergencedf$p.value)))
#percent
(contrast_model_timetotoadmeta_toad_emergencedf$percent_HGFrog_overallHG<-(contrast_model_timetotoadmeta_toad_emergencedf$estimate/contrast_model_timetotoadmeta_toad_emergencedf$emmean_HGFrog_overallHG)*100)
#value
((HGFrog_overallHG_model_timetotoadmeta_toad_emergencedf_value<-contrast_model_timetotoadmeta_toad_emergencedf[c(5),c(1:7,12:13)]))

#MLFrogoverallML
contrast_model_timetotoadmeta_toad_emergencedf$emmean_MLFrog_overallML<-rep(MLFrog_overallML_model_timetotoadmeta_toad_emergencedf,length(contrast_model_timetotoadmeta_toad_emergencedf$p.value))
#percent
contrast_model_timetotoadmeta_toad_emergencedf$percent_MLFrog_overallML<-(contrast_model_timetotoadmeta_toad_emergencedf$estimate/contrast_model_timetotoadmeta_toad_emergencedf$emmean_MLFrog_overallML)*100
#value
(MLFrog_overallML_model_timetotoadmeta_toad_emergencedf_value<-contrast_model_timetotoadmeta_toad_emergencedf[c(6),c(1:7,14:15)])

#write.csv(contrast_model_timetotoadmeta_toad_emergencedf,"contrast_timetotoadmeta_toad_emergencemodel_median_june5_nocaptive.csv")

#bind together
HGLow_controlHG_model_timetotoadmeta_toad_emergencedf_value
MLLow_controlML_model_timetotoadmeta_toad_emergencedf_value
HGFrog_overallHG_model_timetotoadmeta_toad_emergencedf_value
MLFrog_overallML_model_timetotoadmeta_toad_emergencedf_value

(contrast_model_timetotoadmeta_toad_emergencedf_selected<-bind_rows(
  HGLow_controlHG_model_timetotoadmeta_toad_emergencedf_value,
  MLLow_controlML_model_timetotoadmeta_toad_emergencedf_value,
  HGFrog_overallHG_model_timetotoadmeta_toad_emergencedf_value,
  MLFrog_overallML_model_timetotoadmeta_toad_emergencedf_value))



#   timetotoadmeta_toad_emergencemodel_median_june5_nocaptive
#write.csv(contrast_model_timetotoadmeta_toad_emergencedf_selected,"contrast_mvt_timetotoadmeta_toad_emergencemodel_median_june5_nocaptive_REMLTRUE.csv")




#################################
#################### PLOTS ######
################################


labels <- c(ML = "Mainland Toads", HG = "Haida Gwaii Toads")

ylim1 = boxplot.stats(Rawdata_02_03_20_indivrows2_median2_large$MedianDays)$stats[c(1, 5)]



(boxplot_timing_toad_MedianDays_03_05_20 <-Rawdata_02_03_20_indivrows2_median2_large %>%
    mutate(Competition_level = fct_relevel(Treatment_three, 
                                           "Low_control", "Frog","High_control")) %>%
    ggplot(aes(x=Competition_level, y=MedianDays)) +
    geom_boxplot(lwd=0.6)+
    scale_fill_manual("grey")+
    facet_wrap(~Region,labeller = labeller(Region = labels))+
    labs(x="Competition Treatment")+
    labs(y="Median time to toad metamorphosis (days)")+
    scale_x_discrete(breaks=c("Low_control",  "Frog","High_control"),
                     labels=c("Low intraspecific\n80 toads","Interspecific\n16 frogs\n80 toads", "High intraspecific\n160 toads" ))+
    scale_colour_discrete(name="Toad Region",breaks=c("HG", "ML"),
                          labels=c("Haida Gwaii", "Mainland"))+
    theme_bw()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()))


#####################################################################
##############################################################################
########## 4. Toad Mortality #################################################
##############################################################################
##############################################################################

rm(list=ls())

#same data as Rawdata_02_03_20 but only have one mortality point per tank, so it's a short dataframe (don't have repeated measures of mortality per tank)

subset_dataRGRdd_renamedsomecolumns_2_10_20<-read.csv(file="subset_dataRGRdd_renamedsomecolumns_2_10_20.csv",header=T,row.names=NULL,sep=",")

#order levels of factor
subset_dataRGRdd_renamedsomecolumns_2_10_20$Competition_level<- factor(subset_dataRGRdd_renamedsomecolumns_2_10_20$Competition_level, levels = c("Low_control","FrogML","FrogHG","High_control"))
subset_dataRGRdd_renamedsomecolumns_2_10_20<- droplevels(subset_dataRGRdd_renamedsomecolumns_2_10_20)


subset_dataRGRdd_renamedsomecolumns_2_10_20$toadTadpoles_added_as_top.up_june7thand8th_percentagetotal<-(subset_dataRGRdd_renamedsomecolumns_2_10_20$Tadpoles_added_as_top.up_june7thand8th/subset_dataRGRdd_renamedsomecolumns_2_10_20$Total_number_individuals_start)*100

subset_dataRGRdd_renamedsomecolumns_2_10_20$Mortality_05_05_20_tadsadded<-subset_dataRGRdd_renamedsomecolumns_2_10_20$Total_number_individuals_start+subset_dataRGRdd_renamedsomecolumns_2_10_20$Tadpoles_added_as_top.up_june7thand8th-(subset_dataRGRdd_renamedsomecolumns_2_10_20$Running_total_of_adults_emerged+subset_dataRGRdd_renamedsomecolumns_2_10_20$X23.8.18_endexperiment_number_tadpoles+subset_dataRGRdd_renamedsomecolumns_2_10_20$X23.8.18_endexperiment_number_metamorphs+subset_dataRGRdd_renamedsomecolumns_2_10_20$X23.8.18_number_of_metamorphs)

subset_dataRGRdd_renamedsomecolumns_2_10_20$Mortality_05_05_20_tadsadded_percapita<-subset_dataRGRdd_renamedsomecolumns_2_10_20$Mortality_05_05_20_tadsadded/(subset_dataRGRdd_renamedsomecolumns_2_10_20$Total_number_individuals_start +subset_dataRGRdd_renamedsomecolumns_2_10_20$Tadpoles_added_as_top.up_june7thand8th)


#######################################
############# Linear Mixed Effects Model - toad mortality ############
#########################################


Model_mortality_tadsadded_percapita_june9th<-lmer(Mortality_05_05_20_tadsadded_percapita~
                     Region*Competition_level + 
                X25.5.18_experiment_start_tapole_av_individual_weight+
                      Mean_averageTemp+ #fixed effects
                          (1|Population),
   data=subset_dataRGRdd_renamedsomecolumns_2_10_20,na.action=na.exclude,REML=TRUE,control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4))) 

anova(Model_mortality_tadsadded_percapita_june9th)
AIC(Model_mortality_tadsadded_percapita_june9th) #3.34038
visreg(Model_mortality_tadsadded_percapita_june9th)

#capture.output(Model_mortality_tadsadded_percapita_june9th,file="Model_mortality_tadsadded_percapita_june9th_REMLTRUE.txt")
#capture.output(anova(Model_mortality_tadsadded_percapita_june9th),file="anova_Model_mortality_tadsadded_percapita_june9th_REMLTRUE.txt")

#capture.output(Model_mortality_tadsadded_percapita_june9th,file="Model_mortality_tadsadded_percapita_june9th_REMLFALSE.txt")
#capture.output(anova(Model_mortality_tadsadded_percapita_june9th),file="anova_Model_mortality_tadsadded_percapita_june9th_REMLFALSE.txt")



####################################################################
################### CUSTOM CONTRASTS ################
####################################################################

emm3_Model_mortality_tadsadded_percapita=emmeans(Model_mortality_tadsadded_percapita_june9th,specs=~ Competition_level*Region)

#capture.output(emm3_Model_mortality_tadsadded_percapita,file="emm3_Model_mortality_tadsadded_percapita_june9th.txt")



#make into dataframe
emm3_Model_mortality_tadsadded_percapitadf<-as.data.frame(emm3_Model_mortality_tadsadded_percapita)

#make new column
emm3_Model_mortality_tadsadded_percapitadf<-cbind(Model = "Model_mortality_tadsadded_percapita_june9th", emm3_Model_mortality_tadsadded_percapitadf)


#region HG, comp level low control - single out estimate for each level of region/competition to insect into contrast table later to calculate percent
HGLow_controlHG_Model_mortality_tadsadded_percapitadf<-emm3_Model_mortality_tadsadded_percapitadf[1,4]
HGFrogMLHG_Model_mortality_tadsadded_percapitadf<-emm3_Model_mortality_tadsadded_percapitadf[2,4]
HGFrogHGHG_Model_mortality_tadsadded_percapitadf<-emm3_Model_mortality_tadsadded_percapitadf[3,4]
HGHigh_controlHG_Model_mortality_tadsadded_percapitadf<-emm3_Model_mortality_tadsadded_percapitadf[4,4]
MLLow_controlML_Model_mortality_tadsadded_percapitadf<-emm3_Model_mortality_tadsadded_percapitadf[5,4] 
MLFrogMLML_Model_mortality_tadsadded_percapitadf<-emm3_Model_mortality_tadsadded_percapitadf[6,4]
MLFrogHGML_Model_mortality_tadsadded_percapitadf<-emm3_Model_mortality_tadsadded_percapitadf[7,4]
MLHigh_controlML_Model_mortality_tadsadded_percapitadf<-emm3_Model_mortality_tadsadded_percapitadf[8,4]
HGFrog_overallHG_Model_mortality_tadsadded_percapitadf<-(emm3_Model_mortality_tadsadded_percapitadf[2,4]+emm3_Model_mortality_tadsadded_percapitadf[3,4])/2
MLFrog_overallML_Model_mortality_tadsadded_percapitadf<-(emm3_Model_mortality_tadsadded_percapitadf[6,4]+emm3_Model_mortality_tadsadded_percapitadf[7,4])/2


#create custom contrasts
Low_controlHG<-rep(c(1,0), times = c(1,7))
FrogMLHG<-rep(c(0,1,0), times = c(1,1,6))
FrogHGHG<-rep(c(0,1,0), times = c(2,1,5))
High_controlHG<-rep(c(0,1,0), times = c(3,1,4))
Low_controlML<-rep(c(0,1,0), times = c(4,1,3))
FrogMLML<-rep(c(0,1,0), times = c(5,1,2))
FrogHGML<-rep(c(0,1,0), times = c(6,1,1))
High_controlML<-rep(c(0,1), times = c(7,1))

#create joined frog response variable
Frog_overallHG=(FrogMLHG+FrogHGHG)/2
Frog_overallML=(FrogMLML+FrogHGML)/2



# adjust = "mvt" adds multiple tests adjustment, could also use "bonferroni" but its too conservative https://cran.r-project.org/web/packages/emmeans/vignettes/confidence-intervals.html

(contrast_Model_mortality_tadsadded_percapita<-contrast(emm3_Model_mortality_tadsadded_percapita,adjust="mvt",method=list(
  "FrogHGHG - FrogMLHG" = FrogHGHG - FrogMLHG, #not different
  "FrogHGML - FrogMLML" = FrogHGML - FrogMLML, #not different
  "Low_controlHG - Frog_overallHG" = Low_controlHG - Frog_overallHG,
  "Low_controlML - Frog_overallML"=Low_controlML - Frog_overallML, 
  "Frog_overallHG - High_controlHG"=Frog_overallHG - High_controlHG,   
  "Frog_overallML - High_controlML"=Frog_overallML - High_controlML,
  "(Low_controlHG - Frog_overallHG)-(Low_controlML - Frog_overallML)" = 
    ((Low_controlHG - Frog_overallHG)-(Low_controlML - Frog_overallML)),
  "(Frog_overallHG - High_controlHG)-(Frog_overallML - High_controlML)" = 
    ((Frog_overallHG - High_controlHG)-(Frog_overallML - High_controlML))
))) #nothing significant



(contrast_Model_mortality_tadsadded_percapita_bonferroni2<-contrast(emm3_Model_mortality_tadsadded_percapita,adjust="bonferroni",method=list(
  "FrogHGHG - FrogMLHG" = FrogHGHG - FrogMLHG, #not different
  "FrogHGML - FrogMLML" = FrogHGML - FrogMLML, #not different
  "Low_controlHG - Frog_overallHG" = Low_controlHG - Frog_overallHG,
  "Low_controlML - Frog_overallML"=Low_controlML - Frog_overallML, 
  "Frog_overallHG - High_controlHG"=Frog_overallHG - High_controlHG,   
  "Frog_overallML - High_controlML"=Frog_overallML - High_controlML,
  "Low_controlHG - High_controlHG" = Low_controlHG - High_controlHG,
  "Low_controlML - High_controlML"=Low_controlML - High_controlML
)))

#mvt
(contrast_Model_mortality_tadsadded_percapita_mvt<-contrast(emm3_Model_mortality_tadsadded_percapita,adjust="mvt",method=list(
  "FrogHGHG - FrogMLHG" = FrogHGHG - FrogMLHG, #not different
  "FrogHGML - FrogMLML" = FrogHGML - FrogMLML, #not different
  "Low_controlHG - Frog_overallHG" = Low_controlHG - Frog_overallHG,
  "Low_controlML - Frog_overallML"=Low_controlML - Frog_overallML, 
  "Frog_overallHG - High_controlHG"=Frog_overallHG - High_controlHG,   
  "Frog_overallML - High_controlML"=Frog_overallML - High_controlML,
  "Low_controlHG - High_controlHG" = Low_controlHG - High_controlHG,
  "Low_controlML - High_controlML"=Low_controlML - High_controlML
)))


#make into dataframe
(contrast_Model_mortality_tadsadded_percapitadf<-as.data.frame(contrast_Model_mortality_tadsadded_percapita_mvt))

#capture.output(contrast_Model_mortality_tadsadded_percapita_mvt,file="contrast_Model_mortality_tadsadded_percapita_mvt.txt")


#make new column
contrast_Model_mortality_tadsadded_percapitadf<-cbind(Model = "Model_mortality_tadsadded_percapita_june9th", contrast_Model_mortality_tadsadded_percapitadf)


#######################
################## PERCENTS #############
####################



#HGLow_controlHG
contrast_Model_mortality_tadsadded_percapitadf$emmean_HGLow_controlHG<-rep(HGLow_controlHG_Model_mortality_tadsadded_percapitadf,length(contrast_Model_mortality_tadsadded_percapitadf$p.value))
#percent
contrast_Model_mortality_tadsadded_percapitadf$percent_HGLow_controlHG<-(contrast_Model_mortality_tadsadded_percapitadf$estimate/contrast_Model_mortality_tadsadded_percapitadf$emmean_HGLow_controlHG)*100
#value
(HGLow_controlHG_Model_mortality_tadsadded_percapitadf_value<-contrast_Model_mortality_tadsadded_percapitadf[c(3,7),c(1:9)])



#MLLow_controlML
contrast_Model_mortality_tadsadded_percapitadf$emmean_MLLow_controlML<-rep(MLLow_controlML_Model_mortality_tadsadded_percapitadf,length(contrast_Model_mortality_tadsadded_percapitadf$p.value))
#percent
contrast_Model_mortality_tadsadded_percapitadf$percent_MLLow_controlML<-(contrast_Model_mortality_tadsadded_percapitadf$estimate/contrast_Model_mortality_tadsadded_percapitadf$emmean_MLLow_controlML)*100
#value
(MLLow_controlML_Model_mortality_tadsadded_percapitadf_value<-contrast_Model_mortality_tadsadded_percapitadf[c(4,8),c(1:7,10:11)])

#HGFrogoverallHG
contrast_Model_mortality_tadsadded_percapitadf$emmean_HGFrog_overallHG<-rep(HGFrog_overallHG_Model_mortality_tadsadded_percapitadf,length(contrast_Model_mortality_tadsadded_percapitadf$p.value))
#percent
contrast_Model_mortality_tadsadded_percapitadf$percent_HGFrog_overallHG<-(contrast_Model_mortality_tadsadded_percapitadf$estimate/contrast_Model_mortality_tadsadded_percapitadf$emmean_HGFrog_overallHG)*100
#value
(HGFrog_overallHG_Model_mortality_tadsadded_percapitadf_value<-contrast_Model_mortality_tadsadded_percapitadf[c(5),c(1:7,12:13)])

#MLFrogoverallML
contrast_Model_mortality_tadsadded_percapitadf$emmean_MLFrog_overallML<-rep(MLFrog_overallML_Model_mortality_tadsadded_percapitadf,length(contrast_Model_mortality_tadsadded_percapitadf$p.value))
#percent
contrast_Model_mortality_tadsadded_percapitadf$percent_MLFrog_overallML<-(contrast_Model_mortality_tadsadded_percapitadf$estimate/contrast_Model_mortality_tadsadded_percapitadf$emmean_MLFrog_overallML)*100
#value
(MLFrog_overallML_Model_mortality_tadsadded_percapitadf_value<-contrast_Model_mortality_tadsadded_percapitadf[c(6),c(1:7,14:15)])

#write.csv(contrast_Model_mortality_tadsadded_percapitadf,"contrast_Model_mortality_tadsadded_percapitadf.csv")

#bind together
HGLow_controlHG_Model_mortality_tadsadded_percapitadf_value
MLLow_controlML_Model_mortality_tadsadded_percapitadf_value
HGFrog_overallHG_Model_mortality_tadsadded_percapitadf_value
MLFrog_overallML_Model_mortality_tadsadded_percapitadf_value

(contrast_Model_mortality_tadsadded_percapitadf_selected<-bind_rows(
  HGLow_controlHG_Model_mortality_tadsadded_percapitadf_value,
  MLLow_controlML_Model_mortality_tadsadded_percapitadf_value,
  HGFrog_overallHG_Model_mortality_tadsadded_percapitadf_value,
  MLFrog_overallML_Model_mortality_tadsadded_percapitadf_value))

#Model_mortality_tadsadded_percapita_june9th
#write.csv(contrast_Model_mortality_tadsadded_percapitadf_selected,"contrast_Model_mortality_tadsadded_percapitadf_selected_mvt_REMLTRUE.csv")

###############################################
#################### PLOT ###################
##############################################
labels <- c(ML = "Mainland Toads", HG = "Haida Gwaii Toads")


ylim1 = boxplot.stats(subset_dataRGRdd_renamedsomecolumns_2_10_20$Mortality_05_05_20_tadsadded_percapita)$stats[c(1, 5)]

(boxplot_Mortality_05_05_20_tadsadded_percapita_06_05_20 <-subset_dataRGRdd_renamedsomecolumns_2_10_20 %>%
    mutate(Competition_level = fct_relevel(Treatment_three, 
                                           "Low_control", "Frog","High_control")) %>%
    ggplot(aes(x=Competition_level, y=Mortality_05_05_20_tadsadded_percapita)) +
    geom_boxplot(lwd=0.6) +
    scale_fill_manual("grey")+
    coord_cartesian(ylim = ylim1*1.05)+
    facet_wrap(~Region,labeller = labeller(Region = labels))+
    labs(x="Competition Treatment")+
    labs(y="Toad mortality per capita")+
    scale_x_discrete(breaks=c("Low_control",  "Frog","High_control"),
                     labels=c("Low intraspecific\n80 toads","Interspecific\n16 frogs\n80 toads", "High intraspecific\n160 toads" ))+
    scale_colour_discrete(name="Toad Region",breaks=c("HG", "ML"),
                          labels=c("Haida Gwaii", "Mainland"))+
    theme_bw()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()))

