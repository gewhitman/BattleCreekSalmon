#Age validation with FSA

# Prepare the environment -------------------------------------------------
library(FSA)
library(here)
library(tidyverse)
library(ggthemes)
library(cowplot)
library(openxlsx) # read in excel files .xlsx


#Clear the workspace
rm(list = ls())

#Read in data
bc_ages_readin <-read.csv ("data/battle_data.csv")
bc_age_comparison_readin <- read.xlsx(here('data','battlecreek_ages_tsc_20241113.xlsx'), sheet="ages_compare")

#Prepare data
agedata <- bc_ages_readin%>%
  left_join(bc_age_comparison_readin, by="lab_id")%>%
  distinct(lab_id,age, lt_age,ke_age)

## Age Bias plot Reader1 vs Reader2
agedata_1v2 <- agedata %>%dplyr::select (lab_id,lt_age, ke_age)%>%
na.omit(agedata_1v2) 

ab.tA2 <- ageBias(lt_age~ke_age, data=agedata_1v2,
                  ref.lab=" ke_age",nref.lab="lt_age")
plot(ab.tA2)

## Age precision
ap.A <- agePrecision(~lt_age+ke_age,data=agedata_1v2)
summary(ap.A,what="difference")
summary(ap.A,what="precision")

