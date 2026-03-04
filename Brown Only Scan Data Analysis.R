setwd("/Users/anyangyu/Desktop/CAP Code and Results")

library("glmmTMB")
library("lme4")
library("lmerTest")
library("performance")
library("car")
library("MuMIn")
library("DHARMa")
library("see")
library("patchwork")
library("ggpattern")
library("egg")
library("car")
library("tidyverse")
library("emmeans")
library("partR2")
library("ggplot2")
library("ggpubr")
library("egg")
library("here")
library("moments")
library("dplyr")

# Scan Data 
gdata <- read.csv("finalscan.csv")
gdata <- gdata %>% select(-Ethogram, -X.1) # earlier version of the data had two random columns with ethogram and a blank -X.1 labeled 
gdata <- gdata %>% filter(Trial !=2) # remove mixed
  #view(gdata)
  #n_distinct(gdata$Trial)

b_s_data <- gdata %>% filter(Species == "b")

