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
library("readxl")
library("psych")

# FOCAL MODEL BY BROWN ----
data <- read.csv("finalfocal.csv")

# Trial 2 was mixed species group so removing it
data_wo_mixed <- data %>% filter(Trial !=2)
  #view(data_wo_mixed)

# Checking
any(data_wo_mixed$Trial ==2) #checking to see if any of trial 2 left
n_distinct(data_wo_mixed$Trial) #21 trials, which is accurate 

# Brown
b_f_data <- data_wo_mixed%>% filter(Species == "b")
  #view(b_f_data)
n_distinct(b_f_data$Trial) # 12

data <- b_f_data

# ANALYSIS ----
# Finding total time in sight since it varies, used for offset() later.
data <- data %>%
  group_by(Trial, Treatment, Period) %>%
  mutate(total_time_in_sight = sum(Duration[Behavior %in% c("move", "forage", "look", "other")])) %>% 
  ungroup()

# This subsets data for each behavior.
forage_data <- data %>% filter(Behavior == "forage")
look_data <- data %>% filter(Behavior == "look")
move_data <- data %>% filter(Behavior == "move")

# FINAL MODEL ----
Sforage_model <- glmmTMB(Duration ~ Treatment * Period + (1|Trial),
                         family = ziGamma(link = "log"),
                         ziformula = ~ 1,
                         data = forage_data,
                         offset = log(total_time_in_sight+1e-6))

Slook_model <- glmmTMB(Duration ~ Treatment * Period + (1|Trial),
                       family = ziGamma(link = "log"),
                       ziformula = ~ 1,
                       data = look_data,
                       offset = log(total_time_in_sight+1e-6))

Smove_model <- glmmTMB(Duration ~ Treatment * Period + (1|Trial),
                       family = ziGamma(link = "log"),
                       ziformula = ~ 1,
                       data = move_data,
                       offset = log(total_time_in_sight+1e-6))
summary(Sforage_model)
summary(Slook_model) 
summary(Smove_model) 

Sforage_sansint <- glmmTMB(Duration ~ Treatment + Period + (1|Trial),
                           family = ziGamma(link = "log"),
                           ziformula = ~ 1,
                           data = forage_data,
                           offset = log(total_time_in_sight+1e-6))
Smove_sansint <- glmmTMB(Duration ~ Treatment + Period+ (1|Trial),
                         family = ziGamma(link = "log"),
                         ziformula = ~ 1,
                         data = move_data,
                         offset = log(total_time_in_sight+1e-6))
Slook_sansint <- glmmTMB(Duration ~ Treatment + Period + (1|Trial),
                         family = ziGamma(link = "log"),
                         ziformula = ~ 1,
                         data = look_data,
                         offset = log(total_time_in_sight+1e-6))
summary(Sforage_sansint) # p=0.000885 ***
summary(Slook_sansint) # p=0.028 *  
summary(Smove_sansint) # p=0.0369 *  

anova(Sforage_model, Sforage_sansint, test="LRT") # p=0.006092 ** YAY
anova(Slook_model, Slook_sansint, test="LRT") # p=0.7925 eh
anova(Smove_model, Smove_sansint, test="LRT") # p=0.01858 * YAY

# pairwise, only FORAGE was SIGNIFICANT   0.0062 ----
Sforage_emm <- emmeans(Sforage_model, pairwise ~ Period | Treatment)
Sforage_emm 
  #plot(Sforage_emm)
# playback - pre experimental contrast p = 0.0062 !!, none of the control period contrasts are significant as expected
# this should mean CAPUCHINS FORAGE SIG. LESS DURING EXPERIMENTAL PLAYBACK
Slook_emm <- emmeans(Slook_model, pairwise ~ Period | Treatment)
Slook_emm
  #plot(Slook_emm)
Smove_emm <- emmeans(Smove_model, pairwise ~ Period | Treatment)
Smove_emm
  #plot(Smove_emm)




