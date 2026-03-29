getwd()
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

#  ---------------------------  FINAL FOCAL DATA  ---------------------------
# Using final focal data, I'm creating a new version of the dataset with 
# only experimental playback period with the trials
habituation_data <- read.csv("finalfocal.csv")
#view(habituation_data)
# Removing out of sight and other


# Finding total time in sight since it varies, used for offset() later.
habituation_data <- habituation_data %>%
  group_by(Trial, Species, Treatment, Period) %>%
  mutate(total_time_in_sight = 
           sum(Duration[Behavior %in% c("move", "forage", "look", "other")])) %>% 
  ungroup()

# TESTING FOR HOW BEHAVIORS CHANGE
# This subsets data for each behavior.
forage_hab_data <-habituation_data %>% filter(Behavior == "forage",
                                              Treatment == "experimental", Period == "playback")
look_hab_data <- habituation_data %>% filter(Behavior == "look",
                                             Treatment == "experimental", Period == "playback")
move_hab_data <- habituation_data %>% filter(Behavior == "move",
                                             Treatment == "experimental", Period == "playback")
#view(forage_hab_data)

forage_hab_test <- glmmTMB(Duration ~ Trial  + (1|Species),
                           family = ziGamma(link = "log"),
                           ziformula = ~ 1,
                           data = forage_hab_data,
                           offset = log(total_time_in_sight+1e-6))

look_hab_test <- glmmTMB(Duration ~ Trial + (1|Species),
                           family = ziGamma(link = "log"),
                           ziformula = ~ 1,
                           data = look_hab_data,
                           offset = log(total_time_in_sight+1e-6))


move_hab_test <- glmmTMB(Duration ~ Trial + (1|Species),
                           family = ziGamma(link = "log"),
                           ziformula = ~ 1,
                           data = move_hab_data,
                           offset = log(total_time_in_sight+1e-6))
summary(forage_hab_test)
summary(look_hab_test)
summary(move_hab_test)


# NOW TESTING IF THEY GOT MORE OR LESS VISIBLE AS TRIALS CONTINUED
visibility_data <- habituation_data %>% 
  filter (Treatment == "experimental", Period == "playback") %>%
  group_by(Trial, Species) %>%
  summarise(total_time_in_sight = first(total_time_in_sight),.groups="drop")
#?summarise  

visibility_test <- glmmTMB(total_time_in_sight ~ Trial + (1|Species),
                           family = ziGamma(link = "log"),
                           ziformula = ~ 1,
                           data = visibility_data) # no offset, unsure(?)
summary(visibility_test)


# ---------------------------  SCAN FOCAL DATA  ---------------------------

habituation_scandata <- read.csv("finalscan.csv")

habituation_scandata <- habituation_scandata %>% select(-Ethogram, -X.1) 
# earlier version of the data had two random columns with ethogram and a blank -X.1 labeled 


# filter for only playback and experimental again
scandata_playback <- habituation_scandata %>% filter(Period == "Playback",
                                                     Treatment == "Experimental")
View(scandata_playback)

visibility_test_scan <- glmmTMB (Total ~ Trial + (1|Species),
                                family = nbinom2,
                                data = scandata_playback)
summary(visibility_test_scan)

forage_test_scan <- glmmTMB (Forage ~ Trial + (1|Species),
                                 family = nbinom2,
                                 data = scandata_playback)
look_test_scan <- glmmTMB (Look ~ Trial + (1|Species),
                             family = nbinom2,
                             data = scandata_playback)
move_test_scan <- glmmTMB (Move ~ Trial + (1|Species),
                             family = nbinom2,
                             data = scandata_playback)


summary(forage_test_scan)
summary(look_test_scan)
summary(move_test_scan)

