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

#  ---------------------------  TESTING FOCAL DATA  ---------------------------

# Focal Data used is called "Focal Sampling Data - combined BORIS data.xlsx" found in BOX
data <- read.csv("finalfocal.csv")
#View(data)
data <- data %>%
  group_by(Trial, Species, Treatment, Period) %>%
  mutate(total_time_in_sight = sum(Duration[Behavior %in% c("move", "forage", "look", "other")])) %>% 
  ungroup()

# This subsets data for each behavior.
forage_data <- data %>% filter(Behavior == "forage")
look_data <- data %>% filter(Behavior == "look")
move_data <- data %>% filter(Behavior == "move")


Sforage_model <- glmmTMB(Duration ~ Treatment * Period + Species + (1|Trial),
                         family = ziGamma(link = "log"),
                         ziformula = ~ 1,
                         data = forage_data,
                         offset = log(total_time_in_sight+1e-6))

Slook_model <- glmmTMB(Duration ~ Treatment * Period + Species + (1|Trial),
                       family = ziGamma(link = "log"),
                       ziformula = ~ 1,
                       data = look_data,
                       offset = log(total_time_in_sight+1e-6))

Smove_model <- glmmTMB(Duration ~ Treatment * Period + Species + (1|Trial),
                       family = ziGamma(link = "log"),
                       ziformula = ~ 1,
                       data = move_data,
                       offset = log(total_time_in_sight+1e-6))
summary(Sforage_model)
summary(Slook_model) 
summary(Smove_model)  

AIC(Sforage_model, Slook_model, Smove_model)


# ---- observer added to species model as RANDOM effect ----
OSforage_model <- glmmTMB(Duration ~ Treatment * Period + Species + (1|Trial) + (1|Observer),
                          family = ziGamma(link = "log"),
                          ziformula = ~ 1,
                          data = forage_data,
                          offset = log(total_time_in_sight+1e-6))

OSlook_model <- glmmTMB(Duration ~ Treatment * Period + Species + (1|Trial) + (1|Observer),
                        family = ziGamma(link = "log"),
                        ziformula = ~ 1,
                        data = look_data,
                        offset = log(total_time_in_sight+1e-6))

OSmove_model <- glmmTMB(Duration ~ Treatment * Period + Species + (1|Trial) + (1|Observer),
                        family = ziGamma(link = "log"),
                        ziformula = ~ 1,
                        data = move_data,
                        offset = log(total_time_in_sight+1e-6))
summary(OSforage_model) #AIC is higher 1224>1222
summary(OSlook_model) #AIC is higher 1014>1013
summary(OSmove_model) #AIC is higher 1330>1328
AIC(OSforage_model, Sforage_model, OSlook_model, Slook_model, OSmove_model, Smove_model)

# no sig p-value
anova(Sforage_model,OSforage_model)
anova(Slook_model,OSlook_model)
anova(Smove_model,OSmove_model)

# ---- observer added to species model as FIXED effect ----
OFSforage_model <- glmmTMB(Duration ~ Treatment * Period + Species + Observer + (1|Trial),
                          family = ziGamma(link = "log"),
                          ziformula = ~ 1,
                          data = forage_data,
                          offset = log(total_time_in_sight+1e-6))

OFSlook_model <- glmmTMB(Duration ~ Treatment * Period + Species + Observer + (1|Trial),
                        family = ziGamma(link = "log"),
                        ziformula = ~ 1,
                        data = look_data,
                        offset = log(total_time_in_sight+1e-6))

OFSmove_model <- glmmTMB(Duration ~ Treatment * Period + Species + Observer + (1|Trial),
                        family = ziGamma(link = "log"),
                        ziformula = ~ 1,
                        data = move_data,
                        offset = log(total_time_in_sight+1e-6))
summary(OFSforage_model) #AIC is higher 1223>1222
summary(OFSlook_model) #AIC is higher 1011>1013
summary(OFSmove_model) #AIC is higher 1332>1328
AIC(OFSforage_model, Sforage_model, OFSlook_model, Slook_model, OFSmove_model, Smove_model)

look_data %>% count(Observer)

anova(Sforage_model,OFSforage_model)
anova(Slook_model,OFSlook_model) # P value 0.03642 *
anova(Smove_model,OFSmove_model)


#  --------------------------- TESTING SCAN DATA  ---------------------------
gdata <- read.csv("finalscan.csv")
gdata <- gdata %>% select(-Ethogram, -X.1)

nb_forage <- glmmTMB(Forage ~ Species+Period * Treatment + (1|Trial), 
                     family = nbinom2, 
                     data = gdata,
                     offset = log(Total+1))
nb_look <- glmmTMB(Look ~ Species+Period * Treatment + (1|Trial), 
                   family = nbinom2, 
                   data = gdata,
                   offset = log(Total+1))

nb_move <- glmmTMB(Move ~ Species+Period * Treatment + (1|Trial), 
                   family = nbinom2, 
                   data = gdata,
                   offset = log(Total+1))
# looking at TreatmentExperimental
summary(nb_forage) # 1.18e-07 *
summary(nb_move) # 0.0353*
summary(nb_look) # 0.0125*



# ---- observer added to species model as RANDOM effect ----
nb_or_forage <- glmmTMB(Forage ~ Species+Period * Treatment + (1|Trial) + (1|Observer), 
                     family = nbinom2, 
                     data = gdata,
                     offset = log(Total+1))
nb_or_look <- glmmTMB(Look ~ Species+Period * Treatment + (1|Trial) + (1|Observer), 
                   family = nbinom2, 
                   data = gdata,
                   offset = log(Total+1))

nb_or_move <- glmmTMB(Move ~ Species+Period * Treatment + (1|Trial) + (1|Observer), 
                   family = nbinom2, 
                   data = gdata,
                   offset = log(Total+1))

AIC(nb_or_forage, nb_forage, nb_or_look, nb_look, nb_or_move, nb_move)
# Forage with observe as random effect reduced AIC by a decent amount, 621 > 614 (observer added)
# Look and Move both better AIC using original model w/o observer effect

anova(nb_forage, nb_or_forage, test="LRT")
anova(nb_look, nb_or_look, test="LRT")
anova(nb_move, nb_or_move, test="LRT")

# ---- observer added to species model as FIXED effect ----
nb_of_forage <- glmmTMB(Forage ~ Species+Period * Treatment + Observer + (1|Trial), 
                     family = nbinom2, 
                     data = gdata,
                     offset = log(Total+1))
nb_of_look <- glmmTMB(Look ~ Species+Period * Treatment + Observer + (1|Trial), 
                   family = nbinom2, 
                   data = gdata,
                   offset = log(Total+1))

nb_of_move <- glmmTMB(Move ~ Species+Period * Treatment + Observer + (1|Trial), 
                   family = nbinom2, 
                   data = gdata,
                   offset = log(Total+1))
anova(nb_forage, nb_of_forage, test="LRT")
anova(nb_look, nb_of_look, test="LRT")
anova(nb_move, nb_of_move, test="LRT")

AIC(nb_of_forage, nb_forage, nb_of_look, nb_look, nb_of_move, nb_move)
# Forage with observe as random effect reduced AIC by a decent amount, 621 > 606.67 (observer added)

summary(nb_of_forage)
# ObserverTavish                    0.68847    0.14247   4.832 1.35e-06 ***
summary(nb_of_look)
# nothin
summary(nb_of_move)
# ObserverTavish                   -0.160491   0.080211  -2.001   0.0454 * 