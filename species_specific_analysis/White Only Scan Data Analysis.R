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
library("tidyverse")
library("emmeans")
library("partR2")
library("ggplot2")
library("ggpubr")
library("here")
library("moments")
library("dplyr")
rm(list = ls())
# FOR WHITE
    # foraging is still lower during playback when compared to pre and post in experimental
    # looking insignificant
    # movement still higher during playback when compared to pre in experimental

# Scan Data 
gdata <- read.csv("data/finalscan.csv")
gdata <- gdata %>% select(-Ethogram, -X.1) # earlier version of the data had two random columns with ethogram and a blank -X.1 labeled 
gdata <- gdata %>% filter(Trial !=2) # remove mixed
# view(gdata)
# n_distinct(gdata$Trial)

w_s_data <- gdata %>% filter(Species == "w")
#view(w_s_data)

# Since we have count data, run Poisson first.

# ---------- Poisson Test ----------
p_forage <- glmmTMB(Forage ~ Period * Treatment + (1|Trial), 
                    family = poisson, 
                    data = w_s_data,
                    offset = log(Total+1))
res_p_forage <- simulateResiduals(p_forage)
testDispersion(res_p_forage) # no dispersion problem

p_look <- glmmTMB(Look ~ Period * Treatment + (1|Trial), 
                  family = poisson, 
                  data = w_s_data,
                  offset = log(Total+1))
res_p_look <- simulateResiduals(p_look)
testDispersion(res_p_look) # no dispersion problem

p_move <- glmmTMB(Move ~ Period * Treatment + (1|Trial), 
                  family = poisson, 
                  data = w_s_data,
                  offset = log(Total+1))
res_p_move <- simulateResiduals(p_move)
testDispersion(res_p_move) # no dispersion problem

# just incase, I'll still test negative binomial

# ---------- Negative Binomial Test ----------
nb_forage <- glmmTMB(Forage ~ Period * Treatment + (1|Trial), 
                     family = nbinom2, 
                     data = w_s_data,
                     offset = log(Total+1)) # false convergence
nb_look <- glmmTMB(Look ~ Period * Treatment + (1|Trial), 
                   family = nbinom2, 
                   data = w_s_data,
                   offset = log(Total+1))

nb_move <- glmmTMB(Move ~ Period * Treatment + (1|Trial), 
                   family = nbinom2, 
                   data = w_s_data,
                   offset = log(Total+1)) 
AIC (p_forage, nb_forage, p_look, nb_look, p_move, nb_move)
# forage is same, within 1 unit difference, we will use poisson since nb gave convergence error
# look is same, withint 1 AIC unit so I will stick to poisson
# move is 2 AIC units lower for poisson, I will use poisson


# ---------- FINAL SCAN MODELS ----------
summary(p_forage) # 0.00286 ** 
summary(p_look) # 0.929    
summary(p_move) # 0.00817 ** 
# TreatmentExperimental listed

# ---------- remove interaction and ANOVA/emmeans ----------
p_sansint_forage <- glmmTMB(Forage ~  Period + Treatment + (1|Trial), 
                            family = poisson, 
                            data = w_s_data,
                            offset = log(Total+1))
p_sansint_look <- glmmTMB(Look ~  Period + Treatment + (1|Trial), 
                           family = poisson, 
                           data = w_s_data,
                           offset = log(Total+1))
p_sansint_move <- glmmTMB(Move ~  Period + Treatment + (1|Trial), 
                          family = poisson, 
                          data = w_s_data,
                          offset = log(Total+1))
summary(p_sansint_forage) # 0.043397 *
summary(p_sansint_look) # 0.724     
summary(p_sansint_move) # 0.04721 *  

# ANOVA
anova(p_forage, p_sansint_forage, test="LRT") # 0.04136 *
anova(p_look, p_sansint_look, test="LRT") # 0.8797
anova(p_move, p_sansint_move, test="LRT") # 0.08868 . change from full scan data, this isn't significant anymore

forage_emm <- emmeans(p_forage, pairwise ~ Period | Treatment)
summary(forage_emm) # 0.0015 Playback - Post  (z = -3.477), 0.0008 playback-pre (z = -3.639)
# they foraged less during the playback than both post and pre!!

look_emm <- emmeans(p_look, pairwise ~ Period | Treatment)
summary(look_emm) # nothing

move_emm <- emmeans(p_move, pairwise ~ Period | Treatment)
summary(move_emm) # 0.0029 Playback - Pre (3.286), they moved more during playback than pre-playback

# FINAL CONCLUSION
  # For white-fronted capuchins, foraging rates decreased during playback 
  # relative to both pre- and post-playback periods in the experimental treatment,
  # with no comparable effects seen in the control treatment.
  # This follows the original scan data conclusion.

