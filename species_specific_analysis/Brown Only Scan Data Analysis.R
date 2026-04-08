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
# SUMMARY
    # Experimental playback still significantly reduced foraging overall. 
    # Looking was higher overall in experimental treatment with no evidence that this effect varied across periods though.

# Scan Data 
gdata <- read.csv("data/finalscan.csv")
gdata <- gdata %>% select(-Ethogram, -X.1) # earlier version of the data had two random columns with ethogram and a blank -X.1 labeled 
gdata <- gdata %>% filter(Trial !=2) # remove mixed
  #view(gdata)
  #n_distinct(gdata$Trial)

b_s_data <- gdata %>% filter(Species == "b")
  #view(b_s_data)

# Since we have count data, run Poisson first.

# ---------- Poisson Test ----------
p_forage <- glmmTMB(Forage ~ Period * Treatment + (1|Trial), 
                    family = poisson, 
                    data = b_s_data,
                    offset = log(Total+1))
res_p_forage <- simulateResiduals(p_forage)
testDispersion(res_p_forage) # no dispersion problem

p_look <- glmmTMB(Look ~ Period * Treatment + (1|Trial), 
                  family = poisson, 
                  data = b_s_data,
                  offset = log(Total+1))
res_p_look <- simulateResiduals(p_look)
testDispersion(res_p_look) # no dispersion problem

p_move <- glmmTMB(Move ~ Period * Treatment + (1|Trial), 
                  family = poisson, 
                  data = b_s_data,
                  offset = log(Total+1))
res_p_move <- simulateResiduals(p_move)
testDispersion(res_p_move) # no dispersion problem

# just incase, I'll still test negative binomial

# ---------- Negative Binomial Test ----------
nb_forage <- glmmTMB(Forage ~ Period * Treatment + (1|Trial), 
                     family = nbinom2, 
                     data = b_s_data,
                     offset = log(Total+1))
nb_look <- glmmTMB(Look ~ Period * Treatment + (1|Trial), 
                   family = nbinom2, 
                   data = b_s_data,
                   offset = log(Total+1))

nb_move <- glmmTMB(Move ~ Period * Treatment + (1|Trial), 
                   family = nbinom2, 
                   data = b_s_data,
                   offset = log(Total+1)) # THIS SPITS OUT ERROR (I will be using poisson anyways)
# AIC test
AIC (p_forage, nb_forage, p_look, nb_look, p_move, nb_move)
# forage is same, within 1 unit difference, we will use poisson
# look actually is notably lower AIC with nb model
#       ~9 AIC difference is too large to ignore
# move doesn't give me an AIC for nb, figures it spit error


# ---------- FINAL SCAN MODELS ----------
summary(p_forage) # TreatmentExperimental            2.11e-05 ***
summary(nb_look) # TreatmentExperimental             0.00172 **  
summary(p_move) # TreatmentExperimental              0.771     

# ---------- remove interaction and ANOVA/emmeans ----------
p_sansint_forage <- glmmTMB(Forage ~  Period + Treatment + (1|Trial), 
                            family = poisson, 
                            data = b_s_data,
                            offset = log(Total+1))
nb_sansint_look <- glmmTMB(Look ~  Period + Treatment + (1|Trial), 
                           family = nbinom2, 
                           data = b_s_data,
                           offset = log(Total+1))
p_sansint_move <- glmmTMB(Move ~  Period + Treatment + (1|Trial), 
                          family = poisson, 
                          data = b_s_data,
                          offset = log(Total+1))
summary(p_sansint_forage) # 1.84e-05 ***
summary(nb_sansint_look) # 5.12e-05 ***
summary(p_sansint_move) # 0.770    

# ANOVA
anova(p_forage, p_sansint_forage, test="LRT") # 0.005251 **
anova(nb_look, nb_sansint_look, test="LRT") # 0.8086
anova(p_move, p_sansint_move, test="LRT") # 0.06285 . change from full scan data, this isn't significant anymore

forage_emm <- emmeans(p_forage, pairwise ~ Period | Treatment)
summary(forage_emm) # 0.0486 playback-post (z = -4.202), 0.0041 playback-pre (z = -4.837)
# they foraged less during the playback!!

look_emm <- emmeans(nb_look, pairwise ~ Period | Treatment)
summary(look_emm) # 0.0255 playback-post 
#hmm... interesting? looked more duing playback with 2.596 z ration 
#the treatment experimental was sig, but the LRT wasn't sig 

move_emm <- emmeans(p_move, pairwise ~ Period | Treatment)
summary(move_emm) # nothing


# I think all I can say is that foraging decreased?
# We can say that the look more overall in experimental.

