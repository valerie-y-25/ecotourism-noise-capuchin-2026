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
#view(gdata)
gdata <- gdata %>% select(-Ethogram, -X.1) # earlier version of the data had two random columns with ethogram and a blank -X.1 labeled 
# View(gdata)
hist(gdata$Move) # forage and look are right-skewwed
hist(gdata$Forage^.5)
hist(log(gdata$Look)) 

# Since we have count data, run Poisson first.

# ---------- Poisson Test WITH Species ----------
p_forage <- glmmTMB(Forage ~ Species + Period * Treatment + (1|Trial), 
                          family = poisson, 
                          data = gdata,
                          offset = log(Total+1))
res_p_forage <- simulateResiduals(p_forage)
testDispersion(res_p_forage) # no dispersion problem

p_look <- glmmTMB(Look ~ Species + Period * Treatment + (1|Trial), 
                        family = poisson, 
                        data = gdata,
                        offset = log(Total+1))
res_p_look <- simulateResiduals(p_look)
testDispersion(res_p_look) # slight overdispersion, not significant

p_move <- glmmTMB(Move ~ Species + Period * Treatment + (1|Trial), 
                        family = poisson, 
                        data = gdata,
                        offset = log(Total+1))
res_p_move <- simulateResiduals(p_move)
testDispersion(res_p_move) # no dispersion problem

# just incase, I'll still test negative binomial

# ---------- Negative Binomial Test WITH Species ----------
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
                   offset = log(Total+1)) # THIS SPITS OUT ERROR (I will be using poisson anyways)
# AIC test
AIC (p_forage, nb_forage, p_look, nb_look, p_move, nb_move)
# forage is same, within 1 unit difference, we will use poisson
# look actually is notably lower AIC with nb model (p_look 514.2312 VS nb_look 491.3937)
#       we will use nb to account for slight overdispersion even though DHARMa didn't flag it.
#       22.84 AIC difference is too large to ignore
# move doesn't give me an AIC for nb, figures it spit error



# ---------- test to see if Species affects anything ---------- 
p_sansS_forage <- glmmTMB(Forage ~ Period * Treatment + (1|Trial), 
                    family = poisson, 
                    data = gdata,
                    offset = log(Total+1))

nb_sansSlook <- glmmTMB(Look ~ Period * Treatment + (1|Trial), 
                        family = nbinom2, 
                        data = gdata,
                        offset = log(Total+1))

p_sansS_move <- glmmTMB(Move ~ Period * Treatment + (1|Trial), 
                  family = poisson, 
                  data = gdata,
                  offset = log(Total+1))
AIC (p_forage, p_sansS_forage, nb_look, nb_sansSlook, p_move, p_sansS_move)
# all < 2 AIC units, models essentially equivalent, Species is biologically relevant so retained.

# ---------- FINAL SCAN MODELS WITH Species ----------
summary(p_forage) # TreatmentExperimental            -0.829754   0.145350  -5.709 1.14e-08 ***
summary(nb_look) #TreatmentExperimental             0.68080    0.27259   2.498   0.0125 *  
summary(p_move) #TreatmentExperimental             0.183889   0.087374   2.105   0.0353 * 

# ---------- remove interaction and ANOVA/emmeans ----------
p_sansint_forage <- glmmTMB(Forage ~ Species + Period + Treatment + (1|Trial), 
                    family = poisson, 
                    data = gdata,
                    offset = log(Total+1))
nb_sansint_look <- glmmTMB(Look ~ Species + Period + Treatment + (1|Trial), 
                   family = nbinom2, 
                   data = gdata,
                   offset = log(Total+1))
p_sansint_move <- glmmTMB(Move ~ Species + Period + Treatment + (1|Trial), 
                  family = poisson, 
                  data = gdata,
                  offset = log(Total+1))
summary(p_sansint_forage) # 7.56e-07 ***
summary(nb_sansint_look) # 0.0108 *
summary(p_sansint_move) # n/a

# ANOVA
anova(p_forage, p_sansint_forage, test="LRT") # 6.743e-05 ***
anova(nb_look, nb_sansint_look, test="LRT") # 0.4796
anova(p_move, p_sansint_move, test="LRT") # 0.0194 *

forage_emm <- emmeans(p_forage, pairwise ~ Period | Treatment)
summary(forage_emm) # <0.0001 playback-post (z = -4.202), 0.0001 playback-pre (z = -4.837)
# they foraged less during the playback!!

look_emm <- emmeans(nb_look, pairwise ~ Period | Treatment)
summary(look_emm) # nope

move_emm <- emmeans(p_move, pairwise ~ Period | Treatment)
summary(move_emm) # 0.0166 playback-pre (z = 2.746)
# they moved more during the playback!


# ---------- PLOTS ---------- 


ggplot(gdata, aes(x = Period, y = Forage, fill = Treatment)) + geom_boxplot() + 
  labs(title = "Forage Distribution by Period and Treatment", x = "Period", y = "Forage") + theme_minimal() +
  scale_fill_manual(values = c("Control" = "lightblue", "Experimental" = "orange")) +
  scale_x_discrete(limits = c("Pre", "Playback", "Post")) 


ggplot(gdata, aes(x = Period, y = Look, fill = Treatment)) + geom_boxplot() +
  labs(title = "Look Distribution by Period and Treatment", x = "Period", y = "Look") + theme_minimal() +
  scale_fill_manual(values = c("Control" = "lightblue", "Experimental" = "orange")) +
  scale_x_discrete(limits = c("Pre", "Playback", "Post")) 


ggplot(gdata, aes(x = Period, y = Move, fill = Treatment)) + geom_boxplot() +
  labs(title = "Move Distribution by Period and Treatment", x = "Period",y = "Move") +
  theme_minimal() +
  scale_fill_manual(values = c("Control" = "lightblue", "Experimental" = "orange")) +
  scale_x_discrete(limits = c("Pre", "Playback", "Post")) 


# Adding the plots together
p1 <- ggplot(gdata, aes(x = Period, y = Forage, fill = Treatment)) + 
  geom_boxplot() + 
  labs(title = "Forage", x = "Period", y = "Forage") + 
  theme_minimal() +
  scale_fill_manual(values = c("Control" = "lightblue", "Experimental" = "orange"), guide = "none") +  # Remove legend
  scale_x_discrete(limits = c("Pre", "Playback", "Post")) 

p2 <- ggplot(gdata, aes(x = Period, y = Look, fill = Treatment)) + 
  geom_boxplot() + 
  labs(title = "Look", x = "Period", y = "Look") + 
  theme_minimal() +
  scale_fill_manual(values = c("Control" = "lightblue", "Experimental" = "orange"), guide = "none") +  # Remove legend
  scale_x_discrete(limits = c("Pre", "Playback", "Post")) 

p3 <- ggplot(gdata, aes(x = Period, y = Move, fill = Treatment)) + 
  geom_boxplot() + 
  labs(title = "Move", x = "Period", y = "Move") +
  theme_minimal() +
  scale_fill_manual(values = c("Control" = "lightblue", "Experimental" = "orange")) + # Keep legend
  scale_x_discrete(limits = c("Pre", "Playback", "Post")) +
  theme(legend.position = "right")  # Move legend to the right of "Move" plot

combined_plot <- (p1 | p2 | p3) + 
  plot_annotation(title = "Scan Sampling Distribution by Period and Treatment")

print(combined_plot)

#sink()

