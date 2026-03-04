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


# Negative Binomial Test BROWN ----
nb_b_forage <- glmmTMB(Forage ~ Period * Treatment + (1|Trial), 
                          family = nbinom2, 
                          data = b_s_data,
                          offset = log(Total+1))

nb_b_look <- glmmTMB(Look ~ Period * Treatment + (1|Trial), 
                        family = nbinom2, 
                        data = b_s_data,
                        offset = log(Total+1))

nb_b_move <- glmmTMB(Move ~ Period * Treatment + (1|Trial), 
                        family = nbinom2, 
                        data = b_s_data,
                        offset = log(Total+1))
summary(nb_b_forage) # treatment experimental 0.000217 ***
summary(nb_b_look) # treatment experimental 0.00172 **
summary(nb_b_move) # nothin, fails

# move data fails using negative binomial

# Checking for overdispersion
pos_bf <- glm(Forage ~ Period * Treatment + (1|Trial), 
 family = poisson, 
  data = b_s_data)
deviance(pos_bf) / df.residual(pos_bf) #Significant overdispersion!!!!! 5.76

pos_bl <- glm(Look ~ Period * Treatment + (1|Trial), 
 family = poisson, 
 data = b_s_data)
deviance(pos_bl) / df.residual(pos_bl) #Significant overdispersion!!!!! 4.65

pos_bm <- glm(Move ~ Period * Treatment + (1|Trial), 
 family = poisson, 
  data = b_s_data)
deviance(pos_bm) / df.residual(pos_bm) #Significant overdispersion!!!!! 6.86

# changing from nb to poisson (pasted) ----

pois_b_forage <- glmmTMB(Forage ~ Period * Treatment + (1|Trial),
                         family = poisson,
                         data = b_s_data,
                         offset = log(Total + 1))

testDispersion(simulateResiduals(pois_b_forage))
pois_b_look <- glmmTMB(Look ~ Period * Treatment + (1|Trial),
                         family = poisson,
                         data = b_s_data,
                         offset = log(Total + 1))

testDispersion(simulateResiduals(pois_b_look))
pois_b_move <- glmmTMB(Move ~ Period * Treatment + (1|Trial),
                       family = poisson,
                       data = b_s_data,
                       offset = log(Total + 1))
AIC(pois_b_forage,nb_b_forage)
AIC(pois_b_look,nb_b_look)
AIC(pois_b_move,nb_b_move)
res_pois <- simulateResiduals(pois_b_move)
testDispersion(res_pois)

# without interation 
nb_b_forage_sansint <- glmmTMB(Forage ~ Period + Treatment + (1|Trial), 
                       family = nbinom2, 
                       data = b_s_data,
                       offset = log(Total+1))

nb_b_look_sansint <- glmmTMB(Look ~ Period + Treatment + (1|Trial), 
                     family = nbinom2, 
                     data = b_s_data,
                     offset = log(Total+1))

nb_b_move_sansint <- glmmTMB(Move ~ Period + Treatment + (1|Trial), 
                     family = nbinom2, 
                     data = b_s_data,
                     offset = log(Total+1))
summary(nb_b_forage_sansint) # TE 0.000398 ***
summary(nb_b_look_sansint) # 5.12e-05 ***
summary(nb_b_move_sansint) # nothin

# ANOVAs ----
anova(nb_b_forage, nb_b_forage_sansint, test="LRT") # 0.02482 *
anova(nb_b_look, nb_b_look_sansint, test="LRT") # 0.8086
anova(nb_b_move, nb_b_move_sansint, test="LRT") # 0.01938 *

forage_emm <- emmeans(nb_b_forage, pairwise ~ Period | Treatment)
summary(forage_emm) # 0.0009, <.0001 YES

look_emm <- emmeans(nb_b_look, pairwise ~ Period | Treatment)
summary(look_emm) # nope

move_emm <- emmeans(nb_b_move, pairwise ~ Period | Treatment)
