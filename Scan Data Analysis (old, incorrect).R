getwd()
setwd("/Users/anyangyu/Desktop/CAP Code and Results")
#sink("Scan Data Analysis output.txt")

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
# NOTE: major errors here. I started a new R script to preserve this old analysis for reference----

# Scan Data used is called "Scan Sampling Data.xlsx" found in BOX
gdata <- read.csv("finalscan.csv")
#view(gdata)
gdata <- gdata %>% select(-Ethogram, -X.1) # earlier version of the data had two random columns with ethogram and a blank -X.1 labeled 
   # View(gdata)
    hist(gdata$Move) # forage and look are right-skewwed
    hist(gdata$Forage^.5)
    hist(log(gdata$Look)) 


# Checking for overdispersion - this is wrong
      #pos_f <- glm(Forage ~ Period * Treatment + (1|Trial), 
                  # family = poisson, 
                 #  data = gdata)
      #deviance(pos_f) / df.residual(pos_f) #Significant overdispersion!!!!! 5.76
    
      #pos_l <- glm(Look ~ Period * Treatment + (1|Trial), 
            # family = poisson, 
            # data = gdata)
      #deviance(pos_l) / df.residual(pos_l) #Significant overdispersion!!!!! 4.65
    
      #pos_m <- glm(Move ~ Period * Treatment + (1|Trial), 
            # family = poisson, 
           #  data = gdata)
      #deviance(pos_m) / df.residual(pos_m) #Significant overdispersion!!!!! 6.86
    

# Negative Binomial Test without Species
nb_sansSforage <- glmmTMB(Forage ~ Period * Treatment + (1|Trial), 
                     family = nbinom2, 
                     data = gdata,
                     offset = log(Total+1))

nb_sansSlook <- glmmTMB(Look ~ Period * Treatment + (1|Trial), 
                   family = nbinom2, 
                   data = gdata,
                   offset = log(Total+1))

nb_sansSmove <- glmmTMB(Move ~ Period * Treatment + (1|Trial), 
                   family = nbinom2, 
                   data = gdata,
                   offset = log(Total+1))

# Negative Binomial is more appropriate due to overdispersion, maybe use family = nbinom2? 
# ------------ FINAL MODEL USED WITH SPECIES ------------ 
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

AIC(nb_sansSforage, nb_forage, nb_sansSlook, nb_look, nb_sansSmove, nb_move)  
# in conclusion from AIC, WITH species is worse??? - differences are <2 and not large enough to remove Species, Species = more biologically informed model
        #   nb_sansSforage  8 619.0618
        #   nb_forage       9 621.0547
        #   nb_sansSlook    8 489.5815
        #   nb_look         9 491.3937
        #   nb_sansSmove    8 658.8018
        #   nb_move         9 660.7917
anova(nb_forage, nb_sansSforage, test="LRT") # 0.9329
anova(nb_look, nb_sansSlook, test="LRT") # 0.6648
anova(nb_move, nb_sansSmove, test="LRT") # 0.92
# Species isn't affected

summary(nb_forage) # TreatmentExperimental 1.18e-07 *
summary(nb_move) # 0.0353*
summary(nb_look) # 0.0125*
# compare the AIC of both  (Poisson and Negative Binomial)
#AIC(pos_f, nb_forage) # nb is 479 and poisson is 722 so nb is MUCH better
#AIC(pos_l, nb_look) # nb is way lower
#AIC(pos_m, nb_move) # nb is way lower

nbforage_sansint <- glmmTMB(Forage ~ Species+Treatment + Period + (1|Trial),
                          family = nbinom2,
                          data = gdata,
                          offset = log(Total+1))
nbmove_sansint <- glmmTMB(Move ~  Species+Treatment + Period+ (1|Trial),
                        family = nbinom2,
                        data = gdata,
                        offset = log(Total+1))
nblook_sansint <- glmmTMB(Look ~  Species+Treatment + Period + (1|Trial),
                        family = nbinom2,
                        data = gdata,
                        offset = log(Total+1))
summary(nbforage_sansint) # 2.1e-05 
summary(nbmove_sansint) #0.147    
summary(nblook_sansint) # 0.0108

#ANOVAs
anova(nb_forage, nbforage_sansint, test="LRT") # 0.0004371***
anova(nb_look, nblook_sansint, test="LRT") # 0.4796
anova(nb_move, nbmove_sansint, test="LRT") # 0.01938 *

forage_emm <- emmeans(nb_forage, pairwise ~ Period | Treatment)
summary(forage_emm) # 0.0009, <.0001 YES

look_emm <- emmeans(nb_look, pairwise ~ Period | Treatment)
summary(look_emm) # nope

move_emm <- emmeans(nb_move, pairwise ~ Period | Treatment)
summary(move_emm) # 0.0163 YES

# testing for residuals
res_forage <- simulateResiduals(nb_forage)
res_look <- simulateResiduals(nb_look)
res_move <- simulateResiduals(nb_move)


plot(res_forage) # no significance so good
plot(res_look) # no significance, triangles close to red line
plot(res_move) # no significance, triangles mostly on red line



# ------- testing to see if zero inflated is better? - NOPE  -------  
#ZI model causes convergence problems, could add optimizer but this over complicates model
nbzi_forage <- glmmTMB(Forage ~ Period * Treatment + (1|Trial), 
                        family = nbinom2, 
                        ziformula = ~1, 
                        data = gdata,
                       offset = log(Total+1))
nbziforage_sansint <- glmmTMB(Forage ~ Period + Treatment + (1|Trial), 
                       family = nbinom2, 
                       ziformula = ~1, 
                       data = gdata,
                       offset = log(Total+1))

nbzi_move <- glmmTMB(Move ~ Period * Treatment + (1|Trial), 
                        family = nbinom2, 
                        ziformula = ~1, 
                        data = gdata,
                     offset = log(Total+1),
                     control = glmmTMBControl(optimizer = optim, optArgs=list(method="BFGS")))
nbzimove_sansint <- glmmTMB(Move ~ Period + Treatment + (1|Trial), 
                              family = nbinom2, 
                              ziformula = ~1, 
                              data = gdata,
                            offset = log(Total+1))

nbzi_look <- glmmTMB(Look ~ Period * Treatment + (1|Trial), 
                      family = nbinom2, 
                      ziformula = ~1, 
                      data = gdata,
                     offset = log(Total+1))
nbzilook_sansint <- glmmTMB(Look ~ Period + Treatment + (1|Trial), 
                              family = nbinom2, 
                              ziformula = ~1, 
                              data = gdata,
                            offset = log(Total+1))
AIC(nb_forage, nbzi_forage) # nbzi is better??
AIC(nb_move, nbzi_move) # no differences
AIC(nb_look, nbzi_look) # minor differences

summary(nbzi_move) # forage significance only # ADDING ITTERATIONS
        #anova(nbzi_forage, nbziforage_sansint, test="LRT") # 3.927e-05
        #anova(nbzi_look, nbzilook_sansint, test="LRT") #0.1938
        #anova(nbzi_move,nbzimove_sansint, test="LRT")# ???????
        #Zforage_emm <- emmeans(nbzi_forage, pairwise ~ Period | Treatment)
        #summary(forage_emm) # pre-playback <.0001, playback-post 0.0009
        #Zlook_emm <- emmeans(nbzi_look, pairwise ~ Period | Treatment)
        #summary(look_emm) # nope
        #Zmove_emm <- emmeans(nbzi_move, pairwise ~ Period | Treatment)
        #summary(move_emm) # maybe???

# don't use zero inflation cause issues

# PLOTS
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
