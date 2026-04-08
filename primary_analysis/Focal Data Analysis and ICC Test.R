#sink("Focal Data Analysis and ICC Test output.txt")
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

rm(list = ls())

# Focal Data used is called "Focal Sampling Data - combined BORIS data.xlsx" found in BOX
data <- read.csv("data/finalfocal.csv")
      #View(data)
      #colnames(data)
    #hist(data$Duration, breaks = 30, main = "Histogram of Duration", xlab = "Duration", col = "skyblue")
      # duration is right-skewed, data set is not normal!!

# Here, I tried to use sqrt to normalize but doesn't work, neither does log.
    #data$sqrtduration <- sqrt(data$Duration)
    #hist(data$sqrtduration, breaks = 30, main = "Histogram of sqrt(Duration)", xlab = "sqrt(Duration)", col = "skyblue")
 


# Finding total time in sight since it varies, used for offset() later.
  data <- data %>%
      group_by(Trial, Species, Treatment, Period) %>%
      mutate(total_time_in_sight = sum(Duration[Behavior %in% c("move", "forage", "look", "other")])) %>% 
      ungroup()
      # here, I am grouping the data set so that the total time in sight summation applies to each trial, 
      # then, mutate func creates a new column in the dataset called total_time_in_sight which is the sum of all known behaviors
      # we could have used out of sight for this but katie said they do it this way, plus out of sight isn't perfectly correct 

# This subsets data for each behavior.
    forage_data <- data %>% filter(Behavior == "forage")
    # hist(look_data$Duration, breaks = 30, main = "Histogram of Duration", xlab = "Duration", col = "green")
      # use View(whatever dataset here) to see what it looks like or check modifications along the way
      # use hist(dataset$Duration) to view data in a graph, move seems somewhat normal but look and forage are right-skewed
    look_data <- data %>% filter(Behavior == "look")
    move_data <- data %>% filter(Behavior == "move")

# Testing NEGATIVE BINOMAL model first (DOES NOT WORK)
#forage_nbmodel <- glmer.nb(Duration ~ Period + Treatment + Period*Treatment + (1|Trial) , 
                      #data = forage_data) #offset?
#summary(forage_nbmodel) # displays table of regression coefficients	
#testDispersion(forage_nbmodel) # tests for under- and overdispersion 
  # significant p-value < 2.2e-16 indicate poor model fit, VERY small p-value

#simulationOutput <- simulateResiduals(fittedModel = forage_nbmodel, plot = F)
#plot(simulationOutput) # generates QQ plot and other output to check for deviations from model assumptions
#testZeroInflation(simulationOutput) # shows that neg binomial doesn't fit



# Testing zero-inflated GAMMA model because we have a lot of zeros and the zeros are important
# LACKING SPECIES
forage_model <- glmmTMB(Duration ~ Treatment + Period+ Period*Treatment + (1|Trial),
                        family = ziGamma(link = "log"),
                        ziformula = ~ 1,
                        data = forage_data,
                        offset = log(total_time_in_sight+1e-6)) # plus 1 INSTEAD 
look_model <- glmmTMB(Duration ~  Treatment + Period+ Period*Treatment + (1|Trial),
                        family = ziGamma(link = "log"),
                        ziformula = ~ 1,
                        data = look_data,
                        offset = log(total_time_in_sight+1e-6))
move_model <- glmmTMB(Duration ~  Treatment + Period+ Period*Treatment + (1|Trial),
                        family = ziGamma(link = "log"),
                        ziformula = ~ 1,
                        data = move_data,
                        offset = log(total_time_in_sight+1e-6))
# looking at Treatmentexperimental
summary(forage_model) #SIG 0.000119
summary(look_model) #not 0.127
summary(move_model) #SIG 0.00624
AIC(forage_model,look_model,move_model)
    # use ANOVA test to run with and without interaction, look at p-value
    # if close or signif. use pairwise (only if interaction ANOVA p-value is significant)
# without interaction
forage_sansint <- glmmTMB(Duration ~ Treatment + Period + (1|Trial),
                        family = ziGamma(link = "log"),
                        ziformula = ~ 1,
                        data = forage_data,
                        offset = log(total_time_in_sight+1e-6))
move_sansint <- glmmTMB(Duration ~  Treatment + Period+ (1|Trial),
                      family = ziGamma(link = "log"),
                      ziformula = ~ 1,
                      data = move_data,
                      offset = log(total_time_in_sight+1e-6))
look_sansint <- glmmTMB(Duration ~  Treatment + Period + (1|Trial),
                      family = ziGamma(link = "log"),
                      ziformula = ~ 1,
                      data = look_data,
                      offset = log(total_time_in_sight+1e-6))
summary(forage_sansint) # p=0.0166 SIG
summary(look_sansint) # p=0.269 NOT
summary(move_sansint) # p=0.0416 SIG
# anova to see if interaction means anything
anova(forage_model, forage_sansint, test="LRT") # p=0.004349 YAY 
anova(look_model, look_sansint, test="LRT") # p=0.5178 NO
anova(move_model, move_sansint, test="LRT") # p=0.02395 YAY

# pairwise because for forage and move interaction does, only FORAGE was SIGNIFICANT 
forage_emm <- emmeans(forage_model, pairwise ~ Period | Treatment)
forage_emm 
plot(forage_emm)
    # playback - pre experimental contrast p = 0.0355 !!, none of the control period contrasts are significant as expected
    # this should mean CAPUCHINS FORAGE SIG. LESS DURING EXPERIMENTAL PLAYBACK
look_emm <- emmeans(look_model, pairwise ~ Period | Treatment)
look_emm
plot(look_emm)
move_emm <- emmeans(move_model, pairwise ~ Period | Treatment)
move_emm
plot(move_emm) 

#   ----------------- ACTUAL FINAL MODEL USED (WITH SPECIES) -----------------  
# WITH SPECIES (IS BETTER BECAUSE IS BIOLOGICALLY INFORMED MODEL) - however species does not effect behavior (aka white and brown are stat. insig.)

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
summary(Sforage_model) # AIC differs by 0.082, better Yes
summary(Slook_model) # AIC differs by ~1, worse No
summary(Smove_model) # AIC differs by ~2, worse Yes

AIC(forage_model, Sforage_model, look_model, Slook_model, move_model, Smove_model)
# move use without species

#   ----------------- ANOVA AND PAIRWISED (WITH SPECIES) -----------------  

Sforage_sansint <- glmmTMB(Duration ~ Species + Treatment + Period + (1|Trial),
                          family = ziGamma(link = "log"),
                          ziformula = ~ 1,
                          data = forage_data,
                          offset = log(total_time_in_sight+1e-6))
Smove_sansint <- glmmTMB(Duration ~  Species + Treatment + Period+ (1|Trial),
                        family = ziGamma(link = "log"),
                        ziformula = ~ 1,
                        data = move_data,
                        offset = log(total_time_in_sight+1e-6))
Slook_sansint <- glmmTMB(Duration ~  Species + Treatment + Period + (1|Trial),
                        family = ziGamma(link = "log"),
                        ziformula = ~ 1,
                        data = look_data,
                        offset = log(total_time_in_sight+1e-6))
summary(Sforage_sansint) # p=0.00884
summary(Slook_sansint) # p=0.236
summary(Smove_sansint) # p=0.0455

anova(Sforage_model, Sforage_sansint, test="LRT") # p=0.002917 YAY 
anova(Slook_model, Slook_sansint, test="LRT") # p=0.5169 eh
anova(Smove_model, Smove_sansint, test="LRT") # p=0.0355 YAY

# pairwise, only FORAGE was SIGNIFICANT   0.0244
Sforage_emm <- emmeans(Sforage_model, pairwise ~ Period | Treatment)
Sforage_emm 
plot(Sforage_emm)
    # playback - pre experimental contrast p = 0.0244 !!, none of the control period contrasts are significant as expected
    # this should mean CAPUCHINS FORAGE SIG. LESS DURING EXPERIMENTAL PLAYBACK
Slook_emm <- emmeans(Slook_model, pairwise ~ Period | Treatment)
Slook_emm
plot(Slook_emm)
Smove_emm <- emmeans(Smove_model, pairwise ~ Period | Treatment)
Smove_emm
plot(Smove_emm)

# check models for fit, FIRST ONES without species 
simulation_forage <- simulateResiduals(forage_model)
plot(simulation_forage) # not all good in QQ plot
simulation_look <- simulateResiduals(look_model)
plot(simulation_look) # all good in QQ 
simulation_move <- simulateResiduals(move_model)
plot(simulation_move) # DEVIATION SIGNIFICANT

# SECOND ONES WITH SPECIES
simulation_Sforage <- simulateResiduals(Sforage_model)
plot(simulation_Sforage) # not good.
simulation_Slook <- simulateResiduals(Slook_model)
plot(simulation_Slook) # all good in QQ
simulation_Smove <- simulateResiduals(Smove_model)
plot(simulation_Smove) # all good in QQ
  # all models state combined adjusted quantile test significant so...

# PLOT: filter out "out of sight" and "other" behaviors 
filtered_data <- data %>%
  filter(Behavior %in% c("forage", "look", "move"))
# everything is lowercase --> uppercase labels
filtered_data$Behavior <- factor(
  filtered_data$Behavior,
  levels =c("forage", "look", "move"),
  labels = c("Forage", "Look", "Move"),
)

filtered_data$Period <- factor(filtered_data$Period, 
                               levels = c("pre", "playback", "post"),
                               labels = c("Pre", "Playback", "Post"))
filtered_data$Treatment <- factor(filtered_data$Treatment, 
                               levels = c("control", "experimental"),
                               labels = c("Control", "Experimental"))

# create the box-and-whisker plot with proper ordering of Period
focal_plot <- ggplot(filtered_data, aes(x = Period, y = Duration, fill = Treatment)) +
  geom_boxplot(width = 0.7, outlier.size=1.5, size= 0.4) +
  scale_y_continuous(expand = expansion(mult = c(0.05,0.08)))+
  facet_wrap(~Behavior, scales = "free_y", nrow = 1) +
  labs(
       x = "Period",
       y = "Duration (s)") +
  theme_classic(base_family = "Times", base_size = 12) +
  scale_fill_manual(values = c("Control" = "lightblue", "Experimental" = "orange")) +
  theme(
    strip.background = element_blank())  +
  plot_annotation(
    title = "Focal Sampling Distributions", 
    theme = theme(plot.title=element_text(family = "Times",size = 15)))
print(focal_plot)

ggsave(
  "figures/focal_plot.pdf",
  focal_plot,
  width = 8,
  height = 4,
  units = "in",
  dpi = 600
)

#   ----------------- other random tests -----------------  


# here I'm just seeing if the model improves if we add different things, all seem to not help + don't effect duration so
# testing OBSERVER added to species model - WORSE
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

# testing WEATHER added to species model - WORSE
WSforage_model <- glmmTMB(Duration ~ Weather+Treatment * Period + Species + (1|Trial),
                          family = ziGamma(link = "log"),
                          ziformula = ~ 1,
                          data = forage_data,
                          offset = log(total_time_in_sight+1e-6))

WSlook_model <- glmmTMB(Duration ~ Weather+Treatment * Period + Species + (1|Trial),
                        family = ziGamma(link = "log"),
                        ziformula = ~ 1,
                        data = look_data,
                        offset = log(total_time_in_sight+1e-6))

WSmove_model <- glmmTMB(Duration ~ Weather+Treatment * Period + Species + (1|Trial),
                        family = ziGamma(link = "log"),
                        ziformula = ~ 1,
                        data = move_data,
                        offset = log(total_time_in_sight+1e-6))
summary(WSforage_model) # worse
summary(WSlook_model) # worse
summary(WSmove_model) # basically the same but lik .073 higher so worse
AIC(WSforage_model, Sforage_model, WSlook_model, Slook_model, WSmove_model, Smove_model)
# also no sig p-values for weather

# testing DISTANCE added to species model - error in output despite giving me a summary???
DSforage_model <- glmmTMB(Duration ~ Distance+Treatment * Period + Species + (1|Trial),
                          family = ziGamma(link = "log"),
                          ziformula = ~ 1,
                          data = forage_data,
                          offset = log(total_time_in_sight+1e-6))

DSlook_model <- glmmTMB(Duration ~ Distance+Treatment * Period + Species + (1|Trial),
                        family = ziGamma(link = "log"),
                        ziformula = ~ 1,
                        data = look_data,
                        offset = log(total_time_in_sight+1e-6))

DSmove_model <- glmmTMB(Duration ~ Distance+Treatment * Period + Species + (1|Trial),
                        family = ziGamma(link = "log"),
                        ziformula = ~ 1,
                        data = move_data,
                        offset = log(total_time_in_sight+1e-6))
summary(DSforage_model) # worse
summary(DSlook_model) # better by ~1
summary(DSmove_model) # worse
AIC(DSforage_model, Sforage_model, DSlook_model, Slook_model, DSmove_model, Smove_model)
# also no sig p-values for distance

# I don't think anything below is usable or means anything unfortunately.
# testing Total_indiv added to species model 
# covariate
SSforage_model <- glmmTMB(Duration ~ Total_indiv+Treatment * Period + Species + (1|Trial),
                          family = ziGamma(link = "log"),
                          ziformula = ~ 1,
                          data = forage_data,
                          offset = log(total_time_in_sight+1e-6))

SSlook_model <- glmmTMB(Duration ~ Total_indiv+Treatment * Period + Species + (1|Trial),
                        family = ziGamma(link = "log"),
                        ziformula = ~ 1,
                        data = look_data,
                        offset = log(total_time_in_sight+1e-6))

SSmove_model <- glmmTMB(Duration ~ Total_indiv+Treatment * Period + Species + (1|Trial),
                        family = ziGamma(link = "log"),
                        ziformula = ~ 1,
                        data = move_data,
                        offset = log(total_time_in_sight+1e-6))
SSmove_model_sansint <- glmmTMB(Duration ~ Total_indiv+Treatment + Period + Species + (1|Trial),
                        family = ziGamma(link = "log"),
                        ziformula = ~ 1,
                        data = move_data,
                        offset = log(total_time_in_sight+1e-6))
summary(SSforage_model) # none
summary(SSlook_model) # none
summary(SSmove_model) # Total_indiv p=0.00676 **??? z=2.708
anova(SSmove_model,SSmove_model_sansint, test="LRT") #0.02983
SSmove_emm <- emmeans(SSmove_model_sansint, pairwise ~ Period | Treatment)
SSmove_emm # nothing sig  

AIC(SSforage_model, Sforage_model, SSlook_model, Slook_model, SSmove_model, Smove_model)
# about same for forage and look, but move differs 
#   model        df   AIC values
# SSforage_model 11   1221.368 <-- with Total_indiv is slightly lower
# Sforage_model  10   1222.902
# SSlook_model   11   1014.219 <-- with Total_indiv is slightly higher
# Slook_model    10   1013.759
# SSmove_model   11   1322.043 <-- with Total_indiv is 6.864 AIC units
# Smove_model    10   1328.907

cor(forage_data$Total_indiv, forage_data$Duration)  # -0.1157458
cor(look_data$Total_indiv, look_data$Duration) # -0.06067175
cor(move_data$Total_indiv, move_data$Duration) # 0.3646775

#sink()
