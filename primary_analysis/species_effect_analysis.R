library("glmmTMB")
library("lme4")
library("lmerTest")
library("performance")
library("MuMIn")
library("DHARMa")
library("dplyr")

rm(list = ls())
# ----- Main Conclusion ----
# Focal data can remain the same with model as 
#     Duration ~ Treatment * Period + Species + (1|Trial)

# Scan data must change for Look and Move. Forage can remain.
#     Look ~ Species * Period * Treatment + (1|Trial)
#     Move ~ Species * Period * Treatment + (1|Trial)

# ----- FOCAL DATA -----
# No evidence that species modifies the effect of playback in focal.
data <- read.csv("data/finalfocal.csv")
data <- data %>%
  group_by(Trial, Species, Treatment, Period) %>%
  mutate(total_time_in_sight = sum(Duration[Behavior %in% c("move", "forage", "look", "other")])) %>% 
  ungroup()
forage_data <- data %>% filter(Behavior == "forage")
look_data <- data %>% filter(Behavior == "look")
move_data <- data %>% filter(Behavior == "move")
# Forage 
Sforage_model <- glmmTMB(Duration ~ Treatment * Period + Species + (1|Trial),
                         family = ziGamma(link = "log"),
                         ziformula = ~ 1,
                         data = forage_data,
                         offset = log(total_time_in_sight+1e-6))
SIforage_model <- glmmTMB(Duration ~ Treatment * Period * Species + (1|Trial),
                         family = ziGamma(link = "log"),
                         ziformula = ~ 1,
                         data = forage_data,
                         offset = log(total_time_in_sight+1e-6))
anova(Sforage_model, SIforage_model, test="LRT")

# Look
Slook_model <- glmmTMB(Duration ~ Treatment * Period + Species + (1|Trial),
                       family = ziGamma(link = "log"),
                       ziformula = ~ 1,
                       data = look_data,
                       offset = log(total_time_in_sight+1e-6))
SIlook_model <- glmmTMB(Duration ~ Treatment * Period * Species + (1|Trial),
                       family = ziGamma(link = "log"),
                       ziformula = ~ 1,
                       data = look_data,
                       offset = log(total_time_in_sight+1e-6))
anova(Slook_model, SIlook_model, test="LRT")

# Move
Smove_model <- glmmTMB(Duration ~ Treatment * Period + Species + (1|Trial),
                       family = ziGamma(link = "log"),
                       ziformula = ~ 1,
                       data = move_data,
                       offset = log(total_time_in_sight+1e-6))
SImove_model <- glmmTMB(Duration ~ Treatment * Period * Species + (1|Trial),
                       family = ziGamma(link = "log"),
                       ziformula = ~ 1,
                       data = move_data,
                       offset = log(total_time_in_sight+1e-6))
anova(Smove_model, SImove_model, test="LRT")


# ----- SCAN DATA -----
# Some evidence that species modifies the effect of playback in focal.
# Forage remains same. Look and Move indicate species difference.
gdata <- read.csv("data/finalscan.csv") %>% 
  select(-Ethogram, -X.1)
# Forage 
p_forage <- glmmTMB(Forage ~ Species + Period * Treatment + (1|Trial), 
                    family = poisson, 
                    data = gdata,
                    offset = log(Total+1))
pI_forage <- glmmTMB(Forage ~ Species * Period * Treatment + (1|Trial), 
                    family = poisson, 
                    data = gdata,
                    offset = log(Total+1))
anova(p_forage, pI_forage, test="LRT")

# Look
nb_look <- glmmTMB(Look ~ Species + Period * Treatment + (1|Trial), 
                   family = nbinom2, 
                   data = gdata,
                   offset = log(Total+1))
nbI_look <- glmmTMB(Look ~ Species * Period * Treatment + (1|Trial), 
                   family = nbinom2, 
                   data = gdata,
                   offset = log(Total+1))
anova(nb_look, nbI_look, test="LRT")

# Move
p_move <- glmmTMB(Move ~ Species + Period * Treatment + (1|Trial), 
                  family = poisson, 
                  data = gdata,
                  offset = log(Total+1))
pI_move <- glmmTMB(Move ~ Species * Period * Treatment + (1|Trial), 
                  family = poisson, 
                  data = gdata,
                  offset = log(Total+1))
anova(p_move, pI_move, test="LRT") 

