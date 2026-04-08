library("glmmTMB")
library("lme4")
library("lmerTest")
library("performance")
library("MuMIn")
library("DHARMa")
library("see")
library("patchwork")
library("ggpattern")
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

rm(list = ls())

# Scan Data
gdata <- read.csv("data/finalscan.csv")
#view(gdata)
gdata <- gdata %>% select(-Ethogram, -X.1) # earlier version of the data had two random columns with ethogram and a blank -X.1 labeled 

# We determined that species has an effect for Look and Move, but not Forage.
#       Forage can remain with Species as additive.
#       Look ~ Species * Period * Treatment + (1|Trial)
#       Move ~ Species * Period * Treatment + (1|Trial)
# Poisson works best for Forage and Move.
# Neg Binom is best for Look due to significant AIC improvement from Poisson.

# FINAL MODELS ----
# ---------- Forage and Move Poisson ----------
p_forage <- glmmTMB(Forage ~ Species + Period * Treatment + (1|Trial), 
                          family = poisson, 
                          data = gdata,
                          offset = log(Total+1))

p_move <- glmmTMB(Move ~ Species * Period * Treatment + (1|Trial), 
                        family = poisson, 
                        data = gdata,
                        offset = log(Total+1))

# ---------- Look Negative Binomial Test ----------
nb_look <- glmmTMB(Look ~ Species * Period * Treatment + (1|Trial), 
                   family = nbinom2, 
                   data = gdata,
                   offset = log(Total+1))

# ---------- FINAL SCAN MODELS WITH Species ----------
summary(p_forage) # TreatmentExperimental            z -5.709  p 1.14e-08 ***
summary(nb_look)  # TreatmentExperimental              z 3.352   p 0.000802 ***  
summary(p_move)   # TreatmentExperimental             p 0.7644, HOWEVER 
                  #       Speciesw:TreatmentExperimental  p 0.0266 *
        

# ---------- remove interaction and ANOVA/emmeans ----------
p_sansint_forage <- glmmTMB(Forage ~ Species + Period + Treatment + (1|Trial), 
                    family = poisson, 
                    data = gdata,
                    offset = log(Total+1))
nb_sansint_look <- glmmTMB(Look ~ Species * Period + Treatment + (1|Trial), 
                   family = nbinom2, 
                   data = gdata,
                   offset = log(Total+1))
p_sansint_move <- glmmTMB(Move ~ Species * Period + Treatment + (1|Trial), 
                  family = poisson, 
                  data = gdata,
                  offset = log(Total+1))
summary(p_sansint_forage) # 7.56e-07 ***
summary(nb_sansint_look) # 0.00740 ** 
summary(p_sansint_move) # 0.2352    

# ANOVA
anova(p_forage, p_sansint_forage, test="LRT") # 6.743e-05 ***
anova(nb_look, nb_sansint_look, test="LRT") # 0.0656  
anova(p_move, p_sansint_move, test="LRT") # 0.02185 *

forage_emm <- emmeans(p_forage, pairwise ~ Period | Treatment)
summary(forage_emm) # <0.0001 playback-post (z = -4.202), 0.0001 playback-pre (z = -4.837)
# they foraged less during the playback!!

look_emm <- emmeans(nb_look, pairwise ~ Period | Treatment * Species)
summary(look_emm) 
# species brown
#      Playback - Post 0.0157 p value, 2.766 z ratio
# conclusion: Brown Capuchins looked more during the playback period than the post-playback period.

move_emm <- emmeans(p_move, pairwise ~ Period | Treatment * Species)
summary(move_emm) 
# species white
#      Playback - Pre 0.0038 p value, 3.206 z ratio
# conclusion: White-fronted Capuchins moved more during the playback period than the pre-playback period.

# Final interpretation
# Forage: playback reduced foraging across both species.

# Look: Brown showed increase looking during experimental playback compared to post
#       whereas white showed no significant change.

# Move: White showed increase movement during experimental playback compared to pre
#       whereas brown showed no significant change.

# ---------- PLOTS ---------- 

gdata$Species <- factor(gdata$Species,
                        levels = c("b","w"),
                        labels = c("S. apella","C. albifrons"))
gdata$Period <- factor(gdata$Period,
                       levels = c("Pre", "Playback", "Post"))
gdata$Treatment <- factor(gdata$Treatment,
                          levels = c("Control", "Experimental"))

forage_plot <- ggplot(gdata, aes(x = Period, y = Forage, fill = Treatment)) + 
  geom_boxplot(width = 0.7, outlier.size=1.5) + 
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))+
  labs( title = "Forage",
    x = "Period", 
    y = "# of Forage Events") + theme_classic(base_family = "Times", base_size = 12) +
  scale_fill_manual(values = c("Control" = "lightblue", "Experimental" = "orange"))
forage_plot <- forage_plot + 
  plot_annotation(title = "Scan Sampling Distributions: Forage",
                  theme = theme(plot.title=element_text(family = "Times",size = 15)))

print(forage_plot)
ggsave(
  "figures/scan_forage.pdf",
  forage_plot,
  width = 7,
  height = 5,
  units = "in",
  dpi = 600
)

look_plot <- ggplot(gdata, aes(x = Period, y = Look, fill = Treatment)) + 
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))+
  geom_boxplot(width = 0.7, outlier.size=1.5, size= 0.4) + 
  facet_wrap(~Species, nrow=1) +
  labs(title = "Look", 
       x = "Period", 
       y = "# of Look Events") + theme_classic(base_family = "Times", base_size = 12) +
  scale_fill_manual(values = c("Control" = "lightblue", "Experimental" = "orange")) +
  theme(strip.text = element_text(face = "italic"))


move_plot <- ggplot(gdata, aes(x = Period, y = Move, fill = Treatment)) + 
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  geom_boxplot(width = 0.7, outlier.size=1.5, size= 0.4) + 
  facet_wrap(~Species, nrow=1) +
  labs(title = "Move", 
       x = "Period",
       y = "# of Move Events") + theme_classic(base_family = "Times", base_size = 12) +
  scale_fill_manual(values = c("Control" = "lightblue", "Experimental" = "orange")) +
  guides(fill = "none") +
  theme(strip.text = element_text(face = "italic"))
 

combined_look_move <- look_plot / move_plot +
  plot_annotation( title = "Scan Sampling Distributions: Species-specific Look and Move",
                   theme = theme(plot.title=element_text(family = "Times",size = 15)))
ggsave(
  "figures/scan_look_move.pdf",
  combined_look_move,
  width = 7,
  height = 8,
  units = "in",
  dpi = 600
)
print(combined_look_move )
