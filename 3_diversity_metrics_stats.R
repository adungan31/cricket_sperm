library(nlme)
library(emmeans)
library(dplyr)
library(ggplot2)
library(multcomp)
library(multcompView)

###Stats###

##Need to first run "import data," "decontam,", and first half of "diversity_metrics" 
#R scripts to get divMeta file


#Observed ASVs_ANOVA
Obs <- lme(Observed ~ type*HighORLow, random =  ~ 1|maleID,
           data = divMeta, na.action = na.omit)
summary(Obs)
anova(Obs) 

#               numDF denDF  F-value p-value
#(Intercept)        1    79 971.8192  <.0001
#type               4    79   7.9531  <.0001
#HighORLow          1    28   0.4027  0.5309
#type:HighORLow     4    79   1.0173  0.4037

# Obtain least-squares means and confidence intervals
lsm <- emmeans(Obs, "type", adjust = "sidak")

# Get Tukey's letters
lsm_tukey <- cld(lsm, alpha = 0.05)

pairs(emmeans(Obs, "type", adjust="sidak")) #output added to Table S1

#Simpsons_ANOVA
Sim <- lme(Simpson ~ type*HighORLow, random =  ~ 1|maleID,
          data = divMeta, na.action = na.omit)
summary(Sim)
anova(Sim)
#               numDF denDF   F-value p-value
#(Intercept)        1    79 1208.0857  <.0001
#type               4    79    7.7792  <.0001
#HighORLow          1    28    0.2292  0.6358
#type:HighORLow     4    79    0.0621  0.9927

# Obtain least-squares means and confidence intervals
lsm <- emmeans(Sim, "type", adjust = "sidak")

# Get Tukey's letters
lsm_tukey <- cld(lsm, alpha = 0.05)

pairs(emmeans(Sim, "type", adjust="sidak")) #output added to Table S1


#Shannon_ANOVA
Shan <- lme(Shannon ~ type*HighORLow, random =  ~ 1|maleID,
           data = divMeta, na.action = na.omit)
summary(Shan)
anova(Shan)

#               numDF denDF  F-value p-value
#(Intercept)        1    79 954.3981  <.0001
#type               4    79   8.1757  <.0001
#HighORLow          1    28   0.2184  0.6439
#type:HighORLow     4    79   0.0343  0.9977

# Obtain least-squares means and confidence intervals
lsm <- emmeans(Shan, "type", adjust = "sidak")

# Get Tukey's letters
lsm_tukey <- cld(lsm, alpha = 0.05)

pairs(emmeans(Shan, "type", adjust="sidak")) #output added to Table S1

#Chao1_ANOVA
Chao <- lme(Chao1 ~ type*HighORLow, random =  ~ 1|maleID,
            data = divMeta, na.action = na.omit)
summary(Chao)
anova(Chao)

#               numDF denDF  F-value p-value
#(Intercept)        1    79 915.4182  <.0001
#type               4    79   7.8427  <.0001
#HighORLow          1    28   0.4104  0.5270
#type:HighORLow     4    79   1.2969  0.2785


# Obtain least-squares means and confidence intervals
lsm <- emmeans(Chao, "type", adjust = "sidak")

# Get Tukey's letters
lsm_tukey <- cld(lsm, alpha = 0.05)

pairs(emmeans(Chao, "type", adjust="sidak")) #output added to Table S1


rm(Obs)
rm(Shan)
rm(Sim)
rm(Chao)
rm(divMeta)
rm(lsm)
rm(lsm_tukey)
