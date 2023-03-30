## version history -------------------------------------------------------------
# v1 by P.K.N. & B.A.P. 29 Mar 2023

## setup -----------------------------------------------------------------------
rm(list=ls(all.names=TRUE))
setwd("/insert/wd/here")
set.seed(4739)
suppressPackageStartupMessages({library(tidyverse)
  library(bda)
  library(mediation)
  library(ROCR)})
options(stringsAsFactors=FALSE)
theme_set(theme_classic())

## inputs ----------------------------------------------------------------------
outcome <- read.csv("metafile-final.csv")

## data wrangling---------------------------------------------------------------
# make outcomes binary + filter out NAs
outcome$outnum <- c( "ALIVE" = 0, "DEAD" = 1, "Dead" = 1, "DISCHARGE" = 0, 
                     "NADM" = 0)[outcome$outcome]
outcome <- outcome %>% filter(!is.na(outnum))

# make age numeric
outcome$numage = NA
counter = 0
for (pt_age in outcome$new_age){
  counter = counter + 1

  if (grepl("Y", pt_age, fixed = TRUE)) {
    #remove 'Y' --> string to numeric
    outcome$numage[counter] <- strtoi(gsub('Y', '', pt_age))
    
  } else if (grepl("M", pt_age, fixed = TRUE)) {
    #convert months to years
    outcome$numage[counter] = strtoi(gsub('M', '', pt_age))/12
    
  } else if (grepl("W", pt_age, fixed = TRUE)){
    #convert weeks to years
    outcome$numage[counter] = strtoi(gsub('W', '', pt_age))/52}}
rm(counter, pt_age)
# create df with ages and outcomes
ageout <- outcome %>% filter(!is.na(numage)) %>% dplyr::select(numage, outnum)

# make sex binary after filtering NAs
sexout <- outcome %>% filter(!is.na(sex))
sexout$sex_num <- ifelse(sexout$sex=="F", 1, 0)
sexout <- sexout %>% dplyr::select(sex_num, outnum)

# make receipt of treatment binary after filtering NAs
ribavirin <- outcome %>% filter(!is.na(q_treatment))
ribavirin$treatment <- ifelse(ribavirin$q_treatment, 1, 0)
ribavirin <- ribavirin %>% dplyr::select(treatment, outnum)

# replace data points lacking Ct value with NA
outcome$CT_L_gene[grepl("POS",outcome$CT_L_gene)]<- NA
# impute NEG values as Ct = 45 (number of PCR cycles)
outcome$CT_L_gene[outcome$CT_L_gene == "NEG"] <- 45
# convert Ct to numeric
outcome$CT_L_gene <- as.numeric(outcome$CT_L_gene)
# replace data points lacking Ct value with NA
outcome$CT_S_gene[grepl("POS",outcome$CT_S_gene)]<- NA
# impute NEG values as Ct = 45 (number of PCR cycles)
outcome$CT_S_gene[outcome$CT_S_gene == "NEG"] <- 45
# convert Ct to numeric
outcome$CT_S_gene <- as.numeric(outcome$CT_S_gene)
# calculate mean Ct (2-target PCR test)
outcome$meanct <- rowMeans(outcome[,6:7], na.rm=TRUE)
# create df with mean Ct and outcomes
ctmean <-outcome %>% filter(!if_all(c(CT_S_gene, CT_L_gene), is.na))
ctmean <- ctmean %>% dplyr::select(meanct, outnum)

# make malaria status binary after filtering NAs
malaria <- outcome %>% filter(!is.na(malaria))
malaria$malarianum <- ifelse(malaria$malaria, 1, 0)
malaria <- malaria %>% dplyr::select(malarianum, outnum)

# make HIV status binary
outcome$hivnum <- ifelse(outcome$Human_immunodeficiency_virus_1, 1, 0)
hiv <- outcome %>% dplyr::select(hivnum, outnum)

# make pegivirus status binary
outcome$peginum <- ifelse(outcome$Pegivirus_C, 1, 0)
pegivirus <- outcome %>% dplyr::select(peginum, outnum)

# make pregnancy status binary
preg <- outcome %>% filter(sex == "F") %>% 
  dplyr::select(numage, pregnant, outnum)
# assume individuals < 10 or > 60 are not pregnant
preg <- preg %>% 
  mutate(pregnantextra=ifelse(numage < 10 | numage > 60, FALSE, NA)) %>%
  mutate(pregnancy = coalesce(pregnant, pregnantextra)) %>%
  dplyr::select(pregnancy, outnum) %>%
  filter(!is.na(pregnancy))
preg$pregnancy <- ifelse(preg$pregnancy, 1, 0)

## plot variables---------------------------------------------------------------
#plot age
ggplot(ageout, aes(x = outnum,y = numage, group = outnum)) + 
  geom_boxplot() + 
  labs(x = 'Lassa Outcome', y = 'Age (Years)') +
  scale_x_continuous(breaks = c(0,1), label=c("Survived", "Deceased"))

#plot sex
ggplot(sexout, aes(outnum, fill = factor(sex_num))) + 
  geom_bar(position= "fill") + scale_fill_grey() +
  labs(x = 'Lassa Outcome', y = 'Proportion Male') +
  scale_x_continuous(breaks = c(0,1), label=c("Survived", "Deceased")) +
  ylim(0,0.5) + theme(legend.position = "none")

#plot treatment
ggplot(ribavirin, aes(outnum, fill = factor(treatment))) +                          
  geom_bar(position = "fill") + scale_fill_grey() +
  labs(x = 'Lassa Outcome', y = 'Proportion Treated') +
  scale_x_continuous(breaks = c(0,1), label=c("Survived", "Deceased")) +
  ylim(0,0.8) + theme(legend.position = "none")

#plot mean Ct
ggplot(ctmean, aes(x = outnum, y= meanct, group = outnum)) +
  geom_boxplot() +
  labs(x = 'Lassa Outcome', y = 'CT value') +
  scale_x_continuous(breaks = c(0,1), label=c("Survived", "Deceased"))

#plot malaria
ggplot(malaria, aes(outnum, fill = factor(malarianum))) +                          
  geom_bar(position = "fill") + scale_fill_grey() +
  labs(x = 'Lassa Outcome', y = 'Proportion with Malaria') +
  scale_x_continuous(breaks = c(0,1), label=c("Survived", "Deceased")) +
  ylim(0,0.8) + theme(legend.position = "none")

#plot hiv
ggplot(hiv, aes(outnum, fill = factor(hivnum))) +                          
  geom_bar(position = "fill") + scale_fill_grey() +
  labs(x = 'Lassa Outcome', y = 'Proportion with HIV-1') +
  scale_x_continuous(breaks = c(0,1), label=c("Survived", "Deceased")) +
  ylim(0,0.05) + theme(legend.position = "none")

#plot pegivirus
ggplot(pegivirus, aes(outnum, fill = factor(peginum))) +                          
  geom_bar(position = "fill") + scale_fill_grey() +
  labs(x = 'Lassa Outcome', y = 'Proportion with Pegivirus') +
  scale_x_continuous(breaks = c(0,1), label=c("Survived", "Deceased")) +
  ylim(0,0.1) + theme(legend.position = "none")

#plot pregnancy
ggplot(preg, aes(outnum, fill = factor(pregnancy))) +                          
  geom_bar(position = "fill") + scale_fill_grey() +
  labs(x = 'Lassa Outcome', y = 'Proportion of Women Pregnant') +
  scale_x_continuous(breaks = c(0,1), label=c("Survived", "Deceased")) +
  ylim(0,0.15) + theme(legend.position = "none")

## univariate analyses----------------------------------------------------------
# sex
unisex <- glm(outnum ~ sex_num, data = sexout, family = "binomial")
summary(unisex)
sexci <- exp(confint(unisex))[2,]
rm(sexout)

# treatment
unidrug <- glm(outnum ~ treatment, data = ribavirin, family = "binomial")
summary(unidrug)
drugci <- exp(confint(unidrug))[2,]
rm(ribavirin)

# age
uniage <- glm(outnum ~ numage, data = ageout, family = "binomial")
summary(uniage)
ageci <- exp(confint(uniage))[2,]
rm(ageout)

# Ct value
unictmean <- glm(outnum ~ meanct, data = ctmean, family = "binomial")
summary(unictmean)
ctci <- exp(confint(unictmean))[2,]
rm(ctmean)

# malaria
unimalaria <- glm(outnum ~ malarianum, data = malaria, family = "binomial")
summary(unimalaria)
malariaci <- exp(confint(unimalaria))[2,]
rm(malaria)

# hiv
unihiv <- glm(outnum ~ hivnum, data = hiv, family = "binomial")
summary(unihiv)
hivci <- exp(confint(unihiv))[2,]
rm(hiv)

# pegivirus
unipegi <- glm(outnum ~ peginum, data = pegivirus, family = "binomial")
summary(unipegi)
pegici <- exp(confint(unipegi))[2,]
rm(pegivirus)

# pregnancy
unipreg <- glm(outnum ~ pregnancy, data = preg, family = "binomial")
summary(unipreg)
pregci <- exp(confint(unipreg))[2,]
rm(preg)

# create a table with univariate outputs
univariate <- as.data.frame(t(cbind(as.data.frame(coef(summary(unisex))[2,]), 
                    as.data.frame(coef(summary(unidrug))[2,]),
                    as.data.frame(coef(summary(uniage))[2,]),
                    as.data.frame(coef(summary(unipreg))[2,]),
                    as.data.frame(coef(summary(unictmean))[2,]),
                    as.data.frame(coef(summary(unipegi))[2,]),
                    as.data.frame(coef(summary(unimalaria))[2,]),
                    as.data.frame(coef(summary(unihiv))[2,]))))
univariate$OR = exp(univariate$Estimate)
rm(unisex, unidrug, uniage, unipreg, unictmean, unipegi, unimalaria, unihiv)
ci <- rbind(sexci, drugci, ageci, pregci, ctci,pegici, malariaci, hivci)
univariate <- cbind(univariate, ci)
rm(ci, ageci, ctci, drugci, hivci, malariaci, pegici, pregci, sexci)
colnames(univariate) <- c("coefficient","std error", "z-value", "p", "OR",
                          "OR CI 2.5", "OR CI 97.5")

## multivariate analyses--------------------------------------------------------
# make df with univariate predictors at p < 0.25
multidf <- outcome %>% 
  dplyr::select(numage, meanct, q_treatment, Pegivirus_C, outnum) %>%
  filter(!is.na(numage)) %>% 
  filter(!is.na(meanct)) %>%
  filter(!is.na(q_treatment)) %>% 
  rename("Age" = "numage", "Ct" = "meanct", "Ribavirin" = "q_treatment",
         "Pegivirus" = "Pegivirus_C")
multidf$Ribavirin <- ifelse(multidf$Ribavirin, 1, 0)
multidf$Pegivirus <- ifelse(multidf$Pegivirus, 1, 0)

# make df with potential predictors of Ct
pegidf <- outcome %>% 
  dplyr::select(numage, meanct, Pegivirus_C, outnum) %>%
  filter(!is.na(numage)) %>% 
  filter(!is.na(meanct)) %>%
  rename("Age" = "numage", "Ct" = "meanct", "Pegivirus" = "Pegivirus_C")
pegidf$Pegivirus <- ifelse(pegidf$Pegivirus, 1, 0)

# construct linear regression model (Ct as dependent variable)
predCt <- glm(Ct ~ Pegivirus + Age, data = pegidf)
summary(predCt)
# assess residual plot
plot(predCt)

# construct logistic regression model (outcome as dependent variable)
predOutcome <- glm(outnum ~ Ct + Pegivirus + Ribavirin + Age, data = multidf, 
                   family = "binomial")
summary(predOutcome)
# assess ROC curve
pred <- prediction(predict(predOutcome, type="response"), multidf$outnum)
plot(performance(pred, "tpr", "fpr"))
auc_ROCR <- performance(pred, measure = "auc")@y.values[[1]]
rm(pred, auc_ROCR)

## mediation analyses-----------------------------------------------------------
# sobel tests
mediation.test(pegidf$Ct, pegidf$Pegivirus, pegidf$outnum)
mediation.test(pegidf$Ct, pegidf$Age, pegidf$outnum)

# bootstrapping analyses
modCt <- glm(Ct ~ Pegivirus + Age, data = multidf)
pegiboot <- mediate(modCt, predOutcome, sims = 500, boot = TRUE, treat = "Pegivirus",
                   mediator = "Ct", covariates = c("Age", "Ribavirin"))
summary(pegiboot)
plot(pegiboot, treatment = 'treated')

ageboot <- mediate(modCt, predOutcome, sims = 500, boot = TRUE, treat = "Age",
                    mediator = "Ct", covariates = c("Pegivirus", "Ribavirin"))
summary(ageboot)
plot(ageboot, treatment = 'treated')
