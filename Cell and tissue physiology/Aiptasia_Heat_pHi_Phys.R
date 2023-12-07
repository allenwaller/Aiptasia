# Aiptasia physiological analysis
# for: Heat stress disrupts acid-base homeostasis independent of bleaching in the model cnidarian E. diaphana
# Allen-Waller L, KG Jones, MP Martynek, KT Brown, and KL Barott
# compiled by Luella Allen-Waller
# 2023-04-27

# Load packages:
library(ggplot2)
library(csv)
library(tidyverse)
library(ggpubr)
library(readr)
library(Rmisc)
library(lme4)
library(lmerTest)
library(rstatix)
library(gridExtra)
library(data.table)
library(car)
library(patchwork)
library(mgcv)
library(scales)
library(dplyr)
library(MuMIn)
library(devtools)
library(ggbiplot)
library(vegan)


############# Baseline differences between cohorts
# Load data
data <- read.csv("AiptasiaCombinedPhys.csv")
ambphys <- read.csv('AiptasiaCombinedPhys.csv') %>% subset(target.temp.fac == "25")
colors <- read.csv('AiptasiaCombinedPhys.csv') %>% subset(date == "20220202" | date == "20220517")
sizes <- read.csv('AiptasiaSizes.csv')

# Adjust factor names
data$target.temp.fac <- as.factor(data$target.temp.fac)
data$month <- as.factor(data$month)
data$symb.status <- as.factor(data$symb.status)
data$cohort <- dplyr::recode(data$month, May = "HD", Feb = "LD")
colnames(ambphys)[which(names(ambphys) == "month")] <- "Cohort"
ambphys$Cohort <- dplyr::recode(ambphys$Cohort, May = "HD", Feb = "LD")
colors$temp.fac <- as.factor(colors$temp)
colors$date.fac <- as.factor(colors$date)
# Note that zub (Ziploc tub) denotes anemone container

### Initial color
# Compare day 0 color between cohorts
# test whether red color differed between cohorts
red.lm <- lm(red_percent_intensity~cohort, data = colors)
red.lme <- lmer(red_percent_intensity~cohort + (1|zub), data = colors)
AICc(red.lm, red.lme)
# lme is better
# check residuals:
plot(red.lme)
# ANOVA:
anova(red.lme)
# ***cohort (num df = 1, den df = 101.48, F = 18.78, P<0.001)
# Plot:
InitialRedScore <- ggplot(colors, aes(x = cohort, y = 100 - red_percent_intensity, group = cohort)) +
  geom_boxplot(outlier.shape = NA, aes(fill = cohort)) +
  geom_point(size = 2.5, alpha = 0.6, position = position_jitter(width = 0.05),
             aes(shape = cohort)) +
  scale_shape_manual("Cohort", values = c("HD" = 16, "LD" = 1)) +
  scale_fill_manual("Cohort", values = c("HD" = "grey", "LD" = "white")) +
  labs(y = expression(paste("Color score (%)")), x = "") +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 12, color = "black"),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "right", legend.key = element_rect(fill = NA)) +
  coord_cartesian(ylim = c(5,85)) 
InitialRedScore

### Initial symbiont density 
# Compare symbiont counts per anemone in 25ºC controls between cohorts
# We are interested in differences between symbiotic anemones in either cohort, so subset for sym animals:
ambsymphys <- ambphys %>% subset(symb.status == "sym")
symb.lm <- lm(symb.cells.anemone~Cohort, data = ambsymphys)
symb.lme <- lmer(symb.cells.anemone~Cohort + (1|zub), data = ambsymphys)
AICc(symb.lm, symb.lme)
# lme better
plot(symb.lme)
anova(symb.lme)
# num df = 1, den df = 4.5 x 10^37, F = 3.96, p=0.047
# Plot symbiont cells per anemone
Symbs <- ggplot(ambsymphys, aes(x = Cohort, y = symb.cells.anemone / 1000000, group = Cohort)) +
  geom_boxplot(outlier.shape = NA, aes(fill = Cohort)) +
  geom_point(size = 2.5, alpha = 0.6, position = position_jitter(width = 0.05),
             aes(shape = Cohort)) +
  scale_shape_manual("Cohort", values = c("HD" = 16, "LD" = 1)) +
  scale_fill_manual("Cohort", values = c("HD" = "grey", "LD" = "white")) +
  labs(y = expression(atop("Symbiont density", paste((10^6~cells~anemone^-1)))), x = "") +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 12, color = "black"),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "right", legend.key = element_rect(fill = NA)) +
  coord_cartesian(ylim =c(0,3))
Symbs
# cells/ug protein:
symb.t <- t.test(symb.per.prot~Cohort, data = ambsymphys)
symb.t

### Protein
# Compare protein per anemone between symbiotic anemones from both cohorts
sym <- data %>% subset(symb.status == "sym")
prot.lm <- lm(prot.ug.anemone~cohort, data = sym)
prot.lme <- lmer(prot.ug.anemone~cohort + (1|zub), data = sym)
AICc(prot.lm, prot.lme)
# lme is better
plot(prot.lme)
anova(prot.lme)
# num df = 1, den df = 44, F = 75.94, P<0.001
# Plot protein per anemone
Prot <- ggplot(sym, aes(x = cohort, y = prot.ug.anemone/1000, group = cohort)) +
  geom_boxplot(outlier.shape = NA, aes(fill = cohort)) +
  geom_point(size = 2.5, alpha = 0.6, position = position_jitter(width = 0.05),
             aes(shape = cohort)) +
  scale_shape_manual("Cohort", values = c("HD" = 16, "LD" = 1)) +
  scale_fill_manual("Cohort", values = c("HD" = "grey", "LD" = "white")) +
  labs(y = expression(atop("Protein", paste((mg~anemone^-1)))), x = "") +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 12, color = "black"),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "right", legend.key = element_rect(fill = NA)) +
  coord_cartesian(ylim =c(0,1.2)) +
  scale_y_continuous(breaks = c(0.0,0.2,0.4,0.6,0.8,1.0, 1.2)) +
  scale_x_discrete(limits = c("HD", "LD"))
Prot

# test whether apo protein differed between cohorts
apophys <- phys %>% subset(symb.status == "apo")
apoprot.lm <- lm(prot.ug.anemone~Cohort, data = apophys)
apoprot.lme <- lmer(prot.ug.anemone~Cohort + (1|zub), data = apophys)
AICc(apoprot.lm, apoprot.lme)
# lme is better
plot(apoprot.lme)
anova(apoprot.lme)
# ***cohort (F=33.098, P<0.001)
# Visualize:
ApoProt <- ggplot(apophys, aes(x = Cohort, y = prot.ug.anemone /1000, group = Cohort)) +
  geom_boxplot(outlier.shape = NA,aes(fill = Cohort)) +
  geom_point(size = 2.5, alpha = 0.6, position = position_jitter(width = 0.05),
             aes(shape = Cohort)) +
  scale_shape_manual("Cohort", values = c("HD" = 16, "LD" = 1)) +
  scale_fill_manual("Cohort", values = c("HD" = "grey", "LD" = "white")) +
  labs(y = expression(atop("Protein", paste((mg~anemone^-1)))), x = "") +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 12, color = "black"),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "right", legend.key = element_rect(fill = NA)) +
  coord_cartesian(ylim =c(0,0.4)) +
  scale_y_continuous(breaks = c(0.0,0.2,0.4))
ApoProt

###### Size
symsizes <- size %>% subset(symb.status == "sym")
# Test whether oral disk diameters of sym anemones differed between cohorts
size.lm <- lm(diameter.cm~cohort, data = symsizes)
size.lme <- lmer(diameter.cm~cohort + (1|zub), data = symsizes)
AICc(size.lm, size.lme)
# lm is better
plot(size.lm)
# nice
anova(size.lm)
# ***cohort (F=34.568, df = 1, p<0.0001)
# plot oral disk size by cohort
Sizes <- ggplot(symsizes, aes(x = cohort, y = diameter.cm, group = cohort)) +
  geom_boxplot(outlier.shape = NA, aes(fill = cohort)) +
  geom_point(size = 2.5, alpha = 0.6, position = position_jitter(width = 0.05),
             aes(shape = cohort)) +
  scale_shape_manual("Cohort", values = c("HD" = 16, "LD" = 1)) +
  scale_fill_manual("Cohort", values = c("HD" = "grey", "LD" = "white")) +
  labs(y = expression(paste("Oral disk diameter (cm)")), x = "") +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 12, color = "black"),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "right", legend.key = element_rect(fill = NA)) +
  coord_cartesian(ylim = c(0.1, 0.7))
Sizes

### Fig 1: Cohort comparisons
Fig1Phys <- ggarrange(InitialRedScore,Symbs,Sizes,Prot, nrow = 1, ncol = 4, widths = c(0.95,1,1,1.05), common.legend = T)
Fig1Phys


### Verifying that 13C pulse chase worked
# Expectation: enrichment in symbiotic anemones from both cohorts relative to apo controls
# Import raw elemental analysis data with factors from key
C13 <- read.csv("CC7 13C EA-iRMS.csv")
# Adjust factor names
C13$target.temp.fac <- as.factor(C13$target.temp.fac)
C13$fraction <- dplyr::recode(C13$fraction, host = "Host", symbiont = "Symbiont")
C13$symb.status <- dplyr::recode(C13$symb.status, apo = "Apo", sym = "Sym")

# Plot by fraction x symb status (Fig. S3)
HostSymb13C <- ggplot(C13, aes(x = cohort, y = av_at_percent_13C)) +
  geom_boxplot(position=position_dodge(width=0.75),
               aes(color = fraction, linetype = symb.status,
                   group = interaction(cohort, fraction, symb.status)),
               outlier.shape = NA) +
  geom_point(size = 2.5, alpha = 0.6, position = position_dodge(width = 0.75),
             aes(color = fraction, shape = interaction(symb.status,cohort),
                 group = interaction(cohort, fraction, symb.status))) +
  ylab(expression(paste("13C atom-%"))) + 
  xlab("") +
  labs(color = "Fraction", shape = "Symbiotic Status", linetype = "Symbiotic Status") +
  scale_color_manual(values = c("salmon", "#7b5313")) +
  scale_shape_manual(values = c(17,19,2,1)) +
  scale_linetype_manual(values = c("dotted", "solid")) +
  theme_bw()  +
  theme(panel.grid = element_blank(), text = element_text(size = 12), legend.position = "right") +
  coord_cartesian(ylim = c(1.095,1.17))
HostSymb13C

# Testing whether fraction, cohort, symbiotic status affected 13C assimilation
symbhost13C.lmm <- lmer(av_at_percent_13C ~ fraction * cohort * symb.status + (1|zub), data = C13) # rank deficient
symbhost13C.lmm2 <- lmer(av_at_percent_13C ~ fraction * symb.status + (1|zub), data = C13) # rank deficient
symbhost13C.lm2 <- lm(av_at_percent_13C ~ fraction * symb.status, data = C13)
symbhost13C.lm <- lm(av_at_percent_13C ~ fraction * cohort * symb.status, data = C13)
# Test models
AICc(symbhost13C.lmm,symbhost13C.lmm2,symbhost13C.lm2,symbhost13C.lm)
# symbhost13C.lm is best, and simplest
plot(symbhost13C.lm)
anova(symbhost13C.lm)
# ***fraction, *fraction x cohort, *cohort, *symb status
# Posthoc Tukey pairwise comparisons
# By cohort:
symbhost13C.tuk <- emmeans(symbhost13C.lm, list(pairwise ~ fraction*cohort*symb.status), simple = "cohort", adjust = "tukey")
symbhost13C.tuk
# Host fractions had same 13C assimilation between cohorts, but symbiont fractions were different
symbhost13C.tuk2 <- emmeans(symbhost13C.lm, list(pairwise ~ fraction*cohort*symb.status), simple = "fraction", adjust = "tukey")
symbhost13C.tuk2
# Host vs. symbiont fractions from sym animals differed significantly
# Since the effect of symb status didn't differ by cohort, test effect of symb status independently of cohort
# but keep fraction since only host fraction was represented by both apo and sym individuals
symbhost13C.tuk3 <- emmeans(symbhost13C.lm2, list(pairwise ~ fraction*symb.status), simple = "symb.status", adjust = "tukey")
symbhost13C.tuk3
# Symbiont status significantly affected host fraction 13C assimilation.



############# Anemone physiology response to temperature incubation
# First, split dataframe by cohort
low <- data %>% subset(cohort == "LD")
high <- data %>% subset(cohort == "HD")

### Protein content

# High-symbiont-density cohort (HD) protein per anemone
# Models:
hdprot.lme <- lmer(prot.ug.anemone~target.temp.fac + (1|zub), data = high)
hdprot.lm <- lm(prot.ug.anemone~target.temp.fac, data = high)
hdprot.lme2 <- lmer(prot.ug.anemone~target.temp.fac*symb.status + (1|zub), data = high)
hdprot.lm2 <- lm(prot.ug.anemone~target.temp.fac*symb.status, data = high)
hdprot.lme3 <- lmer(prot.ug.anemone~symb.status + (1|zub), data = high)
hdprot.lm3 <- lm(prot.ug.anemone~symb.status, data = high)
# Compete
AICc(hdprot.lme,hdprot.lm,hdprot.lme2,hdprot.lm2,hdprot.lme3,hdprot.lm3)
# hdprot.lme2 is best
plot(hdprot.lme2)
anova(hdprot.lme2)
# significant effect of symbiont status (***) but not of temp (p=0.2)
hdprot.tukey <- emmeans(hdprot.lme2, list(pairwise ~ target.temp.fac*symb.status), simple = "symb.status", adjust = "tukey")
hdprot.tukey
# all different based on symb status
# Plot:
HDProtein <- ggplot(high, aes(x = target.temp.fac, y = prot.ug.anemone, group = interaction(target.temp.fac, symb.status))) +
  geom_boxplot(aes(color = target.temp.fac, linetype = symb.status),
               outlier.shape = NA) +
  geom_point(size = 2.5, alpha = 0.6, position = position_dodge(width = 0.75),
             aes(color = target.temp.fac, shape = symb.status)) +
  scale_color_manual("Target Temperature", values = c("25" = "#0571B0", "27" = "#92C5DE", "29" = "#F4A582", "31" = "#CA0020")) +
  scale_shape_manual("Symbiotic Status", values = c(17, 16)) +
  scale_linetype_manual("Symbiotic Status", values = c("dotted", "solid")) +
  labs(y = expression(atop("Protein", paste((µg~anemone^-1)))), x = "Treatment temperature (ºC)") +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 12, color = "black"),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "right", legend.key = element_rect(fill = NA),
        strip.background = element_blank(), strip.text = element_blank()) +
  coord_cartesian(ylim = c(0,1100))
HDProtein

# Low-symbiont-density cohort (LD) protein per anemone
# Models:
# container n = anemone n, grouping factor levels are not < #observations
ldprot.lm <- lm(prot.ug.anemone~target.temp.fac, data = low)
ldprot.lm2 <- lm(prot.ug.anemone~target.temp.fac*symb.status, data = low)
# Compete
AICc(ldprot.lm,ldprot.lm2)
# ldprot.lm2 is the lower
plot(ldprot.lm2)
anova(ldprot.lm2)
# no effect of temp (p=0.8), but yes symb status (***)
tukey.ldprot<- emmeans(ldprot.lm2, list(pairwise ~ symb.status*target.temp.fac), simple = "symb.status", adjust = "tukey")
tukey.ldprot
# apo and sym are different at every temp
tukey.ldprot2<- emmeans(ldprot.lm2, list(pairwise ~ symb.status*target.temp.fac), simple = "target.temp.fac", adjust = "tukey")
tukey.ldprot2
# no temperature differences
# Plot:
LDProtein <- ggplot(low, aes(x = target.temp.fac, y = prot.ug.anemone, group = interaction(target.temp.fac, symb.status))) +
  geom_boxplot(aes(color = target.temp.fac, linetype = symb.status),
               outlier.shape = NA) +
  geom_point(size = 2.5, alpha = 0.6, position = position_dodge(width = 0.75),
             aes(color = target.temp.fac, shape = symb.status)) +
  scale_color_manual("Target Temperature", values = c("25" = "#0571B0", "27" = "#92C5DE", "29" = "#F4A582", "31" = "#CA0020")) +
  scale_shape_manual("Symbiotic Status", values = c(2, 1)) +
  scale_linetype_manual("Symbiotic Status", values = c("dotted", "solid")) +
  labs(y = expression(atop("Protein", paste((µg~anemone^{-1})))), x = "Treatment temperature (ºC)") +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 12, color = "black"),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "right", legend.key = element_rect(fill = NA),
        strip.background = element_blank(), strip.text = element_blank()) +
  coord_cartesian(ylim =c(0,1100))
LDProtein

### Symbiont density (cells/µg protein)

# HD cohort symbs/protein
# Models:
hdsymb.lme <- lmer(symb.per.prot~target.temp.fac + (1|zub), data = high)
hdsymb.lm <- lm(symb.per.prot~target.temp.fac, data = high)
hdsymb2.lme <- lmer(symb.per.prot~target.temp.fac*symb.status + (1|zub), data = high)
hdsymb2.lm <- lm(symb.per.prot~target.temp.fac*symb.status, data = high)
# Compete
AICc(hdsymb.lme,hdsymb.lm,hdsymb2.lme,hdsymb2.lm)
# hdsymb2.lme is the lowest
plot(hdsymb2.lme)
anova(hdsymb2.lme)
# ***temp X symb stat, ***symb stat, ***temp
# pairwise:
symb.tukey <- emmeans(hdsymb2.lme, list(pairwise ~ target.temp.fac*symb.status), simple = "symb.status", adjust = "tukey")
symb.tukey
# apo and sym significantly different at every temperature
symb.tukey.temp <- emmeans(hdsymb2.lme, list(pairwise ~ target.temp.fac*symb.status), simple = "target.temp.fac", adjust = "tukey")
symb.tukey.temp
# for apos, no one is different, so apos = a
# for syms, all different except 25 vs. 27, 27 vs. 29 (p=0.066), and 29 vs. 31
# sym 25 = b, sym 27 = bc, sym 29 = cd, sym 31 = d
# Plot:
HDSymbdens <- ggplot(high, aes(x = target.temp.fac, y = symb.per.prot/1000, group = interaction(target.temp.fac, symb.status))) +
  geom_boxplot(aes(color = target.temp.fac, linetype = symb.status),
               outlier.shape = NA) +
  geom_point(size = 2.5, alpha = 0.6, position = position_dodge(width = 0.75),
             aes(color = target.temp.fac, shape = symb.status)) +
  scale_color_manual("Target Temperature", values = c("25" = "#0571B0", "27" = "#92C5DE", "29" = "#F4A582", "31" = "#CA0020")) +
  scale_shape_manual("Symbiotic Status", values = c(17, 16)) +
  scale_linetype_manual("Symbiotic Status", values = c("dotted", "solid")) +
  labs(y = expression(atop("Symbiont density", paste((10^3~cells~µg~protein^-1)))), x = "Treatment temperature (ºC)") +
  theme(axis.title = element_text(size = 13), axis.text = element_text(size = 13, color = "black"),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "right", legend.key = element_rect(fill = NA),
        strip.background = element_blank(), strip.text = element_blank())+
  coord_cartesian(ylim = c(0,10))
HDSymbdens

# LD cohort symbs/protein
# Models:
# zub n = anemone n so can't have zub RE to avoid grouping factor levels are not < #observations
ldsymb.lm <- lm(symb.per.prot~target.temp.fac, data = low)
ldsymb.lm2 <- lm(symb.per.prot~target.temp.fac*symb.status, data = low)
# Compete
AICc(ldsymb.lm,ldsymb.lm2)
# ldsymb.lm2 is the lower
plot(ldsymb.lm2)
anova(ldsymb.lm2)
# highly significant effect of symbiotic status (***)
tukey.ldsymb <- emmeans(ldsymb.lm2, list(pairwise ~ symb.status*target.temp.fac), simple = "symb.status", adjust = "tukey")
tukey.ldsymb
# apo and sym are different at 25 and 31
# Plot:
LDSymbdens <- ggplot(low, aes(x = target.temp.fac, y = symb.per.prot/1000, group = interaction(target.temp.fac, symb.status))) +
  geom_boxplot(aes(color = target.temp.fac, linetype = symb.status),
               outlier.shape = NA) +
  geom_point(size = 2.5, alpha = 0.6, position = position_dodge(width = 0.75),
             aes(color = target.temp.fac, shape = symb.status)) +
  scale_color_manual("Target Temperature", values = c("25" = "#0571B0", "27" = "#92C5DE", "29" = "#F4A582", "31" = "#CA0020")) +
  scale_shape_manual("Symbiotic Status", values = c(2, 1)) +
  scale_linetype_manual("Symbiotic Status", values = c("dotted", "solid")) +
  labs(y = expression(atop("Symbiont density", paste((10^3~cells~µg~protein^-1)))), x = "Treatment temperature (ºC)") +
  theme(axis.title = element_text(size = 13), axis.text = element_text(size = 13, color = "black"),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "right", legend.key = element_rect(fill = NA),
        strip.background = element_blank(), strip.text = element_blank())+
  coord_cartesian(ylim= c(0,10))
LDSymbdens


### Symbiont density (cells/anemone)

# HD cohort symbs per anemone
# Models:
hd.totalsymb.lm <- lm(symb.cells.anemone ~ target.temp.fac, data = high)
hd.totalsymb.lm2 <- lm(symb.cells.anemone ~ target.temp.fac*symb.status, data = high)
hd.totalsymb.lm3 <- lm(symb.cells.anemone ~ symb.status, data = high)
hd.totalsymb.lmer <- lmer(symb.cells.anemone ~ target.temp.fac + (1|zub), data = high)
hd.totalsymb.lmer2 <- lmer(symb.cells.anemone ~ target.temp.fac*symb.status + (1|zub), data = high)
hd.totalsymb.lmer3 <- lmer(symb.cells.anemone ~ symb.status + (1|zub), data = high)
AICc(hd.totalsymb.lm,hd.totalsymb.lm2,hd.totalsymb.lm3,hd.totalsymb.lmer,hd.totalsymb.lmer2,hd.totalsymb.lmer3)
# hd.totalsymb.lmer2 is best
plot(hd.totalsymb.lmer2)
anova(hd.totalsymb.lmer2)
# strong interaction as well as both temp and symb status
hd.totalsymb.lmer2.tukey <- emmeans(hd.totalsymb.lmer2, list(pairwise ~ target.temp.fac*symb.status), simple = "target.temp.fac", adjust = "tukey")
hd.totalsymb.lmer2.tukey
# apos and syms differ at every temperature
# apos all the same as each other
# syms differ: all except 25 vs. 27 and 29 vs. 31
# Plot:
HDSymbsPerAnimal <- ggplot(high, aes(x = target.temp.fac, y = symb.cells.anemone/1000000, group = interaction(target.temp.fac, symb.status))) +
  geom_boxplot(aes(color = target.temp.fac, linetype = symb.status),
               outlier.shape = NA) +
  geom_point(size = 2.5, alpha = 0.6, position = position_dodge(width = 0.75),
             aes(color = target.temp.fac, shape = symb.status)) +
  scale_color_manual("Target Temperature", values = c("25" = "#0571B0", "27" = "#92C5DE", "29" = "#F4A582", "31" = "#CA0020")) +
  scale_shape_manual("Symbiotic Status", values = c(17, 16)) +
  scale_linetype_manual("Symbiotic Status", values = c("dotted", "solid")) +
  labs(y = expression(atop("Symbiont total population", paste((10^6~cells~anemone^-1)))), x = "Treatment temperature (ºC)") +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 12, color = "black"),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "right", legend.key = element_rect(fill = NA),
        strip.background = element_blank(), strip.text = element_blank())+
  coord_cartesian(ylim =c(0,2.8))
HDSymbsPerAnimal

# LD cohort symbs per anemone
# Models:
ld.totalsymb.lm <- lm(symb.cells.anemone ~ target.temp.fac, data = low)
ld.totalsymb.lm2 <- lm(symb.cells.anemone ~ target.temp.fac*symb.status, data = low)
ld.totalsymb.lm3 <- lm(symb.cells.anemone ~ symb.status, data = low)
# not enough zubs in the LD experiment so just comparing LMs
AICc(ld.totalsymb.lm,ld.totalsymb.lm2,ld.totalsymb.lm3)
# the one only looking at symb status is best
plot(ld.totalsymb.lm3)
anova(ld.totalsymb.lm3) # *** Symb status
# Verify no differences within symb status:
ld.totalsymb.lm2.tukey <- emmeans(ld.totalsymb.lm2, list(pairwise ~ target.temp.fac*symb.status), adjust = "tukey")
ld.totalsymb.lm2.tukey
# apos and syms all different, but within symb statuses, no differences
ld.totalsymb.lm3.tukey <- emmeans(ld.totalsymb.lm3, list(pairwise ~ symb.status), adjust = "tukey")
ld.totalsymb.lm3.tukey
# yes, p<0.0001 for apo vs. sym
# Plot:
LDSymbsPerAnimal <- ggplot(low, aes(x = target.temp.fac, y = symb.cells.anemone /1000000, group = interaction(target.temp.fac, symb.status))) +
  geom_boxplot(aes(color = target.temp.fac, linetype = symb.status),
               outlier.shape = NA) +
  geom_point(size = 2.5, alpha = 0.6, position = position_dodge(width = 0.75),
             aes(color = target.temp.fac, shape = symb.status)) +
  scale_color_manual("Target Temperature", values = c("25" = "#0571B0", "27" = "#92C5DE", "29" = "#F4A582", "31" = "#CA0020")) +
  scale_shape_manual("Symbiotic Status", values = c(2, 1)) +
  scale_linetype_manual("Symbiotic Status", values = c("dotted", "solid")) +
  labs(y = expression(atop("Symbiont total population", paste((10^6~cells~anemone^-1)))), x = "Treatment temperature (ºC)") +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 12, color = "black"),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "right", legend.key = element_rect(fill = NA),
        strip.background = element_blank(), strip.text = element_blank())+
  coord_cartesian(ylim =c(0,2.8))
LDSymbsPerAnimal

### Symbiont 13C

# HD cohort symb 13C:
# no symb status, and zub = total count, so just one linear model
hdsymb13C.lm <- lm(symb.13C.percent~target.temp.fac, data = high)
plot(hdsymb13C.lm)
anova(hdsymb13C.lm)
# **temp
# pairwise:
symb13C.tukey <- emmeans(hdsymb13C.lm, list(pairwise ~ target.temp.fac), adjust = "tukey")
symb13C.tukey
# same as symbs - sym 25 = a, sym 27 = ab, sym 29 = bc, sym 31 = c
# Plot:
HDSymb13C <- ggplot(high, aes(x = target.temp.fac, y = symb.13C.percent, group = target.temp.fac)) +
  geom_boxplot(aes(color = target.temp.fac),
               outlier.shape = NA) +
  geom_point(size = 2.5, alpha = 0.6, position = position_dodge(width = 0.75),
             aes(color = target.temp.fac)) +
  scale_color_manual("Target Temperature", values = c("25" = "#0571B0", "27" = "#92C5DE", "29" = "#F4A582", "31" = "#CA0020")) +
  labs(y = expression(atop("Symbiont C assimilation", paste('(atom-% '^{13},'C)'))), x = "Treatment temperature (ºC)") +
  theme(axis.title = element_text(size = 13), axis.text = element_text(size = 13, color = "black"),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "right", legend.key = element_rect(fill = NA),
        strip.background = element_blank(), strip.text = element_blank())+
  coord_cartesian(ylim = c(1.1,1.165))
HDSymb13C

# LD cohort symb 13C:
# Model:
# zub n = anemone n, grouping factor levels are not < #observations
ldsymb13C.lm <- lm(symb.13C.percent~target.temp.fac, data = low)
plot(ldsymb13C.lm)
anova(ldsymb13C.lm)
# significant effect of temperature (**)
tukey.ldsymb13C <- emmeans(ldsymb13C.lm, list(pairwise ~ target.temp.fac), adjust = "tukey")
tukey.ldsymb13C
# 25 and 31 were different, 27 and 29 are different, and 29 is different from 31
# Plot:
LDSymb13C <- ggplot(low, aes(x = target.temp.fac, y = symb.13C.percent, group = target.temp.fac)) +
  geom_boxplot(aes(color = target.temp.fac),
               outlier.shape = NA) +
  geom_point(size = 2.5, shape = 1, alpha = 0.6, position = position_dodge(width = 0.75),
             aes(color = target.temp.fac)) +
  scale_color_manual("Target Temperature", values = c("25" = "#0571B0", "27" = "#92C5DE", "29" = "#F4A582", "31" = "#CA0020")) +
  labs(y = expression(atop("Symbiont C assimilation", paste('(atom-% '^{13},'C)'))), x = "Treatment temperature (ºC)") +
  theme(axis.title = element_text(size = 13), axis.text = element_text(size = 13, color = "black"),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "right", legend.key = element_rect(fill = NA),
        strip.background = element_blank(), strip.text = element_blank())+
  coord_cartesian(ylim = c(1.1,1.165))
LDSymb13C

### Host 13C

# HD cohort host 13C
# Model:
# no symb status, and zub = total count, so just one linear model
hdhost13C.lm <- lm(host.13C.percent~target.temp.fac, data = high)
plot(hdhost13C.lm)
anova(hdhost13C.lm)
# p=0.060 for temperature
# pairwise:
host13C.tukey <- emmeans(hdhost13C.lm, list(pairwise ~ target.temp.fac), adjust = "tukey")
host13C.tukey
# significant pairwise difference between 25 and 31
# Plot:
HDHost13C <- ggplot(high, aes(x = target.temp.fac, y = host.13C.percent, group = target.temp.fac)) +
  geom_boxplot(aes(color = target.temp.fac),
               outlier.shape = NA) +
  geom_point(size = 2.5, alpha = 0.6, position = position_dodge(width = 0.75),
             aes(color = target.temp.fac)) +
  scale_color_manual("Target Temperature", values = c("25" = "#0571B0", "27" = "#92C5DE", "29" = "#F4A582", "31" = "#CA0020")) +
  labs(y = expression(atop("Host C assimilation", paste('(atom-% '^{13},'C)'))), x = "Treatment temperature (ºC)") +
  theme(axis.title = element_text(size = 13), axis.text = element_text(size = 13, color = "black"),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "right", legend.key = element_rect(fill = NA),
        strip.background = element_blank(), strip.text = element_blank()) +
  coord_cartesian(ylim = c(1.1,1.135))
HDHost13C

# LD cohort host 13C
# Model:
# zub n = anemone n, grouping factor levels are not < #observations
ldhost13C.lm <- lm(host.13C.percent~target.temp.fac, data = low)
plot(ldhost13C.lm)
anova(ldhost13C.lm)
# no effect
# Plot:
LDHost13C <- ggplot(low, aes(x = target.temp.fac, y = host.13C.percent, group = target.temp.fac)) +
  geom_boxplot(aes(color = target.temp.fac),
               outlier.shape = NA) +
  geom_point(size = 2.5, shape = 1, alpha = 0.6, position = position_dodge(width = 0.75),
             aes(color = target.temp.fac)) +
  scale_color_manual("Target Temperature", values = c("25" = "#0571B0", "27" = "#92C5DE", "29" = "#F4A582", "31" = "#CA0020")) +
  labs(y = expression(atop("Host C assimilation", paste('(atom-% '^{13},'C)'))), x = "Treatment temperature (ºC)") +
  theme(axis.title = element_text(size = 13), axis.text = element_text(size = 13, color = "black"),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "right", legend.key = element_rect(fill = NA),
        strip.background = element_blank(), strip.text = element_blank()) +
  coord_cartesian(ylim = c(1.1,1.135))
LDHost13C

### Fig. 3: Symbionts per protein, symbiont 13C, and host 13C
SymbDensityFunction <- ggarrange(HDSymbdens + rremove("xlab") + rremove("x.text"), 
                                 HDSymb13C + rremove("xlab")+ rremove("x.text"), 
                                 HDHost13C + rremove("xlab")+ rremove("x.text"),
                                 LDSymbdens + rremove("xlab"), 
                                 LDSymb13C + rremove("xlab"), 
                                 LDHost13C + rremove("xlab"),
                                 common.legend = T, legend = "top",
                                 nrow = 2, ncol = 3, heights = 2.4, widths = 4)
SymbDensityFunction

### Fig. S2: Protein per anemone and symbiont cells per anemone
ProtSymbsPerAnemone <- ggarrange(HDProtein, HDSymbsPerAnimal,LDProtein, LDSymbsPerAnimal, nrow = 2, ncol = 2, common.legend = T)
ProtSymbsPerAnemone


############# Intracellular pH

### Symbiocyte pHi

# HD symbiocyte pHi
# Model:
# no symb status and zub = total count so just one linear model
hdsymbpHi.lm <- lm(symb.pHi.med~target.temp.fac, data = high)
plot(hdsymbpHi.lm)
anova(hdsymbpHi.lm)
# p=0.07265 for temp
symbiocytepHi.tukey <- emmeans(hdsymbpHi.lm, list(pairwise ~ target.temp.fac), adjust = "tukey")
symbiocytepHi.tukey
# p=0.087 for 29º vs. 31º, that's all
# Plot:
HDSymbiocytepHi <- ggplot(high, aes(x = target.temp.fac, y = symb.pHi.med, group = target.temp.fac)) +
  geom_boxplot(aes(color = target.temp.fac),
               outlier.shape = NA) +
  geom_point(size = 2.5, alpha = 0.6, position = position_dodge(width = 0.75),
             aes(color = target.temp.fac)) +
  scale_color_manual("Target Temperature", values = c("25" = "#0571B0", "27" = "#92C5DE", "29" = "#F4A582", "31" = "#CA0020")) +
  labs(y = expression(atop("Symbiocyte", paste('Intracellular pH'))), x = "Treatment temperature (ºC)") +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 12, color = "black"),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "right", legend.key = element_rect(fill = NA),
        strip.background = element_blank(), strip.text = element_blank()) +
  coord_cartesian(ylim =c(6.7,7.8))
HDSymbiocytepHi

# LD symbiocyte pHi
# Model:
# zub n = anemone n, grouping factor levels are not < #observations
ldsymb.lm <- lm(symb.pHi.med~target.temp.fac, data = low)
plot(ldsymb.lm)
anova(ldsymb.lm)
# significant effect of temperature (*)
tukey.ldsymbpHi <- emmeans(ldsymb.lm, list(pairwise ~ target.temp.fac), adjust = "tukey")
tukey.ldsymbpHi
# 25 and 31 were different and that is it (27 and 31 borderline)
# Plot:
LDSymbiocytepHi <- ggplot(low, aes(x = target.temp.fac, y = symb.pHi.med, group = target.temp.fac)) +
  geom_boxplot(aes(color = target.temp.fac),
               outlier.shape = NA) +
  geom_point(size = 2.5, alpha = 0.6, shape = 1, position = position_dodge(width = 0.75),
             aes(color = target.temp.fac)) +
  scale_color_manual("Target Temperature", values = c("25" = "#0571B0", "27" = "#92C5DE", "29" = "#F4A582", "31" = "#CA0020")) +
  labs(y = expression(atop("Symbiocyte", paste('Intracellular pH'))), x = "Treatment temperature (ºC)") +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 12, color = "black"),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "right", legend.key = element_rect(fill = NA),
        strip.background = element_blank(), strip.text = element_blank()) +
  coord_cartesian(ylim =c(6.7,7.8))
LDSymbiocytepHi

### Regression: Host 13C vs. symbiocyte pHi
# Model:
symbpHi.all.lm <- lm(symb.pHi.med ~ host.13C.percent, data = data)
symbpHi.all.lm2 <- lm(symb.pHi.med ~ host.13C.percent*cohort, data = data)
symbpHi.all.lm3 <- lm(symb.pHi.med ~ host.13C.percent*target.temp.fac, data = data)
symbpHi.all.lm4 <- lm(symb.pHi.med ~ host.13C.percent*cohort*target.temp.fac, data = data)
AICc(symbpHi.all.lm,symbpHi.all.lm2,symbpHi.all.lm3,symbpHi.all.lm4)
# best is 13C x cohort
anova(symbpHi.all.lm2)
# * sig interaction between cohort and 13C
# so, 13C differentially affected symbiocyte intracellular pH depending on cohort
# HD symbiocyte pHi vs. assimilation lm:
HDsymbpHi.vs.13C <- lm(symb.pHi.med ~ host.13C.percent, data = high)
summary(HDsymbpHi.vs.13C)
# LD symbiocyte pHi vs. assimilation lm:
LDsymbpHi.vs.13C <- lm(symb.pHi.med ~ host.13C.percent, data = low)
summary(LDsymbpHi.vs.13C)
# Plot:
Host13C.vs.SymbpHi <- ggplot(data, aes(y=symb.pHi.med, x = host.13C.percent, color = target.temp.fac, shape = cohort, linetype = month)) +
  geom_point(size = 2.5, alpha = 0.8, aes(color = target.temp.fac)) +
  geom_smooth(method = "lm", formula = y ~ x,
              color = "black", size = 0.5, alpha = 0.3) +
  scale_color_manual(values = c("#0571B0", "#92C5DE", "#F4A582", "#CA0020")) +
  scale_shape_manual("Cohort", values = c(1, 19)) +
  scale_linetype_manual("Cohort", values = c("longdash","solid")) +
  xlab(expression(atop("Host photosynthate assimilation", paste('(atom-% '^{13},'C)')))) + ylab("Symbiocyte Intracellular pH") +
  theme_bw()  + 
  theme(panel.grid = element_blank(), text = element_text(size = 14), axis.text = element_text(color = "black")) +
  coord_cartesian(ylim = c(6.7,7.8))
Host13C.vs.SymbpHi

### Nonsymbiocyte pHi

# HD nonsymbiocyte pHi
# Model:
# zub count = total count so just LMs
hdnonsymbpHi.lm <- lm(nonsymb.pHi.med~target.temp.fac, data = high)
hdnonsymbpHi.lm2 <- lm(nonsymb.pHi.med~target.temp.fac*symb.status, data = high)
AICc(hdnonsymbpHi.lm,hdnonsymbpHi.lm2)
# hdnonsymbpHi.lm is better
plot(hdnonsymbpHi.lm)
anova(hdnonsymbpHi.lm)
# *** temperature 
nonsymbiocytepHi.tukey <- emmeans(hdnonsymbpHi.lm2, list(pairwise ~ target.temp.fac), adjust = "tukey")
nonsymbiocytepHi.tukey
# 25 significantly different from 27; 29 significantly different from 31 and 27
# put differently: 25 is the same as 29 and 31, 27 is the same as 31, 27 is different from 29 and 25
# 25 = ab, 27 = c, 29 = a, 31 = bc
# check that symb statuses don't differ pairwise at any temp:
nonsymbiocytepHi.tukey2 <- emmeans(hdnonsymbpHi.lm2, list(pairwise ~ symb.status*target.temp.fac), simple="symb.status", adjust = "tukey")
nonsymbiocytepHi.tukey2
# confirmed, they do not
# Plot:
HDNonsymbiocytepHi <- ggplot(high, aes(x = target.temp.fac, y = nonsymb.pHi.med, group = interaction(target.temp.fac, symb.status))) +
  geom_boxplot(aes(color = target.temp.fac, linetype = symb.status),
               outlier.shape = NA) +
  geom_point(size = 2.5, alpha = 0.6, position = position_dodge(width = 0.75),
             aes(color = target.temp.fac, shape = symb.status)) +
  scale_color_manual("Target Temperature", values = c("25" = "#0571B0", "27" = "#92C5DE", "29" = "#F4A582", "31" = "#CA0020")) +
  scale_shape_manual("Symbiotic Status", values = c(17, 16)) +
  scale_linetype_manual("Symbiotic Status", values = c("dotted", "solid")) +
  labs(y = expression(atop("Nonsymbiocyte", paste('Intracellular pH'))), x = "Treatment temperature (ºC)") +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 12, color = "black"),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "right", legend.key = element_rect(fill = NA),
        strip.background = element_blank(), strip.text = element_blank()) +
  coord_cartesian(ylim =c(6.7,7.8))
HDNonsymbiocytepHi

# LD nonsymbiocyte pHi
# Model:
# zub n = anemone n, grouping factor levels are not < #observations
ldnonsymb.lm <- lm(nonsymb.pHi.med~target.temp.fac, data = low)
ldnonsymb.lm2 <- lm(nonsymb.pHi.med~target.temp.fac*symb.status, data = low)
# Compete
AICc(ldnonsymb.lm,ldnonsymb.lm2)
# ldsymb.lm is the lower
plot(ldnonsymb.lm)
anova(ldnonsymb.lm)
# highly significant effect of temperature (***)
# other model:
plot(ldnonsymb.lm2)
anova(ldnonsymb.lm2)
# highly significant effect of temperature (***)
# tukey:
tukey.ldnonsymbpHi <- emmeans(ldnonsymb.lm2, list(pairwise ~ target.temp.fac), simple = "target.temp.fac", adjust = "tukey")
tukey.ldnonsymbpHi
# 25 differed from 27 and 31. 29 also differed from 31
# so 25 = a, 27 = bc, 29 = ab, 31 = c
# check that sym and nonsym aren't different at any temp:
tukey.ldnonsymbpHi2 <- emmeans(ldnonsymb.lm2, list(pairwise ~ target.temp.fac*symb.status), simple = "symb.status", adjust = "tukey")
tukey.ldnonsymbpHi2
# correct, no pairwise differences
# Plot:
LDNonsymbiocytepHi <- ggplot(low, aes(x = target.temp.fac, y = nonsymb.pHi.med, group = interaction(target.temp.fac, symb.status))) +
  geom_boxplot(aes(color = target.temp.fac, linetype = symb.status),
               outlier.shape = NA) +
  geom_point(size = 2.5, alpha = 0.6, position = position_dodge(width = 0.75),
             aes(color = target.temp.fac, shape = symb.status)) +
  scale_color_manual("Target Temperature", values = c("25" = "#0571B0", "27" = "#92C5DE", "29" = "#F4A582", "31" = "#CA0020")) +
  scale_shape_manual("Symbiotic Status", values = c(2, 1)) +
  scale_linetype_manual("Symbiotic Status", values = c("dotted", "solid")) +
  labs(y = expression(atop("Nonsymbiocyte", paste('Intracellular pH'))), x = "Treatment temperature (ºC)") +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 12, color = "black"),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "right", legend.key = element_rect(fill = NA),
        strip.background = element_blank(), strip.text = element_blank()) +
  coord_cartesian(ylim =c(6.7,7.8))
LDNonsymbiocytepHi

# Regression: Host 13C vs. nonsymbiocyte pHi
# Model:
nonsymbpHi.all.lm <- lm(nonsymb.pHi.med ~ host.13C.percent, data = data)
nonsymbpHi.all.lm2 <- lm(nonsymb.pHi.med ~ host.13C.percent*cohort, data = data)
nonsymbpHi.all.lm3 <- lm(nonsymb.pHi.med ~ host.13C.percent*target.temp.fac, data = data)
nonsymbpHi.all.lm4 <- lm(nonsymb.pHi.med ~ host.13C.percent*cohort*target.temp.fac, data = data)
AICc(nonsymbpHi.all.lm,nonsymbpHi.all.lm2,nonsymbpHi.all.lm3,nonsymbpHi.all.lm4)
# best is lm2
anova(nonsymbpHi.all.lm2) # ***sig effect of cohort
summary(nonsymbpHi.all.lm2)
# F3,58 =5.781, p=0.0016
# HD nonsymbiocyte pHi vs. assimilation lm:
HDnonsymbpHi.vs.13C <- lm(nonsymb.pHi.med ~ host.13C.percent, data = high)
summary(HDnonsymbpHi.vs.13C)
# LD nonsymbiocyte pHi vs. assimilation lm:
LDnonsymbpHi.vs.13C <- lm(nonsymb.pHi.med ~ host.13C.percent, data = low)
summary(LDnonsymbpHi.vs.13C)
# Plot:
Host13C.vs.NonsymbpHi <- ggplot(data, aes(y=nonsymb.pHi.med, x = host.13C.percent, color = target.temp.fac, shape = cohort, linetype = month)) +
  geom_point(size = 2.5, alpha = 0.8, aes(color = target.temp.fac)) +
  geom_smooth(method = "lm", formula = y ~ x,
              color = "black", size = 0.5, alpha = 0.3) +
  scale_color_manual(values = c("#0571B0", "#92C5DE", "#F4A582", "#CA0020")) +
  scale_shape_manual("Cohort", values = c(1, 19)) +
  scale_linetype_manual("Cohort", values = c("longdash","solid")) +
  xlab(expression(atop("Host photosynthate assimilation", paste('(atom-% '^{13},'C)')))) + ylab("Nonsymbiocyte Intracellular pH") +
  theme_bw()  + 
  theme(panel.grid = element_blank(), text = element_text(size = 14), axis.text = element_text(color = "black"),
        legend.position = "top") +
  coord_cartesian(ylim = c(6.7,7.8))
Host13C.vs.NonsymbpHi


### Fig. 5: 
Total.pHi <- ggarrange(HDSymbiocytepHi + rremove("x.text") + rremove("xlab"), 
                       LDSymbiocytepHi + rremove("y.text") + rremove("ylab") + rremove("x.text") + rremove("xlab"),
                       Host13C.vs.SymbpHi + rremove("y.text") + rremove("ylab") +rremove("x.text") + rremove("xlab"),
                       HDNonsymbiocytepHi + rremove("xlab"), 
                       LDNonsymbiocytepHi + rremove("y.text") + rremove("ylab") + rremove("xlab"),
                       Host13C.vs.NonsymbpHi+ rremove("y.text") + rremove("ylab") + rremove("xlab"), 
                       nrow = 2, ncol = 3, widths = c(1.2,1,1), heights = c(1,1.06),
                       common.legend = T)
Total.pHi



############# Principal components analysis
# Load data
pca <-  read.csv("MDS CC7 Combined Phys.csv")
# Separate by cohort, trim & make sure factors are correct
pca$target.temp.fac<-as.factor(pca$target.temp.fac)
pca$zub<-as.factor(pca$zub)
pca$symb.status<-as.factor(pca$symb.status)
pca$month<-as.factor(pca$month)
lowdens <- pca %>% subset(month == "Feb")
highdens <- pca %>% subset(month == "May")
# Low-density
# select columns needed (exclude respirometry for Low-density, since we did not measure resp for those) and drop rows with NAs
lowdens.multi <- dplyr::select(lowdens, c("month","zub","target.temp.fac","symb.status",
                                          "zav.symb.per.prot","zav.prot.ug.anemone","nonsymb.pHi.med","symb.pHi.med",
                                          "host.13C.percent","symb.13C.percent")) %>% drop_na()
colnames(lowdens.multi)[1:10] <- c("Cohort", "Zub", "Temp", "SymbStatus","Sym", "Protein", "Nonsymb pHi", "Symbiocyte pHi", 
                                   "Host13C", "Symb13C")
# Scale and center data
low.scaled <- scale(lowdens.multi[5:10], center = T, scale = T) 
# Identify factors 
fac.low <- lowdens.multi[1:4]
# Make PCAs
# low:
pca.low <- prcomp(low.scaled, center=FALSE, scale=FALSE)
summary(pca.low)
# PC1 = 44.8%, PC2 = 30.8%
biplot(pca.low)
PC1low <- pca.low$x[,1]
PC2low <- pca.low$x[,2]
# PERMANOVA
mod.low <- adonis2(low.scaled ~ Temp, data = lowdens.multi, method = "euclidian") # PERMANOVA
mod.low
# sig effect of temp *
vec.low <- envfit(pca.low, lowdens.multi[4:10], perm = 1000) #fit physiological vectors onto ordination
vecdf.low <- as.data.frame(vec.low$vectors$arrows * sqrt(vec.low$vectors$r))
vecdf.low$variable <- rownames(vecdf.low)
# Test for differences in betadispersion by temp
bdisp.temp.low <- betadisper(vegdist(low.scaled, method = "euclidian"), fac.low$Temp)
anova(bdisp.temp.low)
# not significant
# Gather info for plotting
pca.df.low <- data.frame(pca.low$x[,1], y = pca.low$x[,2],
                         SymbStatus = as.factor(fac.low$SymbStatus),
                         Zub = as.factor(fac.low$Zub),
                         Temp = as.factor(fac.low$Temp))

# Plot
LowPCA <- ggbiplot(pca.low, ellipse = T, ellipse.prob = 0.95, alpha = 0, 
                   labels=rownames(pca.low), 
                   groups = pca.df.low$Temp, varname.adjust = 1, varname.size = 3.5) + 
  scale_color_manual(name = "Temp", values = c("#0571B0", "#92C5DE", "#F4A582", "#CA0020")) +
  xlab("PC1 (44.8%)") + ylab("PC2 (30.8%)") +
  ggtitle("Low Symbiont Density") +
  geom_vline(xintercept = 0, size = 0.1) + geom_hline(yintercept = 0, size = 0.1)+ 
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 12, color = "black"),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "right", legend.key = element_rect(fill = NA),
        strip.background = element_blank(), strip.text = element_blank()) + 
  coord_equal(ratio = 0.7) +
  coord_cartesian(xlim = c(-3.5,2.5))+
  geom_point(shape = 1, size = 2.5, aes(color = pca.df.low$Temp))
print(LowPCA)

####### High-density
Pmax <- read.csv("CC7 Pmax Fits.csv") %>% dplyr::select(c("anemone", "Pmax.umol.ug.hr"))
hidens.multi <- merge(highdens, Pmax)
# Dropping vials from HD that had values indistinguishable from blanks (below LOD)
hidens.drop.multi <- hidens.multi %>% subset(anemone != "11B" & anemone != "1B" & anemone != "2B" & anemone != "3B"  & anemone != "6B" & anemone != "5B" & anemone != "7B"& anemone != "8B"& anemone != "14B"& anemone != "23B"& anemone != "7A"& anemone != "6A")
Pmax.zub <- hidens.drop.multi %>% group_by(zub) %>% dplyr::summarise(mean(Pmax.umol.ug.hr))
hidens.drop.multi <- merge(hidens.drop.multi,Pmax.zub, by = "zub") 
# take only necessary columns
hidens.drop.multi <- dplyr::select(hidens.drop.multi, c("month","zub","target.temp.fac","symb.status",
                                                        "zav.symb.per.prot","zav.prot.ug.anemone","nonsymb.pHi.med","symb.pHi.med",
                                                        "host.13C.percent","symb.13C.percent","mean(Pmax.umol.ug.hr)","R.umol.ug.L.hr"))%>% drop_na()
colnames(hidens.drop.multi)[1:12] <- c("Cohort", "Zub", "Temp", "SymbStatus","Sym", "Protein", "Nonsymb pH", "Symbiocyte pHi", 
                                       "Host13C", "Symb13C", "Pmax","Respiration")
# Scale and center data
high.scaled <- scale(hidens.drop.multi[5:12], center = T, scale = T) 
# Identify factors 
fac <- hidens.drop.multi[1:4]
# Make PCAs
# high:
pca.hi <- prcomp(high.scaled, center=FALSE, scale=FALSE)
summary(pca.hi)
# PC1 = 42.7%, PC2 = 20.6%
biplot(pca.hi)
PC1hi <- pca.hi$x[,1]
PC2hi <- pca.hi$x[,2]
# PERMANOVA
mod.hi <- adonis2(high.scaled ~ Temp, data = hidens.drop.multi, method = "euclidian") # PERMANOVA
mod.hi
# sig effect of temp ***
vec.hi <- envfit(pca.hi, hidens.drop.multi[4:12], perm = 1000) #fit physiological vectors onto ordination
vecdf.hi <- as.data.frame(vec.hi$vectors$arrows * sqrt(vec.hi$vectors$r))
vecdf.hi$variable <- rownames(vecdf.hi)
# Test for differences in betadispersion by temp
bdisp.temp.hi <- betadisper(vegdist(high.scaled, method = "euclidian"), fac.high$Temp)
anova(bdisp.temp.hi)
# not significant
# Gather info for plotting
pca.df.hi <- data.frame(pca.hi$x[,1], y = pca.hi$x[,2],
                        SymbStatus = as.factor(fac.high$SymbStatus),
                        Zub = as.factor(fac.high$Zub),
                        Temp = as.factor(fac.high$Temp))
# Plot
HiDropPCA <- ggbiplot(pca.hi, ellipse = T, ellipse.prob = 0.95, size = 2.5,
                      labels=rownames(pca.hi), 
                      groups = pca.df.hi$Temp, varname.adjust = 1, varname.size = 3.5) + 
  scale_color_manual(name = "Temp", values = c("#0571B0", "#92C5DE", "#F4A582", "#CA0020")) +
  xlab("PC1 (42.7%)") + ylab("PC2 (20.6%)") +
  ggtitle("High Symbiont Density") + 
  geom_vline(xintercept = 0, size = 0.1) + geom_hline(yintercept = 0, size = 0.1)+ 
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 12, color = "black"),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "right", legend.key = element_rect(fill = NA),
        strip.background = element_blank(), strip.text = element_blank()) + 
  coord_equal(ratio = 0.7) +
  coord_cartesian(xlim = c(-2.5,3.5))
print(HiDropPCA)

# Both together for Figure S7
PCAs <- ggarrange(HiDropPCA,LowPCA,common.legend = T, widths = c(1,1), heights = c(1,1), align = "v")
PCAs

# All animals - apo and symb from both cohorts - 
# with just symb counts, chl, and nonsymbiocyte pH (the variables common to all)
trim <-  dplyr::select(pca, c("month","zub","target.temp.fac","symb.status",
                              "zav.symb.per.prot","zav.prot.ug.anemone","nonsymb.pHi.med")) %>% drop_na()
colnames(trim)[1:7] <- c("Cohort", "Zub", "Temp", "SymbStatus","Sym", "Prot", "pHi")
# Scale and center data
trim.scaled <- scale(trim[5:7], center = T, scale = T) 
# Identify factors 
fac.trim <- trim[1:4]
# Make PCA
pca.trim <- prcomp(trim.scaled, center=FALSE, scale=FALSE)
summary(pca.trim)
# PC1 = 45.3%, PC2 = 32.5%
biplot(pca.hi)
PC1hi <- pca.hi$x[,1]
PC2hi <- pca.hi$x[,2]
# PERMANOVA
mod.trim <- adonis2(trim.scaled ~ SymbStatus * Temp * Cohort, data = trim, method = "euclidian") # PERMANOVA
mod.trim
# Significant interactions: Cohort x SymbStatus (***), Temp x Cohort (*)
# Significant effects: Cohort, Temp, and SymbStatus (***)
mod.trim.notemp <- adonis2(trim.scaled ~ SymbStatus * Cohort, data = trim, method = "euclidian") # PERMANOVA
mod.trim.notemp
# ***S x C, ***S, ***C
# Fit physiological vectors onto ordination:
vec.trim <- envfit(pca.trim, trim[4:7], perm = 1000)
vecdf.trim <- as.data.frame(vec.trim$vectors$arrows * sqrt(vec.trim$vectors$r))
vecdf.trim$variable <- rownames(vecdf.trim)
# Test for differences in betadispersion by temp
bdisp.temp.trim <- betadisper(vegdist(trim.scaled, method = "euclidian"), fac.trim$Temp)
anova(bdisp.temp.trim)
# not significant
# by cohort:
bdisp.cohort.trim <- betadisper(vegdist(trim.scaled, method = "euclidian"), fac.trim$Cohort)
anova(bdisp.cohort.trim)
# not significant
bdisp.symstat.trim <- betadisper(vegdist(trim.scaled, method = "euclidian"), fac.trim$SymbStatus)
anova(bdisp.symstat.trim)
# yes, different betadispersion by symb status - violates normality assumption
# makes sense given much greater homogeneity of variance in aposymbiotic animals
# Gather info for plotting
pca.df.trim <- data.frame(pca.trim$x[,1], y = pca.trim$x[,2],
                          Cohort = as.factor(fac.trim$Cohort),
                          SymbStatus = as.factor(fac.trim$SymbStatus),
                          Zub = as.factor(fac.trim$Zub),
                          Temp = as.factor(fac.trim$Temp))

# Plot Figure S8 - grouping by symb status X cohort
CohortTrimPCABySymbStat <-   ggbiplot(pca.trim, ellipse = F,
                                      labels=rownames(pca.trim), 
                                      groups = interaction(pca.df.trim$Cohort, pca.df.trim$SymbStatus), 
                                      varname.adjust = 1, varname.size = 5)+
  stat_ellipse(level = 0.95, geom = "polygon", alpha = 0.4, aes(fill = interaction(pca.df.trim$Cohort, pca.df.trim$SymbStatus))) +
  scale_shape_manual(name = "Symbiont Status", values = c(2,2,1,16)) +
  scale_color_manual(name = "Symbiont Status", values = c("white", "grey", "darkgoldenrod4", "black")) +
  geom_point(aes(shape=interaction(pca.df.trim$Cohort,pca.df.trim$SymbStatus)), size = 2.4) +
  scale_fill_manual(name = "Symbiont Status", values = c("bisque2", "grey", "darkgoldenrod4", "black")) +
  xlab("PC1 (45.3%)") + ylab("PC2 (32.5%)") +
  geom_vline(xintercept = 0, size = 0.1) + geom_hline(yintercept = 0, size = 0.1)+ 
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 12, color = "black"),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "right", legend.key = element_rect(fill = NA),
        strip.background = element_blank(), strip.text = element_blank())
print(CohortTrimPCABySymbStat)