# Aiptasia respirometry analysis & Pmax extraction
# for: Heat stress disrupts acid-base homeostasis independent of bleaching in a model cnidarian
# Allen-Waller L, KG Jones, M Martynek, KT Brown, and KL Barott
# compiled by Luella Allen-Waller
# 2023-04-27

library(ggplot2)
library(dplyr)
library(coda)
library(reshape)
library(reshape2)
library(scales)
library(zoo)
library(pipeR)
library(lubridate)
library(readr)
library(tidyverse)
library(devtools)
library(ggpubr)
library(patchwork)
library(tidyverse)
library(rstatix)
library(plotrix)
library(patchwork)
library(lme4)
library(lmerTest)
library(utf8)
library(mgcv)
library(emmeans)
library(MuMIn)

# Read in data
# Drop anemones whose respiration rates were below limit of detection (positive or 0 respiration values)
respdatanolow <- read.csv("AiptasiaHeatO2Evolution.csv") %>% 
  subset(anemone != "11B" & anemone != "1B" & anemone != "2B" & anemone != "3B"  & anemone != "6B" & anemone != "5B" &
           anemone != "7B"& anemone != "8B"& anemone != "14B"& anemone != "23B"& anemone != "7A"& anemone != "6A")
respdatanolow$rad <- as.numeric(respdata$light.level.PAR)
respdatanolow$normalized.O2.evol <- as.numeric(respdata$normalized.O2.evol)

# Read in factor data
key <- read.csv("respkey.csv") %>% unique() %>% drop_na(anem.num)
factorkey <- read.csv("MayCC7key.csv")

# merge in factor data
respdatanolow <- merge(respdatanolow,key)
respdatanolow <- merge(respdatanolow,factorkey)
respdatanolow$target.temp.fac <- as.factor(respdatanolow$target.temp)

# normalize blank-subtracted slopes ('normalized.O2.evol') to host protein
respdatanolow$protnormalized.O2.evol <- respdatanolow$normalized.O2.evol / respdatanolow$prot.ug.anemone

######### Photosynthesis-irradiance curves

### Plot raw O2 evolution data
# O2 evolution over light levels, accounting for 4mL chamber volume, per hour instead of per min,
# colored and faceted by temp, polynomial smoothed
ggplot(respdatanolow, mapping=aes(x = light.level.PAR, y = protnormalized.O2.evol*60*4/1000, 
                                  color = target.temp.fac, linetype = symb.status)) +
  geom_hline(yintercept = 0, size = 0.2) + 
  geom_point(aes(group = anem.num, shape = symb.status), size = 1.5, alpha = 0.5) +
  stat_smooth(method = "lm", formula = y ~ poly(x, 2), se = F, size = 0.5, alpha = 0.1, aes(group = anem.num)) +
  scale_color_manual("Target Temp",
                     values = c("25" = "#0571B0", "27" = "#92C5DE", "29" = "#F4A582", "31" = "#CA0020"))+
  scale_shape_manual("Symbiotic Status", values = c("apo"="triangle", "sym"="circle")) +
  scale_linetype_manual("Symbiotic Status", values = c("apo"="dashed", "sym"="solid")) +
  theme_bw()  + 
  theme(panel.grid = element_blank(), text = element_text(size = 14), legend.position = "none") +
  scale_x_continuous(breaks = c(0,200)) +
  labs(y = expression(paste(Oxygen~Evolution~(µmol~O["2"]~hr^-1~µg~anemone~protein^-1))), 
       x = expression(paste(Irradiance~(µmol~m^2~s^1))))+
  facet_grid(rows = vars(symb.status), cols = vars(target.temp.fac), scales = "free", space = "free") +
  scale_y_continuous(n.breaks = 4)
# This is Fig. S5

### Stats on photosynthesis-irradiance
# Test: does light level affect O2 evolution, and does this vary by temperature?
# splitting by symbiont status, since we expect apos and syms to respond differently

# Symbiotic anemones first:
respnolowsym <- respdatanolow %>% subset(symb.status == "sym")
symo2evol.lm <- lm(protnormalized.O2.evol~light.level.PAR, data = respnolowsym)
symo2evol.lm.temp <- lm(protnormalized.O2.evol~target.temp.fac*light.level.PAR, data = respnolowsym)
symo2evol.lme <- lmer(protnormalized.O2.evol~light.level.PAR + (1|zubID), data = respnolowsym)
symo2evol.lme.temp <- lmer(protnormalized.O2.evol~target.temp.fac*light.level.PAR+ (1|zubID), data = respnolowsym)
symo2.evol.gam <- gam(protnormalized.O2.evol~
                        s(light.level.PAR, bs="tp", k=3), data=respnolowsym) 
symo2.evol.gam.temp <- gam(protnormalized.O2.evol~
                             s(light.level.PAR, bs="tp", k=3), by = target.temp.fac, data=respnolowsym) 
symo2.evol.gam.outtemp <- gam(protnormalized.O2.evol ~
                                target.temp.fac +
                                s(light.level.PAR, bs="tp", k=3), data=respnolowsym)
AICc(symo2evol.lm,symo2evol.lm.temp,symo2evol.lme,symo2evol.lme.temp,symo2.evol.gam,symo2.evol.gam.temp,symo2.evol.gam.outtemp)
# symo2.evol.gam.outtemp is best
summary(symo2.evol.gam.outtemp) # model diagnostics & summary
par(mfrow=c(2,2))
gam.check(symo2.evol.gam.outtemp, pch=16)
# k-index < 1, but residual histogram looks OK
# Highly significant effect of light (***p<0.001) and temperature (*p<0.05)

# Aposymbiotic anemones:
respnolowapo <- respdatanolow %>% subset(symb.status == "apo")
apoo2evol.lm <- lm(protnormalized.O2.evol~light.level.PAR, data = respnolowapo)
apoo2evol.lm.temp <- lm(protnormalized.O2.evol~target.temp.fac*light.level.PAR, data = respnolowapo)
apoo2evol.lme <- lmer(protnormalized.O2.evol~light.level.PAR + (1|zubID), data = respnolowapo)
apoo2evol.lme.temp <- lmer(protnormalized.O2.evol~target.temp.fac*light.level.PAR+ (1|zubID), data = respnolowapo)
apoo2.evol.gam <- gam(protnormalized.O2.evol~
                        s(light.level.PAR, bs="tp", k=3), data=respnolowapo) 
apoo2.evol.gam.temp <- gam(protnormalized.O2.evol~
                             s(light.level.PAR, bs="tp", k=3), by = target.temp.fac, data=respnolowapo) 
apoo2.evol.gam.outtemp <- gam(protnormalized.O2.evol ~
                                target.temp.fac +
                                s(light.level.PAR, bs="tp", k=3), data=respnolowapo)
AICc(apoo2evol.lm,apoo2evol.lm.temp,symo2evol.lme,apoo2evol.lme,apoo2.evol.gam,apoo2.evol.gam.temp,apoo2.evol.gam.outtemp)
# apoo2.evol.gam and apoo2.evol.gam.temp perform equally well
# take the simpler one:
gam.check(apoo2.evol.gam, pch = 16)
# again, low k-index, but resids look OK
summary(apoo2.evol.gam)
# significant effect of light (***p<0.001)

######### Test: did temperature significantly affect Pmax?
### Extracting Pmax values for each symbiotic anemone

# get data and fit model for one anemone
a20 <- respdatanolow %>% subset(anemone == "20A")
a20fit <- gam(protnormalized.O2.evol ~ s(light.level.PAR, bs = "tp", k=3), data = a20)
a20fit
# Seems to work ok
fits <- respdatanolow %>%
  group_by(., anemone) %>%
  nest() %>%
  mutate(fit = purrr::map(data, ~mgcv::gam(protnormalized.O2.evol ~ s(light.level.PAR, bs = "tp", k=3), data = .x))) %>%
  mutate(augment_fit = map2(fit, data, ~broom::augment(.x, newdata = .y))) %>%
  unnest(augment_fit)
# Look at output object
dplyr::select(fits, anemone, data, fit)
# Get summary info
info <- fits %>%
  unnest_legacy(fit %>% purrr::map(broom::glance))
info
# Get parameters of GAM fits
params <- fits %>%
  unnest_legacy(fit %>% purrr::map(tidy))
# Subset out just the syms
syminfo <- info %>% subset(symb.status == "sym")
# Add lower and upper bounds of SE
syminfo$lower <- syminfo$.fitted - syminfo$.se.fit
syminfo$upper <- syminfo$.fitted + syminfo$.se.fit
# Plot:
SymAnemoneFits <- ggplot(syminfo, aes(x=light.level.PAR, color = target.temp.fac, group = anemone)) +
  facet_wrap(~target.temp.fac) +
  geom_point(aes(y=protnormalized.O2.evol)) +
  geom_line(aes(y=.fitted)) +
  geom_linerange(aes(ymin = lower, ymax = upper))+
  scale_color_manual("Target Temp",
                     values = c("25" = "#0571B0", "27" = "#92C5DE", "29" = "#F4A582", "31" = "#CA0020")) +
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 16, color = "black"),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "right", legend.key = element_rect(fill = NA)) +
  geom_hline(yintercept=0)
SymAnemoneFits
# ok, GAMs make sense at discrete intervals

# new data frame of predictions - do this to set a sequence to make a smooth curve with your prediction points
new_preds <- respdatanolow %>%  dplyr::do(., data.frame(light.level.PAR = seq(min(.$light.level.PAR), max(.$light.level.PAR), by = 1), stringsAsFactors = FALSE))
# the above sets a specific sequence so you can get a smooth curve from 1-269
# max and min for each curve
max_min <- respdatanolow %>% group_by(anemone) %>%
  dplyr::summarise(., min_light = min(light.level.PAR), max_light = max(light.level.PAR)) %>%
  ungroup()
# create new predictions
preds2 <- fits %>%
  unnest_legacy(fit %>% purrr::map(augment, newdata = new_preds)) %>%
  merge(., max_min, by = "anemone") %>%
  group_by(., anemone) %>%
  dplyr::select(-c("blank.avg.O2.evol", "light.level.PAR...6",".fitted...15", ".se.fit...16")) %>%
  filter(., light.level.PAR...17 > unique(min_light) & light.level.PAR...17 < unique(max_light)) %>%
  dplyr::rename(., pred.evol = .fitted...18) %>% dplyr::rename(., light.level.PAR = light.level.PAR...17) %>% dplyr::rename(., se.fit = .se.fit...19) %>%
  ungroup()
# dataframe now contains 270 predicted O2 evolution rates at increasing light levels for each anemone

# Check that those make sense by plotting over the raw data
sympreds2 <- preds2 %>% subset(symb.status == "sym")
SymAnemoneSmoothFits <- ggplot() +
  facet_wrap(~target.temp.fac) +
  geom_point(data = syminfo, aes(y=protnormalized.O2.evol, x = light.level.PAR, color = target.temp.fac, group = anemone)) +
  geom_line(data = sympreds2, aes(y=pred.evol, x = light.level.PAR, color = target.temp.fac, group = anemone)) +
  scale_color_manual("Target Temp",
                     values = c("25" = "#0571B0", "27" = "#92C5DE", "29" = "#F4A582", "31" = "#CA0020")) +
  theme(axis.title = element_text(size = 16), axis.text = element_text(size = 16, color = "black"),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "right", legend.key = element_rect(fill = NA)) +
  geom_hline(yintercept=0)
SymAnemoneSmoothFits

# Now pull the maximum pred.evol from each anemone and send to a new DF
fit_Pmax <- preds2 %>% group_by(anemone) %>% dplyr::summarise(Pmax.umol.min.ug.L = max(pred.evol))
# Put factor data back in using the keys from earlier
fit_Pmax <- merge(fit_Pmax,key) 
fit_Pmax <- merge(fit_Pmax,factorkey) %>% dplyr::select(-c("anem.num","prot.ug.anemone"))
fit_Pmax$Pmax.umol.ug.hr <- fit_Pmax$Pmax*60*4/1000                                              
# Check:
fit_Pmax
# Save:
write.csv(fit_Pmax, "Aiptasia Heat Pmax Fits.csv")

############ How did Pmax and LEDR respond to temperature?

# Bring back light-enhanced dark respiration from earlier dataframe
LEDR <- respdatanolow %>% subset(light.level.PAR==0) %>% dplyr::select(c(anemone, zubID, symb.status, protnormalized.O2.evol))
colnames(LEDR)[colnames(LEDR) == "protnormalized.O2.evol"] ="LEDR.umol.L.ug.min"
LEDR$LEDR.umol.ug.hr <- LEDR$LEDR.umol.L.ug.min*60*4/1000*(-1)
# Dataframe with both Pmax and LEDR for each anemone
highchronicTPC <- merge(fit_Pmax,LEDR)
highchronicTPC$temp <- as.numeric(highchronicTPC$target.temp)

### Testing how symbiotic metabolism changed with temperature
# Maximum photosynthetic rate:
symTPC <- highchronicTPC %>% subset(symb.status == "sym")
symbGP.TPC.lm <- lm(Pmax.umol.ug.hr ~ temp, data = symTPC)
symbGP.TPC.lm2 <- lm(Pmax.umol.ug.hr ~ poly(temp, 2, raw = T), data = symTPC)
symGP.TPC.lmer <-  lmer(Pmax.umol.ug.hr ~ temp + (1|zubID), data = symTPC)
symGP.TPC.lmer2 <-  lmer(Pmax.umol.ug.hr ~ poly(temp,2,raw = T) + (1|zubID), data = symTPC)
symGP.gam <- gam(Pmax.umol.ug.hr~s(temp, bs="tp", k=3), data=symTPC) 
symGP.gam2 <- gam(Pmax.umol.ug.hr~poly(temp, 2, raw = T), data=symTPC) 
AICc(symbGP.TPC.lm,symbGP.TPC.lm2,symGP.TPC.lmer,symGP.TPC.lmer2,symGP.gam,symGP.gam2)
# symGP.gam is the best
symGP.gam3 <- gam(Pmax.umol.ug.hr~s(temp, bs="tp", k=3), data=symTPC) 
symGP.gam4 <- gam(Pmax.umol.ug.hr~s(temp, bs="tp", k=3) + s(zubID, bs="re", k=3), data=symTPC) 
symGP.gam5 <- gam(Pmax.umol.ug.hr~s(temp, bs="tp", k=4) + s(zubID, bs="re", k=3), data=symTPC) 
symGP.gam6 <- gam(Pmax.umol.ug.hr~s(temp, bs="tp", k=4), data=symTPC) 
AICc(symGP.gam6,symGP.gam5,symGP.gam4,symGP.gam3,symGP.gam)
# symGP.gam5 performs best of these 
summary(symGP.gam5) # model diagnostics & summary
# Temperature significantly affected photosynthesis: F=4.218, p=0.0212, R2-adj = 0.25
par(mfrow=c(2,2))
gam.check(symGP.gam5, pch=16)
# k-index > 1

# Plot Pmax with that model:
Pmax_TPC_gam <- ggplot(data = symTPC, aes(x = temp, y = Pmax.umol.ug.hr)) +
  geom_hline(yintercept=0) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "tp", k=4), se = T, 
              color = "black", size = 0.5, alpha = 0.3) +
  ggtitle("Gross photosynthesis") +
  labs(y = expression(atop("Gross Photosynthesis", paste((µmol~O["2"]~hr^-1~anemone^-1)))), x = expression(paste("Temperature (ºC)"))) +
  scale_x_continuous(n.breaks = 6) +
  theme_bw()  + 
  theme(panel.grid = element_blank(), text = element_text(size = 20)) +
  geom_point(size = 2.5, alpha = 0.6) +
  coord_cartesian(ylim =c(-0.002,0.035)) +
  labs(title = "")
Pmax_TPC_gam
# What is the optimal predicted temperature for Pmax?
# Calculate inflection points using function from Gavin Simpson (June 2014)
script <- getURL("https://gist.githubusercontent.com/gavinsimpson/ca18c9c789ef5237dbc6/raw/095c9be4d3654b5a8c05aaa6b9037ad1bdab53b3/derivSimulCI.R", ssl.verifypeer = FALSE)
eval(parse(text = script))
# for GP gam
# assess the slope by calculating the first derivative of the GAM spline
fd <- derivSimulCI(symGP.gam3)
CI <- lapply(fd[1], function(temp) t(apply(temp$simulations, 1, quantile, probs = c(0.025, 0.975))))
first.zero.slope.index <- min(which(sign(CI$temp[, "2.5%"]) != sign(CI$temp[, "97.5%"])))
fd$eval[first.zero.slope.index]
# output = 25
# this found the first place at which there is no significant difference between the fit and 0 
# (lower bound of CI)
# what about the midpoint of the CI (approaching actual optimum)?
infpt <- lapply(fd[1], function(temp) t(apply(temp$simulations, 1, quantile, probs = c(0.49999, 0.50009))))
zero.slope.index <- min(which(sign(infpt$temp[, "49.999%"]) != sign(infpt$temp[, "50.009%"])))
fd$eval[zero.slope.index] # 49.9% and 49.99% give 26.38693 ºC, but going to 49.999% puts it up to infinity -> too far
26.38693 - 25
# The predicted optimal temperature for Pmax was 26.39 ± 1.39 ºC

# Respiration:
symR.TPC.lm <- lm(LEDR.umol.ug.hr ~ temp, data = symTPC)
symR.TPC.lm2 <- lm(LEDR.umol.ug.hr ~ poly(temp, 2, raw = T), data = symTPC)
symR.TPC.lmer <-  lmer(LEDR.umol.ug.hr ~ temp + (1|zubID), data = symTPC)
symR.TPC.lmer2 <-  lmer(LEDR.umol.ug.hr ~ poly(temp,2,raw = T) + (1|zubID), data = symTPC)
symR.gam <- gam(LEDR.umol.ug.hr~s(temp, bs="tp", k=3), data=symTPC) 
symR.gam2 <- gam(LEDR.umol.ug.hr~poly(temp, 2, raw = T), data=symTPC) 
symR.gam3 <- gam(LEDR.umol.ug.hr~s(temp, bs="tp", k=4), data=symTPC) 
symR.gam4 <- gam(LEDR.umol.ug.hr~s(temp, bs="tp", k=3) + s(zubID, bs="re", k=3), data=symTPC) 
symR.gam5 <- gam(LEDR.umol.ug.hr~s(temp, bs="tp", k=4) + s(zubID, bs="re", k=3), data=symTPC) 
AICc(symR.TPC.lm,symR.TPC.lm2,symR.TPC.lmer,symR.TPC.lmer2,symR.gam,symR.gam2,symR.gam3,symR.gam4,symR.gam5)
# symR.gam3 is the best
summary(symR.gam3)# model diagnostics & summary
par(mfrow=c(2,2))
gam.check(symR.gam3, pch=16)
# k-index = 0.94
# Temperature significantly affected respiration: F=4.46, P=0.030, R2-adj = 0.24
# Plot with that model:
LEDR_TPC_gam <- ggplot(data = symTPC, aes(x = temp, y = LEDR.umol.ug.hr)) +
  geom_hline(yintercept=0) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "tp", k=4), se = T, 
              color = "black", size = 0.5, alpha = 0.3) +
  ggtitle("Respiration") +
  labs(y = expression(atop("Light-Enhanced Dark Respiration", paste((µmol~O["2"]~hr^-1~anemone^-1)))), x = expression(paste("Temperature (ºC)"))) +
  scale_x_continuous(n.breaks = 6) +
  theme_bw()  + 
  theme(panel.grid = element_blank(), text = element_text(size = 20)) +
  geom_point(size = 2.5, alpha = 0.6) +
  coord_cartesian(ylim =c(-0.002,0.035)) +
  labs(title = "")
LEDR_TPC_gam
# What is the optimal predicted temperature for LEDR?
# assess the slope by calculating the first derivative of the GAM spline
fdR <- derivSimulCI(symR.gam3)
CI <- lapply(fdR[1], function(temp) t(apply(temp$simulations, 1, quantile, probs = c(0.025, 0.975))))
first.zero.slope.index <- min(which(sign(CI$temp[, "2.5%"]) != sign(CI$temp[, "97.5%"])))
fdR$eval[first.zero.slope.index]
# output = 26.65829
# this found the first place at which there is no significant difference between the fit and 0 
# (lower bound of CI)
# midpoint of the CI (approaching actual optimum):
infpt <- lapply(fdR[1], function(temp) t(apply(temp$simulations, 1, quantile, probs = c(0.49, 0.51))))
zero.slope.index <- min(which(sign(infpt$temp[, "49%"]) != sign(infpt$temp[, "51%"])))
fdR$eval[zero.slope.index] # 48% and 49% give 27.23116 ºC, but going to 49.1% puts it up to infinity -> too far
27.23116 - 26.65829
# The predicted optimal temperature for Pmax was 27.23 ± 0.57 ºC

### Testing whether aposymbiotic metabolism changed with temperature
aponolow <- highchronicTPC %>% subset(symb.status == "apo")

# Apo Pmax:
apoGP.TPC.lm <- lm(Pmax.umol.ug.hr ~ temp, data = aponolow)
apoGP.TPC.lm2 <- lm(Pmax.umol.ug.hr ~ poly(temp, 2, raw = T), data = aponolow)
apoGP.TPC.lmer <-  lmer(Pmax.umol.ug.hr ~ temp + (1|zubID), data = aponolow)
apoGP.TPC.lmer2 <-  lmer(Pmax.umol.ug.hr ~ poly(temp,2,raw = T) + (1|zubID), data = aponolow)
apoGP.gam <- gam(Pmax.umol.ug.hr~s(temp, bs="tp", k=3), data=aponolow) 
apoGP.gam2 <- gam(Pmax.umol.ug.hr~poly(temp, 2, raw = T), data=aponolow) 
AICc(apoGP.TPC.lm,apoGP.TPC.lm2,apoGP.TPC.lmer,apoGP.TPC.lmer2,apoGP.gam,apoGP.gam2)
# regular gam is the best
apoGP.gam3 <- gam(Pmax.umol.ug.hr~s(temp, bs="tp", k=4), data=aponolow) 
apoGP.gam4 <- gam(Pmax.umol.ug.hr~s(temp, bs="tp", k=3) + s(zubID, bs="re", k=3), data=aponolow) 
apoGP.gam5 <- gam(Pmax.umol.ug.hr~s(temp, bs="tp", k=4) + s(zubID, bs="re", k=3), data=aponolow) 
AICc(apoGP.gam5,apoGP.gam4,apoGP.gam3,apoGP.gam, apoGP.TPC.lm)
# apoGP.gam5 best
gam.check(apoGP.gam5, pch = 16)
# k-indices too low
gam.check(apoGP.gam3, pch = 16)
# k-indices still too low even with simpler gam model
# probably more appropriate to do a linear model, which performs almost as well as GAM anyway
anova(apoGP.TPC.lm)
# No significant effect of temperature on Pmax in apos: F=3.1869, P=0.0902
# Plot with that model:
apo_Pmax_TPC_lm <- ggplot(data = aponolow, aes(x = temp, y = Pmax.umol.ug.hr)) +
  geom_hline(yintercept=0) +
  geom_smooth(method = "lm", formula = y ~ x, se = T, 
              color = "black", size = 0.5, alpha = 0.3) +
  ggtitle("Gross photosynthesis") +
  labs(y = expression(atop("Gross Photosynthesis", paste((µmol~O["2"]~hr^-1~anemone^-1)))), x = expression(paste("Temperature (ºC)"))) +
  scale_x_continuous(n.breaks = 6) +
  theme_bw()  + 
  theme(panel.grid = element_blank(), text = element_text(size = 20)) +
  geom_point(size = 2.5, alpha = 0.6) +
  coord_cartesian(ylim =c(-0.003,0.035)) +
  labs(title = "")
apo_Pmax_TPC_lm

# apo LEDR:
apoR.TPC.lm <- lm(LEDR.umol.ug.hr ~ temp, data = aponolow)
apoR.TPC.lm2 <- lm(LEDR.umol.ug.hr ~ poly(temp, 2, raw = T), data = aponolow)
apoR.TPC.lmer <-  lmer(LEDR.umol.ug.hr ~ temp + (1|zubID), data = aponolow)
apoR.TPC.lmer2 <-  lmer(LEDR.umol.ug.hr ~ poly(temp,2,raw = T) + (1|zubID), data = aponolow)
apoR.gam <- gam(LEDR.umol.ug.hr~s(temp, bs="tp", k=3), data=aponolow) 
apoR.gam2 <- gam(LEDR.umol.ug.hr~poly(temp, 2, raw = T), data=aponolow) 
AICc(apoR.TPC.lm,apoR.TPC.lm2,apoR.TPC.lmer,apoR.TPC.lmer2,apoR.gam,apoR.gam2)
# apoR.TPC.lm and the simplest gam are the best
apoR.gam3 <- gam(LEDR.umol.ug.hr~s(temp, bs="tp", k=4), data=aponolow) 
apoR.gam4 <- gam(LEDR.umol.ug.hr~s(temp, bs="tp", k=3) + s(zubID, bs="re", k=3), data=aponolow) 
apoR.gam5 <- gam(LEDR.umol.ug.hr~s(temp, bs="tp", k=4) + s(zubID, bs="re", k=3), data=aponolow) 
AICc(apoR.TPC.lm,apoR.gam,apoR.gam3,apoR.gam4,apoR.gam5)
# All the same so let's just stick with simple linear model
anova(apoR.TPC.lm)
# Temperature affected aposymbiotic anemone respiration: f=14.223, p=0.0013
# Plot with LM:
apo_LEDR_TPC_lm <- ggplot(data = aponolow, aes(x = temp, y = LEDR.umol.ug.hr)) +
  geom_hline(yintercept=0) +
  geom_smooth(method = "lm", formula = y~x, se = T, 
              color = "black", size = 0.5, alpha = 0.3) +
  ggtitle("Respiration") +
  labs(y = expression(atop("Light-Enhanced Dark Respiration", paste((µmol~O["2"]~hr^-1~anemone^-1)))), x = expression(paste("Temperature (ºC)"))) +
  scale_x_continuous(n.breaks = 6) +
  theme_bw()  + 
  theme(panel.grid = element_blank(), text = element_text(size = 20)) +
  geom_point(size = 2.5, alpha = 0.6) +
  coord_cartesian(ylim =c(-0.003,0.035)) +
  labs(title = "")
apo_LEDR_TPC_lm

### Final metabolic rate plots

# Figure 5: symbiotic metabolic rates
SymTPCs <- ggarrange(Pmax_TPC_gam, LEDR_TPC_gam)
SymTPCs

# Figure S6: aposymbiotic metabolic rates
ApoTPC <- ggarrange(apo_Pmax_TPC_lm, apo_LEDR_TPC_lm)
ApoTPC
