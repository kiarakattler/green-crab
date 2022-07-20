#Green crab functional responses based on Em's code

library(ggplot2)
library(lmerTest)
library(lme4)
library(dplyr)
library(tidyverse)
library(lamW)
library(LambertW)
library(frair)

crab <- read.csv("greencrabforaging.csv")
str(crab)

# male functional response ----------
mcrab <- filter(crab, Clam == "VC", Sex == "M", Sed == "no")

# Test for type II
frair_test(Clam_Consumed ~ Density, data = mcrab)

# Frair fit
outII_male <- frair_fit(Clam_Consumed ~ Density, data = mcrab, response = 'rogersII',
                     start = list(a = 0.2, h = 0.2), fixed = list(T=1))
# A linear fit
outI_male <- frair_fit(Clam_Consumed ~ Density, data = mcrab, response = 'typeI',
                    start = list(a = 0.2), fixed=list(T=1))

# Visualise fits
plot(outII_male, pch=20, col=rgb(0,0,0,0.2), xlim=c(0,16))
lines(outII_male)
lines(outI_male, lty=3)

# male resid
male_attack <- outII_male$coefficients[1] # Get coeffs
male_handling <- outII_male$coefficients[2]

fitsmale <- data.frame(x = mcrab$Density) # Calculate 'a' for each value of x, where x is density of clam

fitsmale$Ne <- fitsmale$x - lambertW0(male_attack * male_handling * fitsmale$x * exp(-male_attack * (1 - male_handling * fitsmale$x)))/(male_attack * male_handling) # calculate expected number of clam eaten, based on frair flexnr equation, using lambert function

fitsmale$actual <- mcrab$Clam_Consumed
fitsmale$resid <- fitsmale$Ne - fitsmale$actual

plot(x = fitsmale$x, y = fitsmale$Ne)
plot(x = fitsmale$Ne, y = fitsmale$resid)
abline(h = 0, lty = 'dotted')

# Have a look at original fits returned by mle2 (*highly* recommended)
summary(outII_male$fit)
summary(outI_male$fit)

# Compare models using AIC
AIC(outI_male$fit,outII_male$fit) #type II is for sure better


# Bootstrap
set.seed(309331)
outII_male_boot <- frair_boot(outII_male, start = NULL, strata=mcrab[,6], nboot=2000,
                           para=TRUE, ncores=NaN, WARN.ONLY=FALSE)

outII_male_boot
confint(outII_male_boot)

# Illustrate bootlines
plot(outII_male_boot, xlim=c(0,16), ylim = c(0, 16), type='n', main='All bootstrapped lines')
lines(outII_male_boot, all_lines=TRUE)
points(outII_male_boot, pch=20)

# Illustrate bootpolys
plot(outII_male_boot, xlim=c(0,16), ylim = c(0, 16), type='n', main='Empirical 95 percent CI')
drawpoly(outII_male_boot, col=hsv(2/6,0.2, 0.8))
points(outII_male_boot, pch=20)
lines(outII_male_boot, all_lines=FALSE)

# nlm to look at asymptote 
mcrab_asymp <- nls(Clam_Consumed ~ SSasymp(Density, Asym, R0, lrc), 
                   data = mcrab) # asymptotic curve with initial guesses
summary(mcrab_asymp)


# female functional response ----------
fcrab <- filter(crab, Clam == "VC", Sex == "F", Sed == "no")

# Test for type II
frair_test(Clam_Consumed ~ Density, data = fcrab)

# Frair fit
outII_female <- frair_fit(Clam_Consumed ~ Density, data = fcrab, response = 'rogersII',
                     start = list(a = 0.2, h = 0.2), fixed = list(T=1))
# A linear fit
outI_female <- frair_fit(Clam_Consumed ~ Density, data = fcrab, response = 'typeI',
                    start = list(a = 0.2), fixed=list(T=1))

# Visualise fits
plot(outII_male, pch=20, col=rgb(0,0,0,0.2), xlim=c(0,16))
lines(outII_female)
lines(outI_female, lty=3)

# female resid
female_attack <- outII_female$coefficients[1] # Get coeffs
female_handling <- outII_female$coefficients[2]

fitsfemale <- data.frame(x = fcrab$Density) # Calculate 'a' for each value of x, where x is density of clam

fitsfemale$Ne <- fitsfemale$x - lambertW0(female_attack * female_handling * fitsfemale$x * exp(-female_attack * (1 - female_handling * fitsfemale$x)))/(female_attack * female_handling) # calculate expected number of clam eaten, based on frair flexnr equation, using lambert function

fitsfemale$actual <- fcrab$Clam_Consumed
fitsfemale$resid <- fitsfemale$Ne - fitsfemale$actual

plot(x = fitsfemale$x, y = fitsfemale$Ne)
plot(x = fitsfemale$Ne, y = fitsfemale$resid)
abline(h = 0, lty = 'dotted')

# Have a look at original fits returned by mle2
summary(outII_female$fit)
summary(outI_female$fit)

# Compare models using AIC
AIC(outI_female$fit,outII_female$fit) #type II is for sure better

# Bootstrap
set.seed(309331)
outII_female_boot <- frair_boot(outII_female, start = NULL, strata=fcrab[,6], nboot=2000,
                           para=TRUE, ncores=NaN, WARN.ONLY=FALSE)

outII_female_boot
confint(outII_female_boot)

# Illustrate bootlines
plot(outII_female_boot, xlim=c(0,16), ylim = c(0, 16), type='n', main='All bootstrapped lines')
lines(outII_female_boot, all_lines=TRUE)
points(outII_female_boot, pch=20)

# Illustrate bootpolys
plot(outII_female_boot, xlim=c(0,16), ylim = c(0, 16), type='n', main='Empirical 95 percent CI')
drawpoly(outII_female_boot, col=hsv(2/6,0.2, 0.8))
points(outII_female_boot, pch=20)
lines(outII_female_boot, all_lines=FALSE)

# nlm to look at asymptote 
fcrab_asymp <- nls(Clam_Consumed ~ SSasymp(Density, Asym, R0, lrc), 
                   data = fcrab) # asymptotic curve with initial guesses
summary(fcrab_asymp)


#Functional response graph -----------------------
par(bg = 'white', fg = 'black')
plot(outII_female_boot, xlim=c(0, 16), ylim = c(0, 16), type='n',
     xlab = "Initial Clam Density",
     ylab="Clams Consumed", 
     cex.lab = 1.5,
     font.lab = 2,
     cex.axis = 1.2,
     cex.main = 1.5)
lines(outII_female_boot, lwd = 3, all_lines=FALSE, col= "#625a94", lty = 2)
lines(outII_male_boot, lwd = 3, all_lines=FALSE, col= "#11c2b5", lty = 1)
drawpoly(outII_female_boot, border = NA, col=adjustcolor("#625a94", alpha.f = 0.4))
drawpoly(outII_male_boot, border = NA, col=adjustcolor("#11c2b5", alpha.f = 0.4))
points(outII_female_boot, pch=17, col=adjustcolor("#625a94", alpha.f = 0.4), cex = 1.4)
points(outII_male_boot, pch=20,  col=adjustcolor("#11c2b5", alpha.f = 0.4), cex = 1.4)
legend(x = "topleft", legend = c("Male", "Female"), col = c("#11c2b5", "#625a94"), lty = c(1, 2), cex = 1.3, pch = c(20, 17))


# comparing asymptotes ----------

male.densities <- c(64)

bootest <- outII_male_boot$bootcoefs
maleboot <- male.densities - lambertW0(bootest[, 1] * bootest[, 2] * male.densities * exp(-bootest[, 1] * (1 - bootest[, 2] * male.densities)))/(bootest[, 1] * bootest[, 2]) 

female.densities <- c(64)

bootest <- outII_female_boot$bootcoefs
femaleboot <- female.densities - lambertW0(bootest[, 1] * bootest[, 2] * female.densities * exp(-bootest[, 1] * (1 - bootest[, 2] * female.densities)))/(bootest[, 1] * bootest[, 2]) 


bootdisttest <- data.frame(mcrab = maleboot, 
                           fcrab = femaleboot)

ggplot(data = bootdisttest)+
  geom_histogram(aes(x = fcrab), alpha = 0.3, fill = 'purple', colour = 'black', binwidth = 2)+
  geom_histogram(aes(x = mcrab), alpha = 0.3, fill = 'cyan', colour = 'black', binwidth = 2)+
  labs(x = 'Expected clams eaten with 16 starting density')

# Testing if there is a difference between the two asymptotes with kolmogorov-smirnov and t.tests
t.test(bootdisttest$fcrab, bootdisttest$mcrab)

stats::ks.test(bootdisttest$fcrab, bootdisttest$mcrab,
        exact = NULL)


# size comparisons --------------------------
#carapace size
mean_crab <- crab %>% 
  group_by(Sex) %>%
  mutate(Clam == "VC", Sed == "no") 
# %>% summarise(mean = mean(CW))

fmean_crab <- filter(mean_crab, Sex == "F")
mmean_crab <- filter(mean_crab, Sex == "M")

#check for normality
shapiro.test(fcrab$CW)
shapiro.test(mcrab$CW)
# can assume normality

#check to see if the variances are the same

cw.ftest <- var.test(CW ~ Sex, data = mean_crab)
cw.ftest
#p value >0.05 so we can assume the variances are equal

ttestcrab <- read_csv("greencrabforaging.csv")
ttestcrab1 <- filter(crab, Clam == "VC", Sed == "no")

ttest.cw <- t.test(CW ~ Sex, data = ttestcrab1, var.equal = TRUE)
ttest.cw
# CW not significantly different between males and females


# look at mean claw size -------
mean_crab_claw <- crab %>% 
  group_by(Sex) %>%
  mutate(Clam == "VC", Sed == "no") 

# %>% summarise(mean = mean(Claw_width))

fmean_crab_claw <- filter(mean_crab_claw, Sex == "F")
mmean_crab_claw <- filter(mean_crab_claw, Sex == "M")

#check for normality
shapiro.test(fmean_crab_claw$Claw_width)
shapiro.test(mmean_crab_claw$Claw_width)
# not normal??

#check to see if the variances are the same


cw.ftest <- var.test(CW ~ Sex, data = ttestcrab1)
cw.ftest

meanCW <- read_csv("meanCW.csv")

ttest.claw <- t.test(Claw_width ~ Sex, data = ttestcrab1, var.equal = TRUE)
ttest.claw
