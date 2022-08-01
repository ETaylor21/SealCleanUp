### drift corrections for SEAL data for NOx
library(tidyverse)
library(ggplot2)
library(dplyr)
library(stringr)
library(chron)
library(lubridate)
library(data.table)

##bring in the raw data from a run
rd <- read.csv("20-09-24_ET_NOx_Aug.csv") ## change here specific for your data

head(rd)
summary(rd)
## break apart by type of nutrient analyzed if needed
## this is redundant if you only have NOx data but if you ran OP or NH4 at the same time it is needed
nitr.raw <- subset(rd, Test=="NO3_NO2")

############################ nitrate ##########################################
## now I need to make a column that has the type of data (unknown vs CCV vs Duplicate vs Standard ect..)
levels(as.factor(nitr.raw$Sample.ID)) ## this is what we have 
nitr.raw <- nitr.raw %>% mutate(Code1 = substr(Sample.ID,1L,3L))
nitr.raw$Code1 <-as.factor(nitr.raw$Code1)
levels(nitr.raw$Code1)[levels(nitr.raw$Code1)=="NO3"] <- 'Sta' #make the N03 to Sta
levels(nitr.raw$Code1)[levels(nitr.raw$Code1)=="Nit"] <- 'Sta' 

keep <- as.factor(c("CCB", "CCV", "Dup", "Sta")) ## create the unknowns
nitr.raw$Code2 <- as.character(ifelse(nitr.raw$Code1 %in% keep, nitr.raw$Code1, "Unk"))

nitr.raw$Code1<- as.character(nitr.raw$Code1)## merge together
nitr.raw$sample.type <- ifelse(nitr.raw$Code2=="Unk", "Unk", nitr.raw$Code1)

names(nitr.raw) ## remove the columns we dont need
nitr <- nitr.raw[ , -which(names(nitr.raw) %in% c("Code1","Code2"))]

################################################################################
## now we are ready to do the drift correction
################################################################################
################ NOx correction of drifts ######################################

## standard curves
stands.nox <- subset(nitr, sample.type=="Sta")
reds <- c("Nitrate Standard", "Nitrite Standard")
red.eff.sta <- stands.nox[stands.nox$Sample.ID %in% reds,] ## take out the nitrate nitrite to use later
stands.nox <- stands.nox[stands.nox$Sample.ID != reds,] ## remove the nitrate nitrite from the standard curve
stands.nox$exp_conc <- c(0, 0.012, 0.02, 0.052, 0.1, 0.2, 0.5, 1,2, 0)
stands.nox$std.type <- c(rep("low.std", 7), rep("clean.std", 2), "low.sd")


## all standards
(all.std.range.n <- range(stands.nox$exp_conc)) ## the range of concentration of standard we will use
all.min.abs.n <- min(stands.nox$Absorbance)
all.max.abs.n <- max(stands.nox$Absorbance)
all.std.curve.n <- lm(exp_conc~Absorbance, data=stands.nox)
all.std.curve.intercept.n <- coef(all.std.curve.n)[1]
all.std.curve.slope.n <- coef(all.std.curve.n)[2]
summary(all.std.curve.n)$r.squared 
coef(all.std.curve.n)
confint(all.std.curve.n, level = 0.95) 

ggplot(data=stands.nox, aes(Absorbance, exp_conc)) + geom_point()+ 
  theme_light() + geom_smooth(method="lm", lty=3) + 
  ggtitle("All standards")+ ylab("Expected Concentration")

## cleaned up standards ######################################################
(clean.std.range.n <- range(stands.nox[stands.nox$exp_conc>0,]$exp_conc)) ## the range of concentration we will use
clean.min.abs.n <- min(stands.nox[stands.nox$exp_conc>0,]$Absorbance)
clean.max.abs.n <- max(stands.nox[stands.nox$exp_conc>0,]$Absorbance)
clean.std.curve.n <- lm(exp_conc~Absorbance, data=stands.nox[stands.nox$exp_conc>0,])
confint(clean.std.curve.n, level = 0.95)  ## shows the 95% CIs of our fit
clean.std.curve.intercept.n <- coef(clean.std.curve.n)[1]
clean.std.curve.slope.n <- coef(clean.std.curve.n)[2]
summary(clean.std.curve.n)$r.squared 
coef(clean.std.curve.n)
confint(clean.std.curve.n, level = 0.95) 


ggplot(data=stands.nox[stands.nox$exp_conc>0,], aes(Absorbance, exp_conc)) + geom_point()+ 
  theme_light() + geom_smooth(method="lm", lty=3) + 
  ggtitle("Cleaned up standards")+ ylab("Expected Concentration")

######## standards in the lower range ########################################
low.stands.nox <-stands.nox[stands.nox$exp_conc >0 & stands.nox$exp_conc <0.2, ]
(low.std.range.n <- range(low.stands.nox$exp_conc)) ## the range of standard we will use
low.min.abs.n <- min(low.stands.nox$Absorbance)
low.max.abs.n <- max(low.stands.nox$Absorbance)
low.std.curve.n <- lm(exp_conc~Absorbance, data=low.stands.nox)
coef(low.std.curve.n)
confint(low.std.curve.n, level = 0.95)  ## shows the 95% CIs of our fit
low.std.curve.intercept.n <- coef(low.std.curve.n)[1]
low.std.curve.slope.n <- coef(low.std.curve.n)[2]
summary(low.std.curve.n)$r.squared 

ggplot(data=low.stands.nox, aes(Absorbance, exp_conc)) + geom_point()+ 
  theme_light() + geom_smooth(method="lm", lty=3) + 
  ggtitle("Low standards")+ ylab("Expected Concentration")

###############################################################################
## CALIBRATION STANDARDS ######################################################
names(stands.nox)
stands.nox$calc_conc <- ifelse(stands.nox$Absorbance < low.max.abs.n,
                               low.std.curve.slope.n * stands.nox$Absorbance + low.std.curve.intercept.n,
                               clean.std.curve.slope.n * stands.nox$Absorbance + clean.std.curve.intercept.n)*stands.nox$Auto.Dil.Factor

stands.nox$calc_per_error <- ((stands.nox$calc_conc - stands.nox$exp_conc) /  stands.nox$exp_conc) *100

#############################################################################
## reduction efficency check

red.eff.sta$calc_conc <- ifelse(red.eff.sta$Absorbance < low.max.abs.n,
                               low.std.curve.slope.n * red.eff.sta$Absorbance + low.std.curve.intercept.n,
                               clean.std.curve.slope.n * red.eff.sta$Absorbance + clean.std.curve.intercept.n)*red.eff.sta$Auto.Dil.Factor

NO3.reduction.efficency <- (red.eff.sta$calc_conc[1]/red.eff.sta$calc_conc[2])*100; NO3.reduction.efficency


##############################################################################
### CHECK STANDARDS ##########################################################
stands.ccv.n<- subset(nitr, sample.type=="CCV")
names(stands.ccv.n)
stands.ccv.n$exp_conc <- rep(0.5, nrow(stands.ccv.n))
stands.ccv.n$std.type <- c(rep("low.std", nrow(stands.ccv.n)))

stands.ccv.n$calc_conc <- ifelse(stands.ccv.n$Absorbance < low.max.abs.n,
                                 low.std.curve.slope.n * stands.ccv.n$Absorbance + low.std.curve.intercept.n,
                                 clean.std.curve.slope.n * stands.ccv.n$Absorbance + clean.std.curve.intercept.n)*stands.ccv.n$Auto.Dil.Factor

stands.ccv.n$calc_per_error <- ((stands.ccv.n$calc_conc - stands.ccv.n$exp_conc) /  stands.ccv.n$exp_conc) *100

### I need to pull out the times of the measurements and convert it to numeric to be able to fit model
#### be careful here and check how your date and time look
stands.ccv.n$Date.and.Time ## you need to tell it where to start and stop cutting
#we want to cut out just the time, which starts at space 12 and end at space 20 in the string or 10 14
#stands.ccv.n <- stands.ccv.n %>% mutate(time = paste(substr(Date.and.Time,10,14), sep=""))
stands.ccv.n <- stands.ccv.n %>% mutate(time = paste(substr(Date.and.Time,12,20), sep=""))
stands.ccv.n$Time <- chron(times= as.character(stands.ccv.n$time)) ## this makes r see the time as a real time and not a charachter
#stands.ccv.n$Time <- chron(times= as.character(stands.ccv.n$time), times = "h:m") 
stands.ccv.n$seconds <- as.numeric(hms(stands.ccv.n$Time))## then to make the regression i convert to seconds
?hms
ccv.curve.n <- lm(calc_conc~seconds, data=stands.ccv.n)
coef(ccv.curve.n)
ccv.curve.intercept.n <- coef(ccv.curve.n)[1]
cc.curve.slope.n <- coef(ccv.curve.n)[2]

ggplot(data=stands.ccv.n, aes(x=seconds, y=calc_conc)) + geom_point()+ 
  theme_light() + geom_smooth(method="lm", lty=3) + 
  ggtitle("Drift Check")+ ylab("Calculated Concentration") + xlab("run time")

stands.ccv.n$expected.ccv <- cc.curve.slope.n * stands.ccv.n$seconds + ccv.curve.intercept.n

stands.ccv.n$meas_minus_exp <- stands.ccv.n$calc_conc[1] - stands.ccv.n$expected.ccv 

stands.ccv.n$drift_corrected_ccv <- stands.ccv.n$calc_conc + stands.ccv.n$meas_minus_exp

ggplot(data=stands.ccv.n, aes(x=seconds, y=drift_corrected_ccv)) + geom_point()+ 
  theme_light() + geom_smooth(method="lm", lty=3) + 
  ggtitle("Drift Check")+ ylab("Drift Corrected Concentration") + xlab("run time")

###########################################################################################
##### now to correct our unknown data and the Duplicates 

unknowns.n <- nitr[nitr$sample.type == "Unk" | nitr$sample.type == "Dup", ]
## us the same notation with the time as above
unknowns.n<- unknowns.n %>% mutate(time = paste(substr(Date.and.Time,12,20), sep=""))
unknowns.n$Time <- chron(times= as.character(unknowns.n$time)) 
unknowns.n$seconds <- as.numeric(hms(unknowns.n$Time))


unknowns.n$calc_conc <- ifelse(unknowns.n$Absorbance < low.max.abs.n,
                               low.std.curve.slope.n * unknowns.n$Absorbance + low.std.curve.intercept.n,
                               clean.std.curve.slope.n * unknowns.n$Absorbance + clean.std.curve.intercept.n)* unknowns.n$Auto.Dil.Factor


ccv.drift_curve.n <- lm(meas_minus_exp~seconds, data=stands.ccv.n)
coef(ccv.drift_curve.n)
ccv.drift.curve.intercept.n <- coef(ccv.drift_curve.n)[1]
cc.drift.curve.slope.n <- coef(ccv.drift_curve.n)[2]


unknowns.n$total.drift <- cc.drift.curve.slope.n * unknowns.n$seconds + ccv.drift.curve.intercept.n 

unknowns.n$actual.conc <- ifelse(unknowns.n$calc_conc > 0.5 * stands.nox$exp_conc[2],
                                 unknowns.n$calc_conc + unknowns.n$total.drift, 0)


### Now I want to make some indentifyer columns and get the data ready for export and
## make some graphs to check the data

## first I will split out the duplicates and we will merge later
## this is because the id variables are in different columns
dups.n <- subset(unknowns.n, sample.type == "Dup")
other.n <- subset(unknowns.n, sample.type == "Unk")

## using library unglue we will break about the sample.id
?unglue_vec ## you can set up the pattern of your data and take out any part
other.n$Sample.ID ## there are three sections {x=site name}_{y=date}_{z=rep} separated by _
other.n$site <- unglue_vec(other.n$Sample.ID,"{x}_{y}_{z}", var = "x")
other.n$sample_date <- unglue_vec(other.n$Sample.ID,"{x}_{y}_{z}", var = "y")
other.n$rep <- unglue_vec(other.n$Sample.ID,"{x}_{y}_{z}", var = "z")

head(dups.n)## here the id is in the sample details
dups.n$site <- unglue_vec(dups.n$Sample.Details,"{x}_{y}_{z}", var = "x")
dups.n$sample_date <- unglue_vec(dups.n$Sample.Details,"{x}_{y}_{z}", var = "y")
dups.n$rep <- unglue_vec(dups.n$Sample.Details,"{x}_{y}_{z}", var = "z")

final.nox <- rbind(other.n, dups.n)
head(final.nox)
final.nox$site <- as.factor(final.nox$site)
levels(final.nox$site)
final.nox$site <- factor(final.nox$site, levels = c( "HAT",   "HOGDN", "HOGUP", "POS"  , "SWBUP" ,  "SWB" ,"TUM" )) ## put in order of urbanization


ggplot(data=final.nox, aes(x=site, y=actual.conc, fill=site)) +
  geom_boxplot() + theme_bw() + geom_jitter(width=0.2, shape=4) + 
  ylab("nox (mg N/L)") + facet_wrap(~sample_date)

## to see how the duplicates line up with their replicate
selected <- levels(as.factor(dups.n$Sample.Details)); selected
levels(as.factor(final.nox$Sample.ID))
red_nh4 <- final.nox[final.nox$Sample.ID %in% selected,]
red_nh42 <- rbind(dups.n, red_nh4)

ggplot(data=red_nh42, aes(x=site, y=actual.conc, color=sample.type)) +
  geom_jitter(data=red_nh42, aes(x=site, y=actual.conc),width=0.3, size=3) + theme_bw() +  
  ylab("nox (mg N/L)") + facet_wrap(~sample_date)

## finally export data in a standard way "nutrient_YYYYMMDD"
final.nox$sample_date
## if there are two or more dates you can subset by date and export individually
# example  final.nox.a <- subset(final.nox, sample_date == "01102020")
write.csv(final.nox, "nox_20200110.csv")



