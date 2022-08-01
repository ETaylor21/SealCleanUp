
### drift corrections for SEAL data
library(tidyverse)
library(ggplot2)
library(dplyr)
library(stringr)
library(chron)
library(lubridate)
library(data.table)
library(unglue)

##bring in the raw data from a run
rd <- read.csv("C:/Users/emily/University of Florida/IFAS-SWS-ReisingerLab - Gainesville Streams/Data/Seal (AQ400) Data/SEAL RAW/20-09-21_ET_NH4OP_JulAug.csv")

head(rd)
summary(rd)
## break apart by type of nutrient analyzed
levels(as.factor(rd$Test))
amon.raw <- subset(rd, Test=="NH4-2mg_L")
phos.raw <- subset(rd, Test=="o-PHOS")

## now I need to make a column that has the type of data #####################################
### this was done in exel previously but I want to skip this step and have it done here ######
################### Ammonium #################################################################
levels(as.factor(amon.raw$Sample.ID)) ## this is what we have now
amon.raw <- amon.raw %>% mutate(Code1 = substr(Sample.ID,1L,3L))
amon.raw$Code1 <-as.factor(amon.raw$Code1)
levels(amon.raw$Code1)[levels(amon.raw$Code1)=="NH4"] <- 'Sta' #make the NH4 to Sta

keep <- as.factor(c("CCB", "CCV", "Dup", "Sta")) ## create the unknowns.n
amon.raw$Code2 <- as.character(ifelse(amon.raw$Code1 %in% keep, amon.raw$Code1, "Unk"))

amon.raw$Code1<- as.character(amon.raw$Code1)## merge together
amon.raw$sample.type <- ifelse(amon.raw$Code2=="Unk", "Unk", amon.raw$Code1)

names(amon.raw) ## remove the columns we dont need
amon <- amon.raw[ , -which(names(amon.raw) %in% c("Code1","Code2"))]


############################ phosphate #########################################
levels(as.factor(phos.raw$Sample.ID)) ## this is what we have now
phos.raw <- phos.raw %>% mutate(Code1 = substr(Sample.ID,1L,3L))
phos.raw$Code1 <-as.factor(phos.raw$Code1)
levels(phos.raw$Code1)[levels(phos.raw$Code1)=="PO4"] <- 'Sta' #make the PO4 to Sta

keep <- as.factor(c("CCB", "CCV", "Dup", "Sta")) ## creat the unknowns.n
phos.raw$Code2 <- as.character(ifelse(phos.raw$Code1 %in% keep, phos.raw$Code1, "Unk"))

phos.raw$Code1<- as.character(phos.raw$Code1)## merge together
phos.raw$sample.type <- ifelse(phos.raw$Code2=="Unk", "Unk", phos.raw$Code1)

names(phos.raw) ## remove the columns we dont need
phos <- phos.raw[ , -which(names(phos.raw) %in% c("Code1","Code2"))]
################################################################################
## now we are ready to do the drift correction, I will first do the NH4 then OP
################################################################################
################ NH4 correction of drifts ######################################

## standard curves
stands.nh4 <- subset(amon, sample.type=="Sta")
stands.nh4$exp_conc <- c(0, 0.011, 0.022, 0.05, 0.1, 0.2, 0.5, 1,2, 0)
stands.nh4$std.type <- c(rep("low.std", 7), rep("clean.std", 2), "low.sd")

names(stands.nh4)

## all standards
(all.std.range.n <- rang- range(stands.nh4$exp_conc)) ## the range of concentrations we will use
all.min.abs.n <- min(stands.nh4$Absorbance)
all.max.abs.n <- max(stands.nh4$Absorbance)
all.std.curve.n <- lm(exp_conc~Absorbance, data=stands.nh4)
all.std.curve.intercept.n <- coef(all.std.curve.n)[1]
all.std.curve.slope.n <- coef(all.std.curve.n)[2]
summary(all.std.curve.n)$r.squared 

ggplot(data=stands.nh4, aes(Absorbance, exp_conc)) + geom_point()+ 
  theme_light() + geom_smooth(method="lm", lty=3) + 
  ggtitle("All standards")+ ylab("Expected Concentration")

## cleaned up standards ######################################################
(clean.std.range.n <- range(stands.nh4[stands.nh4$exp_conc>0,]$exp_conc)) ## the range of concentration we will use
clean.min.abs.n <- min(stands.nh4[stands.nh4$exp_conc>0,]$Absorbance)
clean.max.abs.n <- max(stands.nh4[stands.nh4$exp_conc>0,]$Absorbance)
clean.std.curve.n <- lm(exp_conc~Absorbance, data=stands.nh4[stands.nh4$exp_conc>0,])
clean.std.curve.intercept.n <- coef(clean.std.curve.n)[1]
clean.std.curve.slope.n <- coef(clean.std.curve.n)[2]
summary(clean.std.curve.n)$r.squared 

ggplot(data=stands.nh4[stands.nh4$exp_conc>0,], aes(Absorbance, exp_conc)) + geom_point()+ 
  theme_light() + geom_smooth(method="lm", lty=3) + 
  ggtitle("Cleaned up standards")+ ylab("Expected Concentration")

######## standards in the lower range ########################################
low.stands.nh4 <-stands.nh4[stands.nh4$exp_conc >0 & stands.nh4$exp_conc <2, ]
(low.std.range.n <- range(low.stands.nh4$exp_conc)) ## the range of standard we will use
low.min.abs.n <- min(low.stands.nh4$Absorbance)
low.max.abs.n <- max(low.stands.nh4$Absorbance)
low.std.curve.n <- lm(exp_conc~Absorbance, data=low.stands.nh4)
low.std.curve.intercept.n <- coef(low.std.curve.n)[1]
low.std.curve.slope.n <- coef(low.std.curve.n)[2]
summary(low.std.curve.n)$r.squared 

ggplot(data=low.stands.nh4, aes(Absorbance, exp_conc)) + geom_point()+ 
  theme_light() + geom_smooth(method="lm", lty=3) + 
  ggtitle("Low standards")+ ylab("Expected Concentration")

###############################################################################
## CALIBRATION STANDARDS ######################################################
names(stands.nh4)
stands.nh4$calc_conc <- ifelse(stands.nh4$Absorbance < low.max.abs.n,
                               low.std.curve.slope.n * stands.nh4$Absorbance + low.std.curve.intercept.n,
                               clean.std.curve.slope.n * stands.nh4$Absorbance + clean.std.curve.intercept.n) * stands.nh4$Auto.Dil.Factor

stands.nh4$calc_per_error <- ((stands.nh4$calc_conc - stands.nh4$exp_conc) /  stands.nh4$exp_conc) *100

##############################################################################
### CHECK STANDARDS ##########################################################
stands.ccv.n<- subset(amon, sample.type=="CCV")
names(stands.ccv.n)
stands.ccv.n$exp_conc <- rep(0.5, nrow(stands.ccv.n))
stands.ccv.n$std.type <- c(rep("low.std", nrow(stands.ccv.n)))

stands.ccv.n$calc_conc <- ifelse(stands.ccv.n$Absorbance < low.max.abs.n,
                            low.std.curve.slope.n * stands.ccv.n$Absorbance + low.std.curve.intercept.n,
                            clean.std.curve.slope.n * stands.ccv.n$Absorbance + clean.std.curve.intercept.n)* stands.ccv.n$Auto.Dil.Factor

stands.ccv.n$calc_per_error <- ((stands.ccv.n$calc_conc - stands.ccv.n$exp_conc) /  stands.ccv.n$exp_conc) *100

### I need to pull out the times of the measurements and convert it to numeric to be able to fit model
#### be careful here and check how your date and time look
stands.ccv.n$Date.and.Time ## you need to tell it where to start and stop cutting
#we want to cut out just the time, which starts at space 12 and end at space 20 in the string
stands.ccv.n <- stands.ccv.n %>% mutate(time = paste(substr(Date.and.Time,12,20), sep=""))
stands.ccv.n$Time <- chron(times= as.character(stands.ccv.n$time)) ## this makes r see the time as a real time and not a character
stands.ccv.n$seconds <- as.numeric(hms(stands.ccv.n$Time))## then to make the regression i convert to seconds

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

unknowns.n <- amon[amon$sample.type == "Unk" | amon$sample.type == "Dup", ]
## us the same notation with the time as above
unknowns.n<- unknowns.n %>% mutate(time = paste(substr(Date.and.Time,12,20), sep=""))
unknowns.n$Time <- chron(times= as.character(unknowns.n$time)) 
unknowns.n$seconds <- as.numeric(hms(unknowns.n$Time))


unknowns.n$calc_conc <- ifelse(unknowns.n$Absorbance < low.max.abs.n,
                               low.std.curve.slope.n * unknowns.n$Absorbance + low.std.curve.intercept.n,
                               clean.std.curve.slope.n * unknowns.n$Absorbance + clean.std.curve.intercept.n) * unknowns.n$Auto.Dil.Factor


ccv.drift_curve.n <- lm(meas_minus_exp~seconds, data=stands.ccv.n)  ## this is a repeat of above
coef(ccv.drift_curve.n)
ccv.drift.curve.intercept.n <- coef(ccv.drift_curve.n)[1]
cc.drift.curve.slope.n <- coef(ccv.drift_curve.n)[2]


unknowns.n$total.drift <- cc.drift.curve.slope.n * unknowns.n$seconds + ccv.drift.curve.intercept.n 

unknowns.n$actual.conc <- ifelse(unknowns.n$calc_conc > 0.5 * stands.nh4$exp_conc[2],
                               unknowns.n$calc_conc + unknowns.n$total.drift, 0)


### Now I want to make some indentifyer columns and get the data ready for export and
## make some graphs to check the data

## first I will split out the duplicates and we will merge later
## this is because the id variables are in different columns
dups.n <- subset(unknowns.n, sample.type == "Dup")
other.n <- subset(unknowns.n, sample.type == "Unk")

## using library unglue we will break about the sample.id
?unglue_vec ## you can set up the pattern of your data and take out any part
other.n$Sample.ID ## there are three sections {x=site name}_{y=date}_{z=rep} seperatued by _
other.n$site <- unglue_vec(other.n$Sample.ID,"{x}_{y}_{z}", var = "x")
other.n$sample_date <- unglue_vec(other.n$Sample.ID,"{x}_{y}_{z}", var = "y")
other.n$rep <- unglue_vec(other.n$Sample.ID,"{x}_{y}_{z}", var = "z")

head(dups.n)## here the id is in the sample details
dups.n$site <- unglue_vec(dups.n$Sample.Details,"{x}_{y}_{z}", var = "x")
dups.n$sample_date <- unglue_vec(dups.n$Sample.Details,"{x}_{y}_{z}", var = "y")
dups.n$rep <- unglue_vec(dups.n$Sample.Details,"{x}_{y}_{z}", var = "z")

final.nh4 <- rbind(other.n, dups.n)
head(final.nh4)
final.nh4$site <- as.factor(final.nh4$site)
levels(final.nh4$site)
final.nh4$site <- factor(final.nh4$site, levels = c( "HAT",   "HOGDN", "HOGUP", "POS"  , "SWBUP" ,  "SWB" ,"TUM" )) ## put in order or urbanization

 ggplot(data=final.nh4, aes(x=site, y=actual.conc, fill=site)) +
  geom_boxplot() + theme_bw() + geom_jitter(width=0.2, shape=4) + 
  ylab("NH4") + facet_wrap(~sample_date) # take out this last part if you only have one date
 
 ## to see how the duplicates line up with their replicate
 selected <- levels(as.factor(dups.n$Sample.Details)); selected
 levels(as.factor(final.nh4$Sample.ID))
 red_nh4 <- final.nh4[final.nh4$Sample.ID %in% selected,]
 red_nh42 <- rbind(dups.n, red_nh4)
 ggplot(data=red_nh42, aes(x=site, y=actual.conc, color=sample.type)) +
   geom_jitter(data=red_nh42, aes(x=site, y=actual.conc),width=0.3, size=3) + theme_bw() +  
   ylab("NH4 (mg/L)") + facet_wrap(~sample_date)
 
 ## finally export data in a standard way "nutrient_YYYYMMDD"
final.nh4$sample_date
 ## if there are two or more dates you can subset by date and export individually
 # example  final.nh4.a <- subset(final.nh4, sample_date == "01102020")
 write.csv(final.nh4, "nh4_20200110.csv")

 ################################################################################
 ###################################################################################################
 ################################################################################
 ## now we are ready transition to OP
 ################################################################################
 ################ OP correction of drifts ######################################
 ###############################################################################
 head(phos)
 
 ## standard curves
 stands.op <- subset(phos, sample.type=="Sta")
 stands.op$exp_conc <- c(0, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25, 0.5,1, 0)
 stands.op$std.type <- c(rep("low.std", 7), rep("clean.std", 2), "low.sd")
 
 names(stands.op)
 
 ## all standards
 (all.std.range.p <- range(stands.op$exp_conc)) ## the range of concentraion of standard we will use
 all.min.abs.p <- min(stands.op$Absorbance)
 all.max.abs.p <- max(stands.op$Absorbance)
 all.std.curve.p <- lm(exp_conc~Absorbance, data=stands.op)
 all.std.curve.intercept.p <- coef(all.std.curve.p)[1]
 all.std.curve.slope.p <- coef(all.std.curve.p)[2]
 summary(all.std.curve.p)$r.squared 
 coef(all.std.curve.p)
 confint(all.std.curve.p, level = 0.95) 
 
 ggplot(data=stands.op, aes(Absorbance, exp_conc)) + geom_point()+ 
   theme_light() + geom_smooth(method="lm", lty=3) + 
   ggtitle("All standards")+ ylab("Expected Concentration")
 
 ## cleaned up standards ######################################################
 (clean.std.range.p <- range(stands.op[stands.op$exp_conc>0,]$exp_conc)) ## the range of concentration we will use
 clean.min.abs.p <- min(stands.op[stands.op$exp_conc>0,]$Absorbance)
 clean.max.abs.p <- max(stands.op[stands.op$exp_conc>0,]$Absorbance)
 clean.std.curve.p <- lm(exp_conc~Absorbance, data=stands.op[stands.op$exp_conc>0,])
 clean.std.curve.intercept.p <- coef(clean.std.curve.p)[1]
 clean.std.curve.slope.p <- coef(clean.std.curve.p)[2]
 summary(clean.std.curve.p)$r.squared 
 coef(clean.std.curve.p)
 confint(clean.std.curve.p, level = 0.95) 
 
 ggplot(data=stands.op[stands.op$exp_conc>0,], aes(Absorbance, exp_conc)) + geom_point()+ 
   theme_light() + geom_smooth(method="lm", lty=3) + 
   ggtitle("Cleaned up standards")+ ylab("Expected Concentration")
 
 ######## standards in the lower range ########################################
 low.stands.op <-stands.op[stands.op$exp_conc >0 & stands.op$exp_conc <2, ]
 (low.std.range.p <- range(low.stands.op$exp_conc)) ## the range of standard we will use
 low.min.abs.p <- min(low.stands.op$Absorbance)
 low.max.abs.p <- max(low.stands.op$Absorbance)
 low.std.curve.p <- lm(exp_conc~Absorbance, data=low.stands.op)
 low.std.curve.intercept.p <- coef(low.std.curve.p)[1]
 low.std.curve.slope.p <- coef(low.std.curve.p)[2]
 summary(low.std.curve.p)$r.squared 
 coef(low.std.curve.p)
 confint(low.std.curve.p, level = 0.95) 
 
 ggplot(data=low.stands.op, aes(Absorbance, exp_conc)) + geom_point()+ 
   theme_light() + geom_smooth(method="lm", lty=3) + 
   ggtitle("Low standards")+ ylab("Expected Concentration")
 
 ###############################################################################
 ## CALIBRATION STANDARDS ######################################################
 names(stands.op)
 stands.op$calc_conc <- ifelse(stands.op$Absorbance < low.max.abs.p,
                                low.std.curve.slope.p * stands.op$Absorbance + low.std.curve.intercept.p,
                                clean.std.curve.slope.p * stands.op$Absorbance + clean.std.curve.intercept.p) * stands.op$Auto.Dil.Factor
 
 stands.op$calc_per_error <- ((stands.op$calc_conc - stands.op$exp_conc) /  stands.op$exp_conc) *100
 
 ##############################################################################
 ### CHECK STANDARDS ##########################################################
 stands.ccv.p<- subset(phos, sample.type=="CCV")
 names(stands.ccv.p)
 stands.ccv.p$exp_conc <- rep(0.1, nrow(stands.ccv.p))
 stands.ccv.p$std.type <- c(rep("low.std", nrow(stands.ccv.p)))
 
 stands.ccv.p$calc_conc <- ifelse(stands.ccv.p$Absorbance < low.max.abs.p,
                                  low.std.curve.slope.p * stands.ccv.p$Absorbance + low.std.curve.intercept.p,
                                  clean.std.curve.slope.p * stands.ccv.p$Absorbance + clean.std.curve.intercept.p)* stands.ccv.p$Auto.Dil.Factor
 
 stands.ccv.p$calc_per_error <- ((stands.ccv.p$calc_conc - stands.ccv.p$exp_conc) /  stands.ccv.p$exp_conc) *100
 
 ### I need to pull out the times of the meausrements and convert it to numeric to be able to fit model
 #### be careful here aand check how your date and time look
 stands.ccv.p$Date.and.Time ## you need to tell it where to start and stop cutting
 #we want to cut out just the time, which starts at space 12 and end at space 20 in the string
 stands.ccv.p <- stands.ccv.p %>% mutate(time = paste(substr(Date.and.Time,12,20), sep=""))
 stands.ccv.p$Time <- chron(times= as.character(stands.ccv.p$time)) ## this makes r see the time as a real time and not a charachter
 stands.ccv.p$seconds <- as.numeric(hms(stands.ccv.p$Time))## then to make the regression i convert to seconds

 ccv.curve.p <- lm(calc_conc~seconds, data=stands.ccv.p)
 coef(ccv.curve.p)
 ccv.curve.intercept.p <- coef(ccv.curve.p)[1]
 cc.curve.slope.p <- coef(ccv.curve.p)[2]
 
 ggplot(data=stands.ccv.p, aes(x=seconds, y=calc_conc)) + geom_point()+ 
   theme_light() + geom_smooth(method="lm", lty=3) + 
   ggtitle("Drift Check")+ ylab("Calculated Concentration") + xlab("run time")
 
 stands.ccv.p$expected.ccv <- cc.curve.slope.p * stands.ccv.p$seconds + ccv.curve.intercept.p
 
 stands.ccv.p$meas_minus_exp <- stands.ccv.p$calc_conc[1] - stands.ccv.p$expected.ccv  
 
 stands.ccv.p$drift_corrected_ccv <- stands.ccv.p$calc_conc + stands.ccv.p$meas_minus_exp
 
 ggplot(data=stands.ccv.p, aes(x=seconds, y=drift_corrected_ccv)) + geom_point()+ 
   theme_light() + geom_smooth(method="lm", lty=3) + 
   ggtitle("Drift Check")+ ylab("Drift Corrected Concentration") + xlab("run time")
 
 ###########################################################################################
 ##### now to correct our unknown data and the Duplicates
 
 unknowns.p <- phos[phos$sample.type == "Unk" | phos$sample.type == "Dup", ]
 ## us the same notation with the time as above
 unknowns.p<- unknowns.p %>% mutate(time = paste(substr(Date.and.Time,12,20), sep=""))
 unknowns.p$Time <- chron(times= as.character(unknowns.p$time)) 
 unknowns.p$seconds <- as.numeric(hms(unknowns.p$Time))
 
 
 unknowns.p$calc_conc <- ifelse(unknowns.p$Absorbance < low.max.abs.p,
                                low.std.curve.slope.p * unknowns.p$Absorbance + low.std.curve.intercept.p,
                                clean.std.curve.slope.p * unknowns.p$Absorbance + clean.std.curve.intercept.p)*  unknowns.p$Auto.Dil.Factor
 
 
 ccv.drift_curve.p <- lm(meas_minus_exp~seconds, data=stands.ccv.p)
 coef(ccv.drift_curve.p)
 ccv.drift.curve.intercept.p <- coef(ccv.drift_curve.p)[1]
 cc.drift.curve.slope.p <- coef(ccv.drift_curve.p)[2]
 
 
 unknowns.p$total.drift <- cc.drift.curve.slope.p * unknowns.p$seconds + ccv.drift.curve.intercept.p 
 
 unknowns.p$actual.conc <- ifelse(unknowns.p$calc_conc > 0.5 * stands.op$exp_conc[2],
                                  unknowns.p$calc_conc + unknowns.p$total.drift, 0)
 
 
 ### Now I want to make some indentifyer columns and get the data ready for export and
 ## make some graphs to check the data
 
 ## first I will split out the duplicates and we will merge later
 ## this is because the id variables are in different columns
 dups.p <- subset(unknowns.p, sample.type == "Dup")
 other.p <- subset(unknowns.p, sample.type == "Unk")
 
 ## using library unglue we will breack about the sample.id
 ?unglue_vec ## you can set up the pattern of your data and take out any part
 other.p$Sample.ID ## there are three sections {x=site name}_{y=date}_{z=rep} seperatued by _
 other.p$site <- unglue_vec(other.p$Sample.ID,"{x}_{y}_{z}", var = "x")
 other.p$sample_date <- unglue_vec(other.p$Sample.ID,"{x}_{y}_{z}", var = "y")
 other.p$rep <- unglue_vec(other.p$Sample.ID,"{x}_{y}_{z}", var = "z")
 
 head(dups.p)## here the id is in the sample details
 dups.p$site <- unglue_vec(dups.p$Sample.Details,"{x}_{y}_{z}", var = "x")
 dups.p$sample_date <- unglue_vec(dups.p$Sample.Details,"{x}_{y}_{z}", var = "y")
 dups.p$rep <- unglue_vec(dups.p$Sample.Details,"{x}_{y}_{z}", var = "z")
 
 final.op <- rbind(other.p, dups.p)
 head(final.op)
 final.op$site <- as.factor(final.op$site)
 levels(final.op$site)
 final.op$site <- factor(final.op$site, levels = c( "HAT",   "HOGDN", "HOGUP", "POS"  , "SWBUP" ,  "SWB" ,"TUM" )) ## put in order or urbanization
 
 
 ggplot(data=final.op, aes(x=site, y=actual.conc, fill=site)) +
   geom_boxplot() + theme_bw() + geom_jitter(width=0.2, shape=4) + 
   ylab("OP units)") + facet_wrap(~sample_date) ## takes this last part our if you only have one date
 
 ## to see how the dulplicates line up with their replicate
 selected <- levels(as.factor(dups.p$Sample.Details)); selected
 levels(as.factor(final.op$Sample.ID))
 red_op <- final.op[final.op$Sample.ID %in% selected,]
 red_op2 <- rbind(dups.p, red_op)
 ggplot(data=red_op2, aes(x=site, y=actual.conc, color=sample.type)) +
   geom_jitter(data=red_op2, aes(x=site, y=actual.conc),width=0.3, size=3) + theme_bw() +  
   ylab("op units)") +  facet_wrap(~sample_date) ## takes this last part our if you only have one date
 
 ## finally export data in a standard way "nutrient_YYYYMMDD"
 final.op$sample_date
 ## if there are two or more dates you can subset by date and export individually
 # example  final.op.a <- subset(final.op, sample_date == "01102020")
 write.csv(final.op, "OP_20200110.csv")
 
 
 

