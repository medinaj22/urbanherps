#this is the predictive analysis of species richness in parks in individual cities
#each city was sequentially left out of the model fitting, then the model was
#used to predict species richness in that city

library(ggplot2)
library(lme4)
library(rstan)
library(brms)
rstan_options(auto_write = TRUE)

dat<-ALLCITIES8


#data processing
dat$City<-as.factor(dat$City)
dat$LOC_FINAL<-as.factor(dat$LOC_FINAL)
dat$SHAPE<-dat$perimeter^2/dat$area
dat$forarea<-dat$area*dat$forestcover/100
dat$visitors<-(dat$Junevisits+dat$Julyvisits)/2
dat$ROADDENS<-dat$road_length/dat$area
dat$HERP_SP <- dat$SAL_SP + dat$FROG_SP + dat$SNAKE_SP + dat$LIZARD_SP + dat$TURTLE_SP
dat$AMPH_SP <- dat$SAL_SP + dat$FROG_SP
dat$REP_SP <- dat$SNAKE_SP + dat$LIZARD_SP + dat$TURTLE_SP
dat$logherps<-log(dat$HERP_SP)
summary(dat)

#need to re-scale variables
dat$forest_prop<-dat$forestcover/100
dat$imperv_prop<-dat$imperv/100
dat$fp2<-dat$forest_prop^2
dat$normCONN10<-(dat$CONN1000-mean(dat$CONN1000, na.rm=TRUE))/(sd(dat$CONN1000, na.rm=TRUE))
dat$normCONN5<-(dat$CONN500-mean(dat$CONN500, na.rm=TRUE))/(sd(dat$CONN500, na.rm=TRUE))
dat$normCONN1<-(dat$CONN100-mean(dat$CONN100, na.rm=TRUE))/(sd(dat$CONN100, na.rm=TRUE))
dat$normSHAPE<-(dat$SHAPE-mean(dat$SHAPE, na.rm=TRUE))/(sd(dat$SHAPE, na.rm=TRUE))
dat$normVIS<-(dat$visitors-mean(dat$visitors, na.rm=TRUE))/(sd(dat$visitors, na.rm=TRUE))
dat$logarea<-log(dat$area)
dat$normROAD<-(dat$ROADDENS-mean(dat$ROADDENS, na.rm=TRUE))/(sd(dat$ROADDENS, na.rm=TRUE))


#prior for coefficients
prior1<-set_prior("normal(0,10)", class="b")


#LOOK AT PREDICTIONS BASED ON LEAVE ONE OUT BY CITY
#FIRST CREATE SUBSETS
#test data for each city
DCdata<-subset(dat, City=="DC")
BALTdata<-subset(dat, City=="BALTIMORE")
INDYdata<-subset(dat, City=="Indianapolis")
PHILdata<-subset(dat, City=="Philadelphia")
BRONXdata<-subset(dat, City=="Bronx")
MANHATdata<-subset(dat, City=="Manhattan")
QUEENSdata<-subset(dat, City=="Queens County, NY")
BROOKdata<-subset(dat, City=="Kings County, NY")
SIdata<-subset(dat, City=="Richmond County, NY")
PITTdata<-subset(dat, City=="Pittsburgh ")

#training data with each city left out
looDC<-subset(dat, City!="DC")
looBALT<-subset(dat, City!="BALTIMORE")
looINDY<-subset(dat, City!="Indianapolis")
looPHIL<-subset(dat, City!="Philadelphia")
looBRONX<-subset(dat, City!="Bronx")
looMANHAT<-subset(dat, City!="Manhattan")
looQUEENS<-subset(dat, City!="Queens County, NY")
looBROOK<-subset(dat, City!="Kings County, NY")
looSI<-subset(dat, City!="Richmond County, NY")
looPITT<-subset(dat, City!="Pittsburgh ")

#PREDICTIONS FOR BALTIMORE
BRMnoBALT<-brm(HERP_SP ~ logarea + wetland_types + forest_prop + normCONN5 + (1|City), 
               prior=prior1, family = "negbinomial", data=looBALT)
summary(BRMnoBALT)

fitted_values <- fitted(BRMnoBALT, newdata = BALTdata, re_formula = NA)
head(fitted_values)

BALTpreds<-data.frame(BALTdata$LOC_FINAL, BALTdata$HERP_SP, fitted_values)
BALTpreds

BALTpreds$pred_dev<-abs(BALTpreds$BALTdata.HERP_SP-BALTpreds$Estimate)
mean(BALTpreds$pred_dev)
BALTpreds$pred_dev2<-(BALTpreds$BALTdata.HERP_SP-BALTpreds$Estimate)^2
mean(BALTpreds$pred_dev2)

ggplot(BALTpreds, aes(Estimate, BALTdata.HERP_SP))+
  geom_point(size=5)+
  geom_abline(slope=1, intercept=0)+
  ggtitle("Baltimore") +
  labs(x="Predicted Species", y="Observed Species")+
  xlim(0,8)+
  ylim(0,8)+
  theme(plot.title = element_text(hjust = 0.5, size = 24))+
  theme(axis.title=element_text(size=20))+
  theme(axis.text=element_text(size=18))


#leave out Bronx
BRMnoBRONX<-brm(HERP_SP ~ logarea + wetland_types + forest_prop + normCONN5 + (1|City), 
                prior=prior1, family = "negbinomial", data=looBRONX)
summary(BRMnoBRONX)

fitted_values <- fitted(BRMnoBRONX, newdata = BRONXdata, re_formula = NA)
head(fitted_values)

BRONXpreds<-data.frame(BRONXdata$LOC_FINAL, BRONXdata$HERP_SP, fitted_values)
BRONXpreds

BRONXpreds$pred_dev<-abs(BRONXpreds$BRONXdata.HERP_SP-BRONXpreds$Estimate)
mean(BRONXpreds$pred_dev)

BRONXpreds$pred_dev2<-(BRONXpreds$BRONXdata.HERP_SP-BRONXpreds$Estimate)^2
mean(BRONXpreds$pred_dev2)

ggplot(BRONXpreds, aes(Estimate, BRONXdata.HERP_SP))+
  geom_point(size=5)+
  geom_abline(slope=1, intercept=0)+
  ggtitle("Bronx") +
  labs(x="Predicted species", y="Observed species")+
  xlim(0,10)+
  ylim(0,10)+
  theme(plot.title = element_text(hjust = 0.5, size = 24))+
  theme(axis.title=element_text(size=20))+
  theme(axis.text=element_text(size=18))


#leave out Brooklyn
BRMnoBROOK<-brm(HERP_SP ~ logarea + wetland_types + forest_prop + normCONN5 + (1|City), 
                prior=prior1, family = "negbinomial", data=looBROOK)
summary(BRMnoBROOK)

fitted_values <- fitted(BRMnoBROOK, newdata = BROOKdata, re_formula = NA)
head(fitted_values)

BROOKpreds<-data.frame(BROOKdata$LOC_FINAL, BROOKdata$HERP_SP, fitted_values)
BROOKpreds

BROOKpreds$pred_dev<-abs(BROOKpreds$BROOKdata.HERP_SP-BROOKpreds$Estimate)
mean(BROOKpreds$pred_dev)
BROOKpreds$pred_dev2<-(BROOKpreds$BROOKdata.HERP_SP-BROOKpreds$Estimate)^2
mean(BROOKpreds$pred_dev2)

ggplot(BROOKpreds, aes(Estimate, BROOKdata.HERP_SP))+
  geom_point(size=5)+
  geom_abline(slope=1, intercept=0)+
  ggtitle("Brooklyn") +
  labs(x="Predicted species", y="Observed species")+
  xlim(0,12)+
  ylim(0,12)+
  theme(plot.title = element_text(hjust = 0.5, size = 24))+
  theme(axis.title=element_text(size=20))+
  theme(axis.text=element_text(size=18))



#PREDICTIONS FOR DC
BRMnoDC<-brm(HERP_SP ~ logarea + wetland_types + forest_prop + normCONN5 + (1|City), 
             family = "negbinomial", prior=prior1, data=looDC)
summary(BRMnoDC)

fitted_values <- fitted(BRMnoDC, newdata = DCdata, re_formula = NA)
head(fitted_values)

DCpreds<-data.frame(DCdata$LOC_FINAL, DCdata$HERP_SP, fitted_values)
DCpreds
DCpreds$pred_dev<-abs(DCpreds$DCdata.HERP_SP-DCpreds$Estimate)
mean(DCpreds$pred_dev)
DCpreds$pred_dev2<-(DCpreds$DCdata.HERP_SP-DCpreds$Estimate)^2
mean(DCpreds$pred_dev2)


ggplot(DCpreds, aes(Estimate, DCdata.HERP_SP))+
  geom_point(size=5)+
  geom_abline(slope=1, intercept=0)+
  ggtitle("DC") +
  labs(x="Predicted species", y="Observed species")+
  xlim(0,20)+
  ylim(0,20)+
  theme(plot.title = element_text(hjust = 0.5, size = 24))+
  theme(axis.title=element_text(size=20))+
  theme(axis.text=element_text(size=18))


#PREDICTIONS FOR INDIANAPOLIS
BRMnoINDY<-brm(HERP_SP ~ logarea + wetland_types + forest_prop + normCONN5 + (1|City), 
               prior=prior1, family = "negbinomial", data=looINDY)
summary(BRMnoINDY)

fitted_values <- fitted(BRMnoINDY, newdata = INDYdata, re_formula = NA)
head(fitted_values)

INDYpreds<-data.frame(INDYdata$LOC_FINAL, INDYdata$HERP_SP, fitted_values)
INDYpreds

INDYpreds$pred_dev<-abs(INDYpreds$INDYdata.HERP_SP-INDYpreds$Estimate)
mean(INDYpreds$pred_dev)
INDYpreds$pred_dev2<-(INDYpreds$INDYdata.HERP_SP-INDYpreds$Estimate)^2
mean(INDYpreds$pred_dev2)

ggplot(INDYpreds, aes(Estimate, INDYdata.HERP_SP))+
  geom_point(size=5)+
  geom_abline(slope=1, intercept=0)+
  ggtitle("Indianapolis") +
  labs(x="Predicted species", y="Observed species")+
  xlim(0,20)+
  ylim(0,20)+
  theme(plot.title = element_text(hjust = 0.5, size = 24))+
  theme(axis.title=element_text(size=20))+
  theme(axis.text=element_text(size=18))



#leave out Manhattan
BRMnoMANHAT<-brm(HERP_SP ~ logarea + wetland_types + forest_prop + normCONN5 + (1|City), 
                 prior=prior1, family = "negbinomial", data=looMANHAT)
summary(BRMnoMANHAT)

fitted_values <- fitted(BRMnoMANHAT, newdata = MANHATdata, re_formula = NA)
head(fitted_values)

MANHATpreds<-data.frame(MANHATdata$LOC_FINAL, MANHATdata$HERP_SP, fitted_values)
MANHATpreds

MANHATpreds$pred_dev<-abs(MANHATpreds$MANHATdata.HERP_SP-MANHATpreds$Estimate)
mean(MANHATpreds$pred_dev)

MANHATpreds$pred_dev2<-(MANHATpreds$MANHATdata.HERP_SP-MANHATpreds$Estimate)^2
mean(MANHATpreds$pred_dev2)

ggplot(MANHATpreds, aes(Estimate, MANHATdata.HERP_SP))+
  geom_point(size=5)+
  geom_abline(slope=1, intercept=0)+
  ggtitle("Manhattan") +
  labs(x="Predicted species", y="Observed species")+
  xlim(0,10)+
  ylim(0,10)+
  theme(plot.title = element_text(hjust = 0.5, size = 24))+
  theme(axis.title=element_text(size=20))+
  theme(axis.text=element_text(size=18))



#PREDICTIONS FOR PHILADELPHIA
BRMnoPHIL<-brm(HERP_SP ~ logarea + wetland_types + forest_prop + normCONN5 + (1|City), 
               prior=prior1, family = "negbinomial", data=looPHIL)
summary(BRMnoPHIL)

fitted_values <- fitted(BRMnoPHIL, newdata = PHILdata, re_formula = NA)
head(fitted_values)

PHILpreds<-data.frame(PHILdata$LOC_FINAL, PHILdata$HERP_SP, fitted_values)
PHILpreds

PHILpreds$pred_dev<-abs(PHILpreds$PHILdata.HERP_SP-PHILpreds$Estimate)
mean(PHILpreds$pred_dev)
PHILpreds$pred_dev2<-(PHILpreds$PHILdata.HERP_SP-PHILpreds$Estimate)^2
mean(PHILpreds$pred_dev2)

ggplot(PHILpreds, aes(Estimate, PHILdata.HERP_SP))+
  geom_point(size=5)+
  geom_abline(slope=1, intercept=0)+
  ggtitle("Philadelphia") +
  labs(x="Predicted species", y="Observed species")+
  xlim(0,15)+
  ylim(0,15)+
  theme(plot.title = element_text(hjust = 0.5, size = 24))+
  theme(axis.title=element_text(size=20))+
  theme(axis.text=element_text(size=18))



#PREDICTIONS FOR PITTSBURGH
BRMnoPITT<-brm(HERP_SP ~ logarea + wetland_types + forest_prop + normCONN5 + (1|City), 
               prior=prior1, family = "negbinomial", data=looPITT)
summary(BRMnoPITT)

fitted_values <- fitted(BRMnoPITT, newdata = PITTdata, re_formula = NA)
head(fitted_values)

PITTpreds<-data.frame(PITTdata$LOC_FINAL, PITTdata$HERP_SP, fitted_values)
PITTpreds

PITTpreds$pred_dev<-abs(PITTpreds$PITTdata.HERP_SP-PITTpreds$Estimate)
mean(PITTpreds$pred_dev)
PITTpreds$pred_dev2<-(PITTpreds$PITTdata.HERP_SP-PITTpreds$Estimate)^2
mean(PITTpreds$pred_dev2)

ggplot(PITTpreds, aes(Estimate, PITTdata.HERP_SP))+
  geom_point(size=5)+
  geom_abline(slope=1, intercept=0)+
  ggtitle("Pittsburgh") +
  labs(x="Predicted species", y="Observed species")+
  xlim(0,15)+
  ylim(0,15)+
  theme(plot.title = element_text(hjust = 0.5, size = 24))+
  theme(axis.title=element_text(size=20))+
  theme(axis.text=element_text(size=18))



#leave out Queens
BRMnoQUEENS<-brm(HERP_SP ~ logarea + wetland_types + forest_prop + normCONN5 + (1|City), 
                 prior=prior1, family = "negbinomial", data=looQUEENS)
summary(BRMnoQUEENS)



#leave out staten island
BRMnoSI<-brm(HERP_SP ~ logarea + wetland_types + forest_prop + normCONN5 + (1|City), 
             prior=prior1, family = "negbinomial", data=looSI)
summary(BRMnoSI)

fitted_values <- fitted(BRMnoSI, newdata = SIdata, re_formula = NA)
head(fitted_values)

SIpreds<-data.frame(SIdata$LOC_FINAL, SIdata$HERP_SP, fitted_values)
SIpreds

SIpreds$pred_dev<-abs(SIpreds$SIdata.HERP_SP-SIpreds$Estimate)
mean(SIpreds$pred_dev)

SIpreds$pred_dev2<-(SIpreds$SIdata.HERP_SP-SIpreds$Estimate)^2
mean(SIpreds$pred_dev2)
          
ggplot(SIpreds, aes(Estimate, SIdata.HERP_SP))+
  geom_point(size=5)+
  geom_abline(slope=1, intercept=0)+
  ggtitle("Staten Island") +
  labs(x="Predicted species", y="Observed species")+
  xlim(0,15)+
  ylim(0,15)+
  theme(plot.title = element_text(hjust = 0.5, size = 24))+
  theme(axis.title=element_text(size=20))+
  theme(axis.text=element_text(size=18))



