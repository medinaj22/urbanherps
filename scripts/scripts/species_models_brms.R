#models for presence/absence of individual species

dat<-species_compilation
library(ggplot2)
library(lme4)
library(rstan)
library(brms)


#data processing
dat$City<-as.factor(dat$City)
dat$ROADDENS<-dat$road_length/dat$area
dat$SHAPE<-dat$perimeter^2/dat$area
summary(dat)

#rescaling variables
dat$forest_prop<-dat$forestcover/100
dat$normCONN10<-(dat$CONN1000-mean(dat$CONN1000, na.rm=TRUE))/(sd(dat$CONN1000, na.rm=TRUE))
dat$normCONN5<-(dat$CONN500-mean(dat$CONN500, na.rm=TRUE))/(sd(dat$CONN500, na.rm=TRUE))
dat$normCONN1<-(dat$CONN100-mean(dat$CONN100, na.rm=TRUE))/(sd(dat$CONN100, na.rm=TRUE))
dat$normSHAPE<-(dat$SHAPE-mean(dat$SHAPE, na.rm=TRUE))/(sd(dat$SHAPE, na.rm=TRUE))
dat$logarea<-log(dat$area)
dat$normROAD<-(dat$ROADDENS-mean(dat$ROADDENS, na.rm=TRUE))/(sd(dat$ROADDENS, na.rm=TRUE))

#prior for coefficients
prior1<-set_prior("normal(0,10)", class="b")

#Redback salamanders
dat$RB10<-ifelse(dat$RBS_OBS>0,1,0) #creates new 1/0 column for RBS


RBmod<-glmer(RB10~logarea + forest_prop + wetland_types + normCONN5 + (1|City), 
           family="binomial", data=dat)
summary(RBmod)

RBocc<-mean(dat$RB10)  #proportion of parks with RBS

RB_brms<-brm(RB10~logarea + forest_prop + wetland_types + normCONN5 + (1|City), 
             prior=prior1, family="bernoulli", data=dat)
summary(RB_brms)
mcmc_plot(RB_brms, pars=c("b_logarea", "b_wetland_types", "b_forest_prop", "b_normCONN5"))

#two-lined salamanders
dat$twolined10<-ifelse(dat$LINED_OBS>0,1,0) #creates new 1/0 column for 2LS
TLmod<-glmer(twolined10~logarea + forest_prop + wetland_types + normCONN5 + (1|City), 
             family="binomial", data=dat)
summary(TLmod)

TLocc<-mean(dat$twolined10)

TL_brms<-brm(twolined10~logarea + forest_prop + wetland_types + normCONN5 + (1|City), 
             prior=prior1, family="bernoulli", data=dat)
summary(TL_brms)
mcmc_plot(TL_brms, pars=c("b_logarea", "b_wetland_types", "b_forest_prop", "b_normCONN5"))

hypothesis(TL_brms, "logarea>0")

#painted turtles
dat$PAINT10<-ifelse(dat$PAINT_OBS>0,1,0) #creates new 1/0 column for PAINT
PAINTmod<-glmer(PAINT10~logarea + forest_prop + wetland_types + normCONN5 + (1|City), 
             family="binomial", data=dat)
summary(PAINTmod)

PAINT_brms<-brm(PAINT10~logarea + forest_prop + wetland_types + normCONN5 + (1|City), 
             prior=prior1, family="bernoulli", data=dat)
summary(PAINT_brms)
mcmc_plot(PAINT_brms, pars=c("b_logarea", "b_wetland_types", "b_forest_prop", "b_normCONN5"))

PAINTocc<-mean(dat$PAINT10)

#snapping turtles

dat$SNAP10<-ifelse(dat$SNAP_OBS>0,1,0) #creates new 1/0 column for SNAP
SNAPmod<-glmer(SNAP10~logarea + forest_prop + wetland_types + normCONN5 + (1|City), 
               family="binomial", data=dat)
summary(SNAPmod)

SNAPocc<-mean(dat$SNAP10)

SNAP_brms<-brm(SNAP10~logarea + forest_prop + wetland_types + normCONN5 + (1|City), 
                prior=prior1, family="bernoulli", data=dat)
summary(SNAP_brms)
mcmc_plot(SNAP_brms, pars=c("b_logarea", "b_wetland_types", "b_forest_prop", "b_normCONN5"))
hypothesis(SNAP_brms, "wetland_types>0")
hypothesis(SNAP_brms, "normCONN5>0")


#american toads

dat$AMTO10<-ifelse(dat$AMTO_OBS>0,1,0) #creates new 1/0 column for AMTO
AMTOmod<-glmer(AMTO10~logarea + forest_prop + wetland_types + normCONN5 + (1|City), 
               family="binomial", data=dat)
summary(AMTOmod)

AMTOocc<-mean(dat$AMTO10)

AMTO_brms<-brm(AMTO10~logarea + forest_prop + wetland_types + normCONN5 + (1|City), 
               prior=prior1, family="bernoulli", data=dat)
summary(AMTO_brms)
mcmc_plot(AMTO_brms, pars=c("b_logarea", "b_wetland_types", "b_forest_prop", "b_normCONN5"))

hypothesis(AMTO_brms, "normCONN5>0")

#fowler's toads

dat$FOWL10<-ifelse(dat$FOWL_OBS>0,1,0) #creates new 1/0 column for FOWL
FOWLmod<-glmer(FOWL10~logarea + forest_prop + wetland_types + normCONN5 + (1|City), 
               family="binomial", data=dat)
summary(FOWLmod)

FOWLocc<-mean(dat$FOWL10, na.rm=TRUE)

FOWL_brms<-brm(FOWL10~logarea + forest_prop + wetland_types + normCONN5 + (1|City), 
               prior=prior1, family="bernoulli", data=dat)
summary(FOWL_brms)
mcmc_plot(FOWL_brms, pars=c("b_logarea", "b_wetland_types", "b_forest_prop", "b_normCONN5"))

#spring peeper

dat$PEEP10<-ifelse(dat$PEEP_OBS>0,1,0) #creates new 1/0 column for PEEP
PEEPmod<-glmer(PEEP10~logarea + forest_prop + wetland_types + normCONN5 + (1|City), 
               family="binomial", data=dat)
summary(PEEPmod)

PEEPocc<-mean(dat$PEEP10)

PEEP_brms<-brm(PEEP10~logarea + forest_prop + wetland_types + normCONN5 + (1|City), 
               prior=prior1, family="bernoulli", data=dat)
summary(PEEP_brms)
mcmc_plot(PEEP_brms, pars=c("b_logarea", "b_wetland_types", "b_forest_prop", "b_normCONN5"))

#gray treefrog (combines Hyla chrysoscelis and H. versicolor)

dat$GRAYTR10<-ifelse(dat$GRAYTR_OBS>0,1,0) #creates new 1/0 column for GRAYTR
GRAYTRmod<-glmer(GRAYTR10~logarea + forest_prop + wetland_types + normCONN5 + (1|City), 
               family="binomial", data=dat)
summary(GRAYTRmod)

GRAYTRocc<-mean(dat$GRAYTR10)

GRAYTR_brms<-brm(GRAYTR10~logarea + forest_prop + wetland_types + normCONN5 + (1|City), 
               prior=prior1, family="bernoulli", data=dat)
summary(GRAYTR_brms)
mcmc_plot(GRAYTR_brms, pars=c("b_logarea", "b_wetland_types", "b_forest_prop", "b_normCONN5"))

hypothesis(GRAYTR_brms, "logarea>0")
hypothesis(GRAYTR_brms, "wetland_types>0")

#American bullfrog

dat$BULL10<-ifelse(dat$BULL_OBS>0,1,0) #creates new 1/0 column for BULL
BULLmod<-glmer(BULL10~logarea + forest_prop + wetland_types + normCONN5 + (1|City), 
                 family="binomial", data=dat)
summary(BULLmod)

BULLocc<-mean(dat$BULL10, na.rm=TRUE)

BULL_brms<-brm(BULL10~logarea + forest_prop + wetland_types + normCONN5 + (1|City), 
                 prior=prior1, family="bernoulli", data=dat)
summary(BULL_brms)
mcmc_plot(BULL_brms, pars=c("b_logarea", "b_wetland_types", "b_forest_prop", "b_normCONN5"))
hypothesis(BULL_brms, "logarea>0")

#green frog

dat$GRFR10<-ifelse(dat$GRFR_OBS>0,1,0) #creates new 1/0 column for GRFR
GRFRmod<-glmer(GRFR10~logarea + forest_prop + wetland_types + normCONN5 + (1|City), 
               family="binomial", data=dat)
summary(GRFRmod)

GRFRocc<-mean(dat$GRFR10)

GRFR_brms<-brm(GRFR10~logarea + forest_prop + wetland_types + normCONN5 + (1|City), 
               prior=prior1, family="bernoulli", data=dat)
summary(GRFR_brms)
mcmc_plot(GRFR_brms, pars=c("b_logarea", "b_wetland_types", "b_forest_prop", "b_normCONN5"))

#five-line skink

dat$SKINK10<-ifelse(dat$SKINK_OBS>0,1,0) #creates new 1/0 column for SKINK
SKINKmod<-glmer(SKINK10~logarea + forest_prop + wetland_types + normCONN5 + (1|City), 
               family="binomial", data=dat)
summary(SKINKmod)

SKINKocc<-mean(dat$SKINK10)

SKINK_brms<-brm(SKINK10~logarea + forest_prop + wetland_types + normCONN5 + (1|City), 
                prior=prior1, family="bernoulli", data=dat)
summary(SKINK_brms)
mcmc_plot(SKINK_brms, pars=c("b_logarea", "b_wetland_types", "b_forest_prop", "b_normCONN5"))
hypothesis(SKINK_brms, "normCONN5>0")

#garter snake

dat$GART10<-ifelse(dat$GART_OBS>0,1,0) #creates new 1/0 column for GART
GARTmod<-glmer(GART10~logarea + forest_prop + wetland_types + normCONN5 + (1|City), 
                family="binomial", data=dat)
summary(GARTmod)

GARTocc<-mean(dat$GART10)

GART_brms<-brm(GART10~logarea + forest_prop + wetland_types + normCONN5 + (1|City), 
                prior=prior1, family="bernoulli", data=dat)
summary(GART_brms)
mcmc_plot(GART_brms, pars=c("b_logarea", "b_wetland_types", "b_forest_prop", "b_normCONN5"))

#eastern rat snake

dat$ERAT10<-ifelse(dat$ERAT_OBS>0,1,0) #creates new 1/0 column for ERAT
ERATmod<-glmer(ERAT10~logarea + forest_prop + wetland_types + normCONN5 + (1|City), 
               family="binomial", data=dat)
summary(ERATmod)

ERATocc<-mean(dat$ERAT10)

ERAT_brms<-brm(ERAT10~logarea + forest_prop + wetland_types + normCONN5 + (1|City), 
               prior=prior1, family="bernoulli", data=dat)
summary(ERAT_brms)
mcmc_plot(ERAT_brms, pars=c("b_logarea", "b_wetland_types", "b_forest_prop", "b_normCONN5"))
hypothesis(ERAT_brms, "logarea>0")
hypothesis(ERAT_brms, "wetland_types>0")
hypothesis(ERAT_brms, "forest_prop<0")

#dekay's brown snake

dat$DEKAY10<-ifelse(dat$DEKAY_OBS>0,1,0) #creates new 1/0 column for DEKAY
DEKAYmod<-glmer(DEKAY10~logarea + forest_prop + wetland_types + normCONN5 + (1|City), 
               family="binomial", data=dat)
summary(DEKAYmod)

DEKAYocc<-mean(dat$DEKAY10)

DEKAY_brms<-brm(DEKAY10~logarea + forest_prop + wetland_types + normCONN5 + (1|City), 
               prior=prior1, family="bernoulli", data=dat)
summary(DEKAY_brms)
mcmc_plot(DEKAY_brms, pars=c("b_logarea", "b_wetland_types", "b_forest_prop", "b_normCONN5"))
hyo

#northern water snake

dat$NO_WAT10<-ifelse(dat$NO_WAT_OBS>0,1,0) #creates new 1/0 column for NO_WAT
NO_WATmod<-glmer(NO_WAT10~logarea + forest_prop + wetland_types + normCONN5 + (1|City), 
               family="binomial", data=dat)
summary(NO_WATmod)

NO_WATocc<-mean(dat$NO_WAT10, na.rm=TRUE)


NO_WAT_brms<-brm(NO_WAT10~logarea + forest_prop + wetland_types + normCONN5 + (1|City), 
               prior=prior1, family="bernoulli", data=dat)
summary(NO_WAT_brms)
mcmc_plot(NO_WAT_brms, pars=c("b_logarea", "b_wetland_types", "b_forest_prop", "b_normCONN5"))
