#this is the main set of analyses of species richness (herps, amphibians..
#reptiles, uncommon species, and null model)

#predictive analyses for individual cities are in the "LOObycity" script

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

#CHECK OUT CORRELATIONS AMONG VARIABLES:

cordat<-cbind(dat$logarea, dat$forest_prop, dat$normCONN5, dat$normCONN10, 
              dat$imperv_prop, dat$normROAD, dat$normSHAPE, dat$wetland_types)
colnames(cordat)<-c("logarea", "forest_prop", "normCONN5", "normCONN10", 
                    "imperv_prop", "normROAD", "normSHAPE", "wetland_types")
cordat
cor(cordat)

#correlation among months in visitation data

cor(dat$Junevisits, dat$Julyvisits, use="pairwise.complete.obs")

#prior for coefficients
prior1<-set_prior("normal(0,10)", class="b")


#full model with respect to predictor variables
BRMfull<-brm(HERP_SP ~ logarea + wetland_types + forest_prop + imperv_prop + normCONN5 + normSHAPE + normROAD + (1|City), 
             family = "negbinomial", prior=prior1, data=dat)
summary(BRMfull, prob=0.90)


#extract fitted values (predictions) from the model and compare to actual values
fitted_values <- fitted(BRMfull)
head(fitted_values)
pred_dat <- as.data.frame(cbind(Y = standata(BRMfull)$Y, fitted_values))
pred_dat$dev<-abs(pred_dat$Y-pred_dat$Estimate)
pred_dat$dev2<-(pred_dat$Y-pred_dat$Estimate)^2
mean(pred_dat$dev)  #this is just the mean diff between predicted and actual values
mean(pred_dat$dev2)

#NOW DROP IMPERV, SHAPE, AND ROAD AND THIS IS THE BASIC MODEL

BRM1<-brm(HERP_SP ~ logarea + wetland_types + forest_prop + normCONN5 + (1|City), 
          family = "negbinomial", prior=prior1, data=dat)
summary(BRM1)

plot(BRM1)
plot(conditional_effects(BRM1))
mcmc_plot(BRM1, pars=c("b_logarea", "b_wetland_types", "b_forest_prop", "b_normCONN5"))

plot(conditional_effects(BRM1))

fitted_values <- fitted(BRM1)
head(fitted_values)
pred_dat <- as.data.frame(cbind(Y = standata(BRM1)$Y, fitted_values))
pred_dat$dev<-abs(pred_dat$Y-pred_dat$Estimate)
mean(pred_dat$dev)  #this is just the mean diff between predicted and actual values
pred_dat$dev2<-(pred_dat$Y-pred_dat$Estimate)^2  #SQUARED DEVIATION
mean(pred_dat$dev2)

#plot of predicted versus actual richness for primary model
#THIS IS FIGURE 2
ggplot(pred_dat, aes(Estimate, Y))+
  geom_point(size=3)+
  geom_abline(slope=1, intercept=0)+
  ggtitle("Herp Richness")+
  labs(x="Expected Species", y="Observed Species")+
  theme(plot.title = element_text(hjust = 0.5))

#plot of area vs richness across cities
#this is figure 1

ggplot(dat, aes(logarea, HERP_SP))+
  geom_point(size=3)+
  labs(x="Log Area", y="Herp Richness")+
  facet_wrap(~City)

#for comparison here's a null model:
BRM0<-brm(HERP_SP ~ 1 + (1|City), family="negbinomial", prior=prior1, data=dat)
summary(BRM0)

fitted_values <- fitted(BRM0)
head(fitted_values)
pred_dat <- as.data.frame(cbind(Y = standata(BRM0)$Y, fitted_values))
pred_dat$dev<-abs(pred_dat$Y-pred_dat$Estimate)
pred_dat$dev2<-(pred_dat$Y-pred_dat$Estimate)^2
mean(pred_dat$dev)  #this is just the mean diff between predicted and actual values
mean(pred_dat$dev2)

#THIS ONE INCLUDES VISITATION RATE, BUT ONLY A SUBSET OF THE DATA
BRMvis<-brm(HERP_SP ~ logarea + wetland_types + forest_prop + normCONN5 + normVIS + (1|City), 
            family = "negbinomial", prior=prior1, data=dat)
summary(BRMvis)

hypothesis(BRMvis, "normVIS<0")

#Model1 but with amphibians only - strong effect of forest_prop
BRMamph<-brm(AMPH_SP ~ logarea + wetland_types + forest_prop + normCONN5 + (1|City), 
             family = "negbinomial", prior=prior1, data=dat,
             control=list(adapt_delta=0.95))
summary(BRMamph)

mcmc_plot(BRMamph, pars=c("b_logarea", "b_wetland_types", "b_forest_prop", "b_normCONN5"))

fitted_values <- fitted(BRMamph)
head(fitted_values)
pred_dat <- as.data.frame(cbind(Y = standata(BRMamph)$Y, fitted_values))
pred_dat$dev<-abs(pred_dat$Y-pred_dat$Estimate)
pred_dat$dev2<-(pred_dat$Y-pred_dat$Estimate)^2
mean(pred_dat$dev)  #this is just the mean diff between predicted and actual values
mean(pred_dat$dev2)

ggplot(pred_dat, aes(Estimate, Y))+
  geom_point(size=2)+
  geom_abline(slope=1, intercept=0)+
  ggtitle("AMPHIBIANS")+
  labs(x="expected species", y="actual species")+
  theme(plot.title = element_text(hjust = 0.5))

#for comparison, null model for amphibians
AMPHNULL<-brm(AMPH_SP ~ 1 + (1|City), family="negbinomial", prior=prior1, data=dat)
summary(AMPHNULL)

fitted_values <- fitted(AMPHNULL)
head(fitted_values)
pred_dat <- as.data.frame(cbind(Y = standata(AMPHNULL)$Y, fitted_values))
pred_dat$dev<-abs(pred_dat$Y-pred_dat$Estimate)
pred_dat$dev2<-(pred_dat$Y-pred_dat$Estimate)^2
mean(pred_dat$dev)  #this is just the mean diff between predicted and actual values
mean(pred_dat$dev2)

#Model with reptiles only
BRMrep<-brm(REP_SP ~ logarea + wetland_types + forest_prop + normCONN5 + (1|City), 
            family = "negbinomial", prior=prior1, data=dat)
summary(BRMrep)

mcmc_plot(BRMrep, pars=c("b_logarea", "b_wetland_types", "b_forest_prop", "b_normCONN5"))

REPNULL<-brm(REP_SP ~ 1 + (1|City), family="negbinomial", prior=prior1, data=dat)
summary(REPNULL)

fitted_values <- fitted(REPNULL)
head(fitted_values)
pred_dat <- as.data.frame(cbind(Y = standata(REPNULL)$Y, fitted_values))
pred_dat$dev<-abs(pred_dat$Y-pred_dat$Estimate)
pred_dat$dev2<-(pred_dat$Y-pred_dat$Estimate)^2
mean(pred_dat$dev)  #this is just the mean diff between predicted and actual values
mean(pred_dat$dev2)

fitted_values <- fitted(BRMrep)
head(fitted_values)
pred_dat <- as.data.frame(cbind(Y = standata(BRMrep)$Y, fitted_values))
pred_dat$dev<-abs(pred_dat$Y-pred_dat$Estimate)
pred_dat$dev2<-(pred_dat$Y-pred_dat$Estimate)^2
mean(pred_dat$dev)  #this is just the mean diff between predicted and actual values
mean(pred_dat$dev2)

ggplot(pred_dat, aes(Estimate, Y))+
  geom_point(size=2)+
  geom_abline(slope=1, intercept=0)+
  ggtitle("REPTILES") +
  labs(x="expected species", y="actual species")+
  theme(plot.title = element_text(hjust = 0.5))

#MODEL1 BUT JUST UNCOMMON SPECIES

hist(dat$RARE_SP)
BRMrare<-brm(RARE_SP ~ logarea + wetland_types + forest_prop + normCONN5 + (1|City), 
             family = "negbinomial", prior=prior1, data=dat)
summary(BRMrare)

mcmc_plot(BRMrare, pars=c("b_logarea", "b_wetland_types", "b_forest_prop", "b_normCONN5"))

fitted_values <- fitted(BRMrare)
head(fitted_values)
pred_dat <- as.data.frame(cbind(Y = standata(BRMrare)$Y, fitted_values))
pred_dat$dev<-abs(pred_dat$Y-pred_dat$Estimate)
pred_dat$dev2<-(pred_dat$Y-pred_dat$Estimate)^2
mean(pred_dat$dev)  #this is just the mean diff between predicted and actual values
mean(pred_dat$dev2)

ggplot(pred_dat, aes(Estimate, Y))+
  geom_point(size=2)+
  geom_abline(slope=1, intercept=0)+
  ggtitle("RARE SPECIES") +
  labs(x="expected species", y="actual species")+
  theme(plot.title = element_text(hjust = 0.5))

