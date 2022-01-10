#supplementary analyses
#these are analyses that we carried out but did not include in the manuscript.
#most were non-informative, some were checks to make sure we weren't missing
#anything in the primary set of analyses

#same data processing as main script

library(ggplot2)
library(lme4)
library(rstan)
library(brms)
rstan_options(auto_write = TRUE)

dat<-ALLCITIES8

#primary model but with flat prior (default in brms)
BRMflatprior<-brm(HERP_SP ~ logarea + wetland_types + forest_prop + normCONN5 + (1|City), 
          family = "negbinomial", data=dat)
summary(BRMflatprior)

#prior for coefficients
prior1<-set_prior("normal(0,10)", class="b")

#primary model but with 1000 m connectivity (instead of 500 m in other models)
BRMconn10<-brm(HERP_SP ~ logarea + wetland_types + forest_prop + normCONN10 + (1|City), 
          family = "negbinomial", prior=prior1, data=dat)
summary(BRMconn10)


#include a quadratic term for forest proportion
#(too much forest might crowd out reptiles)
BRMf2<-brm(HERP_SP ~ logarea + wetland_types + forest_prop + fp2 + normCONN5 + (1|City), 
           family = "negbinomial", prior=prior1, data=dat)
summary(BRMf2)
#random slope models - these mostly don't converge - uncertainty too high to do cross-city comparisons

#difference in area effect across cities - MAJOR CONVERGENCE ISSUES
BRM6<-brm(HERP_SP ~ logarea + wetland_types + forest_prop + normCONN5 + (logarea|City), 
          family = "negbinomial", prior=prior1, data=dat,
          control = list(adapt_delta = 0.95))
summary(BRM6)
ranef(BRM6)

fitted_values <- fitted(BRM6)
head(fitted_values)
pred_dat <- as.data.frame(cbind(Y = standata(BRM6)$Y, fitted_values))
pred_dat$dev<-abs(pred_dat$Y-pred_dat$Estimate)
mean(pred_dat$dev)  

#this has a random slope for connectivity (connectivity effect varies across cities)
BRM7<-brm(HERP_SP ~ normCONN5 + logarea + wetland_types + forest_prop + (normCONN5|City),
          family="negbinomial", prior=prior1, data=dat,
          control = list(adapt_delta = 0.95))
summary(BRM7)
ranef(BRM7)

fitted_values <- fitted(BRM7)
head(fitted_values)
pred_dat <- as.data.frame(cbind(Y = standata(BRM7)$Y, fitted_values))
pred_dat$dev<-abs(pred_dat$Y-pred_dat$Estimate)
mean(pred_dat$dev)  


#random slope for forest_prop
BRM8<-brm(HERP_SP ~ logarea + wetland_types + normCONN5 + forest_prop + (forest_prop|City),
          family="negbinomial", prior=prior1, data=dat,
          control = list(adapt_delta = 0.95))
summary(BRM8)
ranef(BRM8)

fitted_values <- fitted(BRM8)
head(fitted_values)
pred_dat <- as.data.frame(cbind(Y = standata(BRM8)$Y, fitted_values))
pred_dat$dev<-abs(pred_dat$Y-pred_dat$Estimate)
mean(pred_dat$dev)  

#random slope for wetland types
BRM9<-brm(HERP_SP ~ logarea + normCONN + wetland_types + forest_prop + (wetland_types|City),
          family="poisson", prior=prior1, data=dat,
          control = list(adapt_delta = 0.95))
summary(BRM9)
ranef(BRM9)

fitted_values <- fitted(BRM9)
head(fitted_values)
pred_dat <- as.data.frame(cbind(Y = standata(BRM9)$Y, fitted_values))
pred_dat$dev<-abs(pred_dat$Y-pred_dat$Estimate)
mean(pred_dat$dev)  


#run each variable against richness by itself:

BRMarea<-brm(HERP_SP~logarea + (1|City), prior=prior1, family="negbinomial", data=dat)
summary(BRMarea)

plot(marginal_effects(BRMarea))

#extract fitted values (predictions) from the model and compare to actual values
fitted_values <- fitted(BRMarea)
head(fitted_values)
pred_dat <- as.data.frame(cbind(Y = standata(BRMarea)$Y, fitted_values))
pred_dat$dev<-abs(pred_dat$Y-pred_dat$Estimate)
pred_dat$dev2<-(pred_dat$Y-pred_dat$Estimate)^2
mean(pred_dat$dev)  #this is just the mean diff between predicted and actual values
mean(pred_dat$dev2)

BRMfor<-brm(HERP_SP~forest_prop + (1|City), prior=prior1, family="negbinomial", data=dat)
summary(BRMfor)

#extract fitted values (predictions) from the model and compare to actual values
fitted_values <- fitted(BRMfor)
head(fitted_values)
pred_dat <- as.data.frame(cbind(Y = standata(BRMfor)$Y, fitted_values))
pred_dat$dev<-abs(pred_dat$Y-pred_dat$Estimate)
pred_dat$dev2<-(pred_dat$Y-pred_dat$Estimate)^2
mean(pred_dat$dev)  #this is just the mean diff between predicted and actual values
mean(pred_dat$dev2)

BRMCONN<-brm(HERP_SP~normCONN5 + (1|City), prior=prior1, family="negbinomial", data=dat)
summary(BRMCONN)

#extract fitted values (predictions) from the model and compare to actual values
fitted_values <- fitted(BRMCONN)
head(fitted_values)
pred_dat <- as.data.frame(cbind(Y = standata(BRMCONN)$Y, fitted_values))
pred_dat$dev<-abs(pred_dat$Y-pred_dat$Estimate)
pred_dat$dev2<-(pred_dat$Y-pred_dat$Estimate)^2
mean(pred_dat$dev)  #this is just the mean diff between predicted and actual values
mean(pred_dat$dev2)

BRMwet<-brm(HERP_SP~ wetland_types + (1|City), prior=prior1, family="negbinomial", data=dat)
summary(BRMwet)

#extract fitted values (predictions) from the model and compare to actual values
fitted_values <- fitted(BRMwet)
head(fitted_values)
pred_dat <- as.data.frame(cbind(Y = standata(BRMwet)$Y, fitted_values))
pred_dat$dev<-abs(pred_dat$Y-pred_dat$Estimate)
pred_dat$dev2<-(pred_dat$Y-pred_dat$Estimate)^2
mean(pred_dat$dev)  #this is just the mean diff between predicted and actual values
mean(pred_dat$dev2)

BRMimp<-brm(HERP_SP~ imperv_prop + (1|City), prior=prior1, family="negbinomial", data=dat)
summary(BRMimp)

#extract fitted values (predictions) from the model and compare to actual values
fitted_values <- fitted(BRMimp)
head(fitted_values)
pred_dat <- as.data.frame(cbind(Y = standata(BRMimp)$Y, fitted_values))
pred_dat$dev<-abs(pred_dat$Y-pred_dat$Estimate)
pred_dat$dev2<-(pred_dat$Y-pred_dat$Estimate)^2
mean(pred_dat$dev)  #this is just the mean diff between predicted and actual values
mean(pred_dat$dev2)

BRMroad<-brm(HERP_SP~ normROAD + (1|City), prior=prior1, family="negbinomial", data=dat)
summary(BRMroad)

#extract fitted values (predictions) from the model and compare to actual values
fitted_values <- fitted(BRMroad)
head(fitted_values)
pred_dat <- as.data.frame(cbind(Y = standata(BRMroad)$Y, fitted_values))
pred_dat$dev<-abs(pred_dat$Y-pred_dat$Estimate)
pred_dat$dev2<-(pred_dat$Y-pred_dat$Estimate)^2
mean(pred_dat$dev)  #this is just the mean diff between predicted and actual values
mean(pred_dat$dev2)

BRMshape<-brm(HERP_SP~ normSHAPE + (1|City), family="negbinomial", prior=prior1, data=dat)
summary(BRMshape)

#extract fitted values (predictions) from the model and compare to actual values
fitted_values <- fitted(BRMshape)
head(fitted_values)
pred_dat <- as.data.frame(cbind(Y = standata(BRMshape)$Y, fitted_values))
pred_dat$dev<-abs(pred_dat$Y-pred_dat$Estimate)
pred_dat$dev2<-(pred_dat$Y-pred_dat$Estimate)^2
mean(pred_dat$dev)  #this is just the mean diff between predicted and actual values
mean(pred_dat$dev2)

#univariate model for visitation rate
BRMvis1<-brm(HERP_SP~ normVIS + (1|City), family="negbinomial", prior=prior1, data=dat)
summary(BRMvis1)

#add data set with city level variables
dat<-ALLCITIESwithCITYLEVEL

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

#primary model but with city-level impervious surface 
BRMwithIMP<-brm(HERP_SP ~ logarea + wetland_types + forest_prop + normCONN5 + city_imperv + (1|City), 
               family = "negbinomial", prior=prior1, data=dat)
summary(BRMwithIMP)

#primary model but with city-level forest cover
BRMwithFOR<-brm(HERP_SP ~ logarea + wetland_types + forest_prop + normCONN5 + city_forest + (1|City), 
                family = "negbinomial", prior=prior1, data=dat)
summary(BRMwithFOR)






