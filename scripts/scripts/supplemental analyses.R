#supplementary analyses

#first set are analyses that were added in revision as sensitivity analyses (i.e)
#relaxing or added different assumptions

#second set are analyses that we carried out but did not include in the manuscript.
#most were non-informative, some were checks to make sure we weren't missing
#anything in the primary set of analyses

#same data processing as main script

library(ggplot2)
library(lme4)
library(rstan)
library(brms)
rstan_options(auto_write = TRUE)

dat<-ALLCITIES_revision


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

#added in revision - analysis of non-native species
dat$exoSP<-dat$exoAMPH+dat$exoREPS  #totals for exotic species
hist(dat$exoSP)

hist(dat$exoSP)
BRMexot<-brm(exoSP ~ logarea + wetland_types + forest_prop + normCONN5 + (1|City), 
             family = "negbinomial", prior=prior1, data=dat)
summary(BRMexot)

mcmc_plot(BRMexot, pars=c("b_logarea", "b_wetland_types", "b_forest_prop", "b_normCONN5"))

fitted_values <- fitted(BRMexot)
head(fitted_values)
pred_dat <- as.data.frame(cbind(Y = standata(BRMexot)$Y, fitted_values))
pred_dat$dev<-abs(pred_dat$Y-pred_dat$Estimate)
pred_dat$dev2<-(pred_dat$Y-pred_dat$Estimate)^2
mean(pred_dat$dev)  #this is just the mean diff between predicted and actual values
mean(pred_dat$dev2)

#for comparison here's a null model:
BRM0<-brm(exoSP ~ 1 + (1|City), family="negbinomial", data=dat)
summary(BRM0)

fitted_values <- fitted(BRM0)
head(fitted_values)
pred_dat <- as.data.frame(cbind(Y = standata(BRM0)$Y, fitted_values))
pred_dat$dev<-abs(pred_dat$Y-pred_dat$Estimate)
pred_dat$dev2<-(pred_dat$Y-pred_dat$Estimate)^2
mean(pred_dat$dev)  #this is just the mean diff between predicted and actual values
mean(pred_dat$dev2)

#added in revision - add variable for landscape vs. natural parks
dat$landscaped<-relevel(dat$landscaped, ref="NO")
BRMrev<-brm(HERP_SP ~ logarea + wetland_types + forest_prop + normCONN5 + landscaped + (1|City), 
            family = "negbinomial", prior=prior1, data=dat)
summary(BRMrev)

plot(BRMrev)
plot(conditional_effects(BRMrev))
mcmc_plot(BRMrev, pars=c("b_logarea", "b_wetland_types", "b_forest_prop", "b_normCONN5"))

plot(conditional_effects(BRMrev))

fitted_values <- fitted(BRMrev)
head(fitted_values)
pred_dat <- as.data.frame(cbind(Y = standata(BRMrev)$Y, fitted_values))
pred_dat$dev<-abs(pred_dat$Y-pred_dat$Estimate)
mean(pred_dat$dev)  #this is just the mean diff between predicted and actual values
pred_dat$dev2<-(pred_dat$Y-pred_dat$Estimate)^2  #SQUARED DEVIATION
mean(pred_dat$dev2)


#for comparison here's a null model:
BRM0<-brm(HERP_SP ~ 1 + (1|City), family="negbinomial", data=dat)
summary(BRM0)

fitted_values <- fitted(BRM0)
head(fitted_values)
pred_dat <- as.data.frame(cbind(Y = standata(BRM0)$Y, fitted_values))
pred_dat$dev<-abs(pred_dat$Y-pred_dat$Estimate)
pred_dat$dev2<-(pred_dat$Y-pred_dat$Estimate)^2
mean(pred_dat$dev)  #this is just the mean diff between predicted and actual values
mean(pred_dat$dev2)

#subsets for landscaped and unlandscapes parks
datNAT<-subset(dat, landscaped!="YES") #includes both unlandscaped and mixed
datSCAPE<-subset(dat, landscaped=="YES")

BRMnat<-brm(HERP_SP ~ logarea + wetland_types + forest_prop + normCONN5 + (1|City), 
            family = "negbinomial", prior=prior1, data=datNAT)
summary(BRMnat)


fitted_values <- fitted(BRMnat)
head(fitted_values)
pred_dat <- as.data.frame(cbind(Y = standata(BRMnat)$Y, fitted_values))
pred_dat$dev<-abs(pred_dat$Y-pred_dat$Estimate)
mean(pred_dat$dev)  #this is just the mean diff between predicted and actual values
pred_dat$dev2<-(pred_dat$Y-pred_dat$Estimate)^2  #SQUARED DEVIATION
mean(pred_dat$dev2)

#plot of predicted versus actual richness for primary model
ggplot(pred_dat, aes(Estimate, Y))+
  geom_point(size=3)+
  geom_abline(slope=1, intercept=0)+
  ggtitle("Herp Richness")+
  labs(x="Expected Species", y="Observed Species")+
  theme(plot.title = element_text(hjust = 0.5))

#compare to a null model

BRM0nat<-brm(HERP_SP ~ 1 + (1|City), family="negbinomial", data=datNAT)
summary(BRM0nat)

fitted_values <- fitted(BRM0nat)
head(fitted_values)
pred_dat <- as.data.frame(cbind(Y = standata(BRM0nat)$Y, fitted_values))
pred_dat$dev<-abs(pred_dat$Y-pred_dat$Estimate)
pred_dat$dev2<-(pred_dat$Y-pred_dat$Estimate)^2
mean(pred_dat$dev)  #this is just the mean diff between predicted and actual values
mean(pred_dat$dev2)

#plot of area vs richness across cities

ggplot(datNAT, aes(logarea, HERP_SP))+
  geom_point(size=3)+
  labs(x="Log Area", y="Herp Richness")+
  facet_wrap(~City)

#now repeat for landscaped sites

brmSCAPE<-brm(HERP_SP ~ logarea + wetland_types + forest_prop + normCONN5 + (1|City), 
              family = "negbinomial", prior=prior1, data=datSCAPE)
summary(brmSCAPE)


fitted_values <- fitted(brmSCAPE)
head(fitted_values)
pred_dat <- as.data.frame(cbind(Y = standata(brmSCAPE)$Y, fitted_values))
pred_dat$dev<-abs(pred_dat$Y-pred_dat$Estimate)
mean(pred_dat$dev)  #this is just the mean diff between predicted and actual values
pred_dat$dev2<-(pred_dat$Y-pred_dat$Estimate)^2  #SQUARED DEVIATION
mean(pred_dat$dev2)

#plot of predicted versus actual richness for primary model
ggplot(pred_dat, aes(Estimate, Y))+
  geom_point(size=3)+
  geom_abline(slope=1, intercept=0)+
  ggtitle("Herp Richness")+
  labs(x="Expected Species", y="Observed Species")+
  theme(plot.title = element_text(hjust = 0.5))

#comparison to null model


BRM0scape<-brm(HERP_SP ~ 1 + (1|City), family="negbinomial", data=datSCAPE)
summary(BRM0scape)

fitted_values <- fitted(BRM0scape)
head(fitted_values)
pred_dat <- as.data.frame(cbind(Y = standata(BRM0scape)$Y, fitted_values))
pred_dat$dev<-abs(pred_dat$Y-pred_dat$Estimate)
pred_dat$dev2<-(pred_dat$Y-pred_dat$Estimate)^2
mean(pred_dat$dev)  #this is just the mean diff between predicted and actual values
mean(pred_dat$dev2)

#plot of area vs richness across cities

ggplot(datSCAPE, aes(logarea, HERP_SP))+
  geom_point(size=3)+
  labs(x="Log Area", y="Herp Richness")+
  facet_wrap(~City)


#re-analyze with exotic species included


dat$HERP_TOT<-dat$HERP_SP+dat$exoREPS+dat$exoAMPH

BRMtot<-brm(HERP_TOT ~ logarea + wetland_types + forest_prop + normCONN5 + (1|City), 
            family = "negbinomial", prior=prior1, data=dat)
summary(BRMtot)

mcmc_plot(BRMtot, pars=c("b_logarea", "b_wetland_types", "b_forest_prop", "b_normCONN5"))

fitted_values <- fitted(BRMtot)
head(fitted_values)
pred_dat <- as.data.frame(cbind(Y = standata(BRMtot)$Y, fitted_values))
pred_dat$dev<-abs(pred_dat$Y-pred_dat$Estimate)
pred_dat$dev2<-(pred_dat$Y-pred_dat$Estimate)^2
mean(pred_dat$dev)  #this is just the mean diff between predicted and actual values
mean(pred_dat$dev2)

#for comparison here's a null model:
BRM0<-brm(HERP_TOT ~ 1 + (1|City), family="negbinomial", data=dat)
summary(BRM0)

fitted_values <- fitted(BRM0)
head(fitted_values)
pred_dat <- as.data.frame(cbind(Y = standata(BRM0)$Y, fitted_values))
pred_dat$dev<-abs(pred_dat$Y-pred_dat$Estimate)
pred_dat$dev2<-(pred_dat$Y-pred_dat$Estimate)^2
mean(pred_dat$dev)  #this is just the mean diff between predicted and actual values
mean(pred_dat$dev2)

#what happens if you use species pool as an offset in addition to the city random effect
BRMpool<-brm(HERP_SP ~ logarea + wetland_types + forest_prop + normCONN5 + offset(log(species.pool)) + (1|City), 
             family = "negbinomial", prior=prior1, data=dat)
summary(BRMpool)

fitted_values <- fitted(BRMpool)
head(fitted_values)
pred_dat <- as.data.frame(cbind(Y = standata(BRMpool)$Y, fitted_values))
pred_dat$dev<-abs(pred_dat$Y-pred_dat$Estimate)
mean(pred_dat$dev)  #this is just the mean diff between predicted and actual values
pred_dat$dev2<-(pred_dat$Y-pred_dat$Estimate)^2  #SQUARED DEVIATION
mean(pred_dat$dev2)

#plot of predicted versus actual richness for primary model
ggplot(pred_dat, aes(Estimate, Y))+
  geom_point(size=3)+
  geom_abline(slope=1, intercept=0)+
  ggtitle("Herp Richness")+
  labs(x="Expected Species", y="Observed Species")+
  theme(plot.title = element_text(hjust = 0.5))

#plot of area vs richness across cities

ggplot(dat, aes(logarea, HERP_SP))+
  geom_point(size=3)+
  labs(x="Log Area", y="Herp Richness")+
  facet_wrap(~City)

#for comparison here's a null model:
BRM0<-brm(HERP_SP ~ 1 + (1|City), family="negbinomial", data=dat)
summary(BRM0)

fitted_values <- fitted(BRM0)
head(fitted_values)
pred_dat <- as.data.frame(cbind(Y = standata(BRM0)$Y, fitted_values))
pred_dat$dev<-abs(pred_dat$Y-pred_dat$Estimate)
pred_dat$dev2<-(pred_dat$Y-pred_dat$Estimate)^2
mean(pred_dat$dev)  #this is just the mean diff between predicted and actual values
mean(pred_dat$dev2)

#primary model but with NYC as one city instead of 5 burroughs

BRMcity2<-brm(HERP_SP ~ logarea + wetland_types + forest_prop + normCONN5 + (1|City2), 
            family = "negbinomial", prior=prior1, data=dat)
summary(BRMcity2)


#additional analyses that were carried out to but not included in paper
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





