
#required packages
require(readr)
require(lme4)
require(MCMCglmm)
require(dplyr)
require(ggplot2)
require(ggeffects)
require(arm)

install.packages("prettydoc")
#load data-
data<-read.csv("data/nest_defense_test_final.csv")
names(data)

# Make year,ID, ID_Series factors
data$Year2<-as.factor(data$Year2)
data$ID_Series<-as.factor(data$ID_Series)
data$ID<-as.factor(data$ID)

summary(data$MinDis)
data$MinDis100<-ifelse(data$MinDisRaw>100,100,data$MinDisRaw)
data$LogMin<-log(data$MinDis100+1)

###1. Sources of variation in nest defense----
###first- do univariate analyses of each of the three measures to i) evaluate if all three behaviours are represatations of nest defense, and ii) choose best measure to use in main text
#make is such that minimum distance is capped at 100 meters (i.e., >100 is not a response)
data2<-subset(data, data$MinDisRaw<=100)
names(data2)

hist(data2$MinDisRaw)
data2$logMin<-log(data2$MinDisRaw+1)
hist(data2$logMin)

#Minimum distance---- 
require(lme4)
m1<-lmer(logMin~ Sex + NestStage + Year2+ (1|ID)+(1|Observer) ,
         data=data2) 
#summary
summary(m1)
#model check
plot(m1)
hist(resid(m1)) # residuals
lattice::qqmath(m1) #normality of errors

sm1<-sim(m1)
smfixef=sm1@fixef
smranef=sm1@ranef
smfixef=as.mcmc(smfixef)
posterior.mode(smfixef)
HPDinterval(smfixef)

bID<-sm1@ranef$ID
bvar<-as.vector(apply(bID, 1, var)) ##between individual variance posterior distribution
bvar<-as.mcmc(bvar)
posterior.mode(bvar )## mode of the distribution
HPDinterval(bvar)

bObs<-sm1@ranef$Observer
ObsVar<-as.vector(apply(bObs, 1, var)) ##between individual variance posterior distribution
ObsVar<-as.mcmc(ObsVar)
posterior.mode(ObsVar )## mode of the distribution
HPDinterval(ObsVar)

rvar<-sm1@sigma^2
rvar<-as.mcmc(rvar)
posterior.mode(rvar)
HPDinterval(rvar)

r <- bvar/(bvar + rvar)        ###individual
posterior.mode(r)
HPDinterval(r)


###Min distance (truncated to max of 100)----
data$MinDis100<-ifelse(data$MinDisRaw>100,100,data$MinDisRaw)
data$logMinDis100<-log(data$MinDis100+1)
hist(data$logMinDis100)
m1.1<-lmer(logMinDis100~ Sex + NestStage + Year2+ (1|ID)+(1|Observer) ,
         data=data) 
#summary
summary(m1.1)
#model check
plot(m1.1)
hist(resid(m1.1)) # residuals
lattice::qqmath(m1.1) #normality of errors

sm1<-sim(m1.1)
smfixef=sm1@fixef
smranef=sm1@ranef
smfixef=as.mcmc(smfixef)
posterior.mode(smfixef)
HPDinterval(smfixef)

bID<-sm1@ranef$ID
bvar<-as.vector(apply(bID, 1, var)) ##between individual variance posterior distribution
bvar<-as.mcmc(bvar)
posterior.mode(bvar )## mode of the distribution
HPDinterval(bvar)

bObs<-sm1@ranef$Observer
ObsVar<-as.vector(apply(bObs, 1, var)) ##between individual variance posterior distribution
ObsVar<-as.mcmc(ObsVar)
posterior.mode(ObsVar )## mode of the distribution
HPDinterval(ObsVar)

rvar<-sm1@sigma^2
rvar<-as.mcmc(rvar)
posterior.mode(rvar)
HPDinterval(rvar)

r <- bvar/(bvar + rvar)        ###individual
posterior.mode(r)
HPDinterval(r)

###Number of stoops----
names(data)
m2<-glmer(Dives~ Sex + NestStage + Year2+ (1|ID)+ (1|Observer) ,
         data=data, family="poisson") 
#summary
summary(m2)

sm2<-sim(m2)
smfixef=sm2@fixef
smranef=sm2@ranef
smfixef=as.mcmc(smfixef)
posterior.mode(smfixef)
HPDinterval(smfixef)

bID<-sm2@ranef$ID
bvar<-as.vector(apply(bID, 1, var)) ##between individual variance posterior distribution
bvar<-as.mcmc(bvar)
posterior.mode(bvar )## mode of the distribution
HPDinterval(bvar)

bObs<-sm2@ranef$Observer
ObsVar<-as.vector(apply(bObs, 1, var)) ##between individual variance posterior distribution
ObsVar<-as.mcmc(ObsVar)
posterior.mode(ObsVar )## mode of the distribution
HPDinterval(ObsVar)


r <- bvar/(bvar + 1)        ###individual
posterior.mode(r)
HPDinterval(r)

###FID----
data3<-subset(data, data$FID>0)
hist(data3$FID)
data3$logFID<-log(data3$FID+1)
hist(data3$logFID)

m3<-lmer(logFID~ Sex + NestStage + Year2+ (1|ID)+(1|Observer) ,
         data=data3) 
#summary
summary(m3)
#model check
plot(m3)
hist(resid(m3)) # residuals
lattice::qqmath(m3) #normality of errors

sm3<-sim(m3)
smfixef=sm3@fixef
smranef=sm3@ranef
smfixef=as.mcmc(smfixef)
posterior.mode(smfixef)
HPDinterval(smfixef)

bID<-sm3@ranef$ID
bvar<-as.vector(apply(bID, 1, var)) ##between individual variance posterior distribution
bvar<-as.mcmc(bvar)
posterior.mode(bvar )## mode of the distribution
HPDinterval(bvar)

bObs<-sm3@ranef$Observer
ObsVar<-as.vector(apply(bObs, 1, var)) ##between individual variance posterior distribution
ObsVar<-as.mcmc(ObsVar)
posterior.mode(ObsVar )## mode of the distribution
HPDinterval(ObsVar)

rvar<-sm1@sigma^2
rvar<-as.mcmc(rvar)
posterior.mode(rvar)
HPDinterval(rvar)

r <- bvar/(bvar + rvar)        ###individual
posterior.mode(r)
HPDinterval(r)




###2. Multivariate model with 9 nest defense traits (across contexts)----
data3c<-read.csv("data/nest_defense_test_final_20220318.csv")
names(data3c)

data_3c<-subset(data3c, data3c$FID>0)


data_3c$logLFID<-log(data_3c$LFID+1) 
data_3c$logIFID<-log(data_3c$IFID+1)
data_3c$logPFID<-log(data_3c$PFID+1)

data_3c$logLMinDis<-log(data_3c$LMinDis+1)
data_3c$logIMinDis<-log(data_3c$IMinDis+1)
data_3c$logPMinDis<-log(data_3c$PMinDis+1)




### Multivariate: FID across breeding contexts----
prior3var=list(R=list(V=diag(3),nu=3.002),G=list(G1=list(V=diag(3), nu=3.002)))

m3FID<-MCMCglmm(cbind(logLFID, logIFID,logPFID)~(trait-1),
             random=~us(trait):ID,rcov=~idh(trait):units,
             family=c("gaussian","gaussian", "gaussian"), prior=prior3var,
             nitt=1300000,thin=1000,burnin=300000, data=data_3c,verbose=TRUE)

summary(m3FID)
plot(m3FID)

#among-individual covariance
c1 <- posterior.cor(m3FID$VCV[,1:9])
round(apply(c1,2,mean),2)
round(apply(c1,2, quantile, c(0.025, 0.975)),2)


##Multivariate: minimum approach distance across breeding contexts----
m3Min<-MCMCglmm(cbind(logLMinDis, logIMinDis,logPMinDis)~(trait-1),
                random=~us(trait):ID,rcov=~idh(trait):units,
                family=c("gaussian","gaussian", "gaussian"), prior=prior3var,
                nitt=1300000,thin=1000,burnin=300000, data=data_3c,verbose=TRUE)

summary(m3Min)
plot(m3Min)

#among-individual covariance
c1 <- posterior.cor(m3Min$VCV[,1:9])
round(apply(c1,2,mean),2)
round(apply(c1,2, quantile, c(0.025, 0.975)),2)


#Multivariate: number of dives across breeding contexts----
m3Dives<-MCMCglmm(cbind(Ldives, Idives,Pdives)~(trait-1),
                random=~us(trait):ID,rcov=~idh(trait):units,
                family=c("poisson","poisson", "poisson"), prior=prior3var,
                nitt=1300000,thin=1000,burnin=300000, data=data_3c,verbose=TRUE)

summary(m3Dives)
plot(m3Dives)

#among-individual covariance
c1 <- posterior.cor(m3Dives$VCV[,1:9])
round(apply(c1,2,mean),2)
round(apply(c1,2, quantile, c(0.025, 0.975)),2)



#3. Bivariate model with minimum distance and number of dives----

prior_2var= list(R = list(V = diag(2), nu = 1.002),
                 G = list(G1 = list(V = diag(2), nu = 2,
                                    alpha.mu = rep(0,2),
                                    alpha.V = diag(25^2,2))))
names(data)

model.2var <- MCMCglmm(
  cbind(LogMin, Dives) ~ trait - 1,
  random = ~ us(trait):ID,
  rcov = ~ us(trait):units,
  family = c("gaussian", "poisson"),
  data = data,
  prior = prior_2var,
  verbose = TRUE,
  #nitt = 1300000, thin = 1000, burnin = 300000
  )

plot(model.2var)
summary(model.2var)

#among-individual covariance
c1 <- posterior.cor(model.2var$VCV[,1:4])
round(apply(c1,2,mean),2)
round(apply(c1,2, quantile, c(0.025, 0.975)),2)

#within individual covariance
c2 <- posterior.cor(model.2var$VCV[,5:8])
round(apply(c2,2,mean),2)
round(apply(c2,2, quantile, c(0.025, 0.975)),2)

#4. Univariate model to estimate short-versus long-term repeatability----
require(lme4)

#no sex:nest stage interaction... Table S3
m1.0<-lmer(LogMin~ (Sex + NestStage)^2 + Year2+ (1|ID) + (1|ID_Series) ,
         data=data) 
summary(m1.0)
anova(m1.0)

smod<-sim(m1.0,1000)
posterior.mode(as.mcmc(smod@fixef))
HPDinterval(as.mcmc(smod@fixef))
##Between individual variance
bID<-smod@ranef$ID[,,1]
bvar<-as.vector(apply(bID, 1, var)) ##ID variance posterior distribution
require(MCMCglmm)
bvar<-as.mcmc(bvar)
posterior.mode(bvar )## mode of the distribution
HPDinterval(bvar)
##Between individual variance, ID_Series
bID2<-smod@ranef$ID_Series[,,1]
bvar2<-as.vector(apply(bID2, 1, var)) ##ID_Series variance posterior distribution
require(MCMCglmm)
bvar2<-as.mcmc(bvar2)
posterior.mode(bvar2 )## mode of the distribution
HPDinterval(bvar2)
###residual variance
rvar<-smod@sigma^2
rvar<-as.mcmc(rvar)
posterior.mode(rvar)
HPDinterval(rvar)
### Long-term repeatability
r<-bvar/(rvar+bvar +bvar2 )
posterior.mode(r)
HPDinterval(r)  ##repeatability 
###Short-term
r1<-(bvar +bvar2  )/(rvar+bvar +bvar2 )
posterior.mode(r1)
HPDinterval(r1)  ##repeatability



#Main effects similar without including interactions. Results for main text Table 1----
m1<-lmer(LogMin~ Sex + NestStage + Year2+ (1|ID) + (1|ID_Series) ,
         data=data) 
#summary
summary(m1)
#model check
plot(m1.1)
hist(resid(m1.1)) # residuals
lattice::qqmath(m1.1) #normality of errors


#simulated parameters + repeatability
require(dplyr)
library(arm)
require(MCMCglmm)
smod<-sim(m1,1000)
posterior.mode(as.mcmc(smod@fixef))
HPDinterval(as.mcmc(smod@fixef))
##Between individual variance
bID<-smod@ranef$ID[,,1]
bvar<-as.vector(apply(bID, 1, var)) ##ID variance posterior distribution
require(MCMCglmm)
bvar<-as.mcmc(bvar)
posterior.mode(bvar )## mode of the distribution
HPDinterval(bvar)
##Between individual variance, ID_Series
bID2<-smod@ranef$ID_Series[,,1]
bvar2<-as.vector(apply(bID2, 1, var)) ##ID_Series variance posterior distribution
require(MCMCglmm)
bvar2<-as.mcmc(bvar2)
posterior.mode(bvar2 )## mode of the distribution
HPDinterval(bvar2)
###residual variance
rvar<-smod@sigma^2
rvar<-as.mcmc(rvar)
posterior.mode(rvar)
HPDinterval(rvar)
### Long-term repeatability
r<-bvar/(rvar+bvar +bvar2 )
posterior.mode(r)
HPDinterval(r)  ##repeatability 
###Short-term
r1<-(bvar +bvar2  )/(rvar+bvar +bvar2 )
posterior.mode(r1)
HPDinterval(r1)  ##repeatability

#Plot----
d18<-subset(data, data$Year=="2018")
d19<-subset(data, data$Year=="2019")

#plot
p18<-ggplot(d18, aes(x = NestStage, y = MinDisRaw, factor=NestStage, fill=Sex, color=Sex)) + geom_boxplot() +
  scale_y_continuous("Nest Defense (Cube-root transformed)",limits = c(0, 100)) +
  scale_x_discrete("Nest Stage", labels=c('Egg-laying', 'Incubating', 'Provisioning'))+ 
  guides(fill="none")+
  guides(color="none")+
  theme_classic() 

p18

p19<-ggplot(d19, aes(x = NestStage, y = MinDisRaw, factor=NestStage, fill=Sex, color=Sex)) + geom_boxplot() +
  scale_y_continuous("",limits = c(0, 100)) +
  scale_x_discrete("Nest Stage", labels=c('Egg-laying', 'Incubating', 'Provisioning')) +
  theme_classic() 

p19


require(cowplot)
plot_grid(p18,p19, labels = "AUTO")

###5. Models of (dis-)assortative mating + relative fitness----
#data
data2018<-read.csv("data/data_2018.csv")
data2019<-read.csv("data/data_2019.csv")

names(data2018)
data2018$Male<-log(data2018$male_raw+1)*-1
data2018$Female<-log(data2018$female_raw+1)*-1

names(data2019)
data2019$Male<-log(data2019$male_raw+1)*-1
data2019$Female<-log(data2019$female_raw+1)*-1

#prior for three traits- from Houslay tutorial
prior_E_B_fit_1px = list(R = list(V = diag(c(1,1,0.0001),3,3), nu = 1.002, fix = 3),
                         G = list(G1 = list(V = diag(3), nu = 3,
                                            alpha.mu = rep(0,3),
                                            alpha.V = diag(25^2,3,3))))
#set residual variance for fitness near 0
prior_E_B_fit_1px$R$V[3,3]<-0.0001

#2018 3 trait model ----
mod.12 <- MCMCglmm(cbind(Male, Female, rfit) ~ (trait-1),  
                   random = ~us(trait):NestID ,
                   rcov = ~us(trait):units, 
                   family = c("gaussian", "gaussian", "gaussian"),
                   data=data2018, 
                   prior = prior_E_B_fit_1px, 
                   verbose = TRUE,
                   nitt=590000,thin=500,burnin=90000
                   )

plot(mod.12)

summary(mod.12)

#auto-correlation
autocorr.diag(mod.12$VCV)
autocorr(mod.12$Sol)


# posteriors
posteriors_3<-as.mcmc(mod.12$VCV)
posterior.mode(posteriors_3)

# among-pair correlations
pair.correlation_3<-posteriors_3[,"traitFemale:traitMale.NestID"]/
  sqrt(posteriors_3[,"traitFemale:traitFemale.NestID"]*
         posteriors_3[,"traitMale:traitMale.NestID"])

posterior.mode(pair.correlation_3)
HPDinterval(pair.correlation_3)
# within-pair correlations (residuals)
residual.correlation_3<-posteriors_3[,"traitFemale:traitMale.units"]/
  sqrt(posteriors_3[,"traitFemale:traitFemale.units"]*
         posteriors_3[,"traitMale:traitMale.units"])

posterior.mode(residual.correlation_3)
HPDinterval(residual.correlation_3)
#male selection gradient
Male_sel_18<- posteriors_3[,"traitrfit:traitMale.NestID"]/
  (sqrt(posteriors_3[,"traitrfit:traitrfit.NestID"])*
     sqrt(posteriors_3[,"traitMale:traitMale.NestID"]))
posterior.mode(Male_sel_18)
HPDinterval(Male_sel_18)

#female selection gradient
Female_sel_18<- posteriors_3[,"traitrfit:traitFemale.NestID"]/
  (sqrt(posteriors_3[,"traitrfit:traitrfit.NestID"])*
     sqrt(posteriors_3[,"traitFemale:traitFemale.NestID"]))
posterior.mode(Female_sel_18)
HPDinterval(Female_sel_18)


#2019 three trait  model

#prior for three trait- Houslay tutorial 
prior_E_B_fit_1px = list(R = list(V = diag(c(1,1,0.0001),3,3), nu = 1.002, fix = 3),
                         G = list(G1 = list(V = diag(3), nu = 3,
                                            alpha.mu = rep(0,3),
                                            alpha.V = diag(25^2,3,3))))

#model
mod.13 <- MCMCglmm(cbind(Male, Female, rfit) ~ (trait-1),  
                   random = ~us(trait):NestID ,
                   rcov = ~us(trait):units, 
                   family = c("gaussian", "gaussian", "gaussian"),
                   data=data2019, 
                   prior = prior_E_B_fit_1px, 
                   verbose = TRUE,
                   nitt=590000,thin=500,burnin=90000)

plot(mod.13)
summary(mod.13)

#auto-correlation
autocorr.diag(mod.13$VCV)
autocorr(mod.13$Sol)

# posteriors
posteriors1<-as.mcmc(mod.13$VCV)
posterior.mode(posteriors1)

# among-pair correlations
pair.correlation_4<-posteriors1[,"traitFemale:traitMale.NestID"]/
  sqrt(posteriors1[,"traitFemale:traitFemale.NestID"]*
         posteriors1[,"traitMale:traitMale.NestID"])

posterior.mode(pair.correlation_4)
HPDinterval(pair.correlation_4)
# within-pair correlations (residuals)
residual.correlation4<-posteriors1[,"traitFemale:traitMale.units"]/
  sqrt(posteriors1[,"traitFemale:traitFemale.units"]*
         posteriors1[,"traitMale:traitMale.units"])

posterior.mode(residual.correlation4)
HPDinterval(residual.correlation4)
#male selection gradient
Male_sel_19<- posteriors1[,"traitrfit:traitMale.NestID"]/
  (sqrt(posteriors1[,"traitrfit:traitrfit.NestID"])*
     sqrt(posteriors1[,"traitMale:traitMale.NestID"]))
posterior.mode(Male_sel_19)
HPDinterval(Male_sel_19)
#female selection gradient
Female_sel_19<- posteriors1[,"traitrfit:traitFemale.NestID"]/
  (sqrt(posteriors1[,"traitrfit:traitrfit.NestID"])*
     sqrt(posteriors1[,"traitFemale:traitFemale.NestID"]))
posterior.mode(Female_sel_19)
HPDinterval(Female_sel_19)


###comparing 2018 and 2019
#among
y2018a<-ifelse(pair.correlation_3<0.04291904,1,0)
sum(y2018a)

y2018b<-ifelse(pair.correlation_3>0.04291904,1,0)
sum(y2018b)

y2019a<-ifelse(pair.correlation_4>(-0.334157 ),1,0)
sum(y2019a)

y2019b<-ifelse(pair.correlation_4<(-0.334157 ),1,0)
sum(y2019b)

p<-(sum(y2018a)/1000)*(sum(y2019a)/1000)
p

#within
y2018w<-ifelse(residual.correlation_3>0.04291904,1,0)
sum(y2018w)
y2018x<-ifelse(residual.correlation_3<0.04291904,1,0)
sum(y2018x)
y2019w<-ifelse(residual.correlation4<0.3588893,1,0)
sum(y2019w)
y2019x<-ifelse(residual.correlation4>0.3588893,1,0)
sum(y2019x)

p<-sum(y2018w)/1000*sum(y2019w)/1000
p

##compare sexes within years
m2018<-ifelse(Male_sel_18<0.09374586,1,0)
m2018b<-ifelse(Male_sel_18<0.09374586,0,1)
f2018<-ifelse(Female_sel_18>(-0.4853156),1,0)
f2018b<-ifelse(Female_sel_18>(-0.4853156),0,1)

sum(m2018)
sum(m2018b)
sum(f2018)
sum(f2018b)

p<-sum(m2018)/1000*sum(f2018)/1000
p


m2019<-ifelse(Male_sel_19>-0.2069016,1,0)
m2019b<-ifelse(Male_sel_19>-0.2069016,0,1)
f2019<-ifelse(Female_sel_19<0.08781835,1,0)
f2019b<-ifelse(Female_sel_19<0.08781835,0,1)

sum(m2019)
sum(m2019b)
sum(f2019)
sum(f2019b)

p<-sum(m2019)/1000*sum(f2019)/1000
p

##compare within sex across years
#females
f2018<-ifelse(Female_sel_18>(-0.2069016),1,0)
f2018b<-ifelse(Female_sel_18>(-0.2069016),0,1)
f2019<-ifelse(Female_sel_19<0.09374586,1,0)
f2019b<-ifelse(Female_sel_19<0.09374586,0,1)

sum(f2018)
sum(f2018b)
sum(f2019)
sum(f2019b)

p<-sum(f2018)/1000*sum(f2019)/1000
p

#males
m2018<-ifelse(Male_sel_18<(0.08781835),1,0)
m2018b<-ifelse(Male_sel_18<(0.08781835),0,1)
m2019<-ifelse(Male_sel_19>(-0.4853156),1,0)
m2019b<-ifelse(Male_sel_19>(-0.4853156),0,1)

sum(m2018)
sum(m2018b)
sum(m2019)
sum(m2019b)

p<-sum(m2018)/1000*sum(m2019)/1000
p

### create data frame and within-subject center for plot


#create data frame subtracting subjects mean from each observation (i.e., cw2= within-subject centered)
require(dplyr)
d2018 <- data2018 %>% group_by(NestID) %>% 
  mutate(
    male_mean = mean(Male, na.rm = T),
    female_mean= mean(Female, na.rm = T)
  ) %>% ungroup() %>% 
  mutate(
    male_c2 = scale(Male, center = T, scale = F),
    female_c2 = scale(Female, center = T, scale = F),
    
    male_cw2 = Male - male_mean,
    female_cw2 = Female - female_mean,
    
    male_cb2 = male_c2 - male_cw2,
    female_cb2 = female_c2 - female_cw2
  )

d2019 <- data2019 %>% group_by(NestID) %>% 
  mutate(
    male_mean = mean(Male, na.rm = T),
    female_mean= mean(Female, na.rm = T)
  ) %>% ungroup() %>% 
  mutate(
    male_c2 = scale(Male, center = T, scale = F),
    female_c2 = scale(Female, center = T, scale = F),
    
    male_cw2 = Male - male_mean,
    female_cw2 = Female - female_mean,
    
    male_cb2 = male_c2 - male_cw2,
    female_cb2 = female_c2 - female_cw2
  )


#load package
require(ggplot2)

# Calculate regression lines from the model fit (co)variances
between_slope2018 <- mean(mod.12$VCV[,"traitMale:traitFemale.NestID"]/
                            mod.12$VCV[,"traitMale:traitMale.NestID"])
between_slope2019 <- mean(mod.13$VCV[,"traitMale:traitFemale.NestID"]/
                            mod.13$VCV[,"traitMale:traitMale.NestID"])
within_slope2018 <- mean(mod.12$VCV[,"traitMale:traitFemale.units"]/
                           mod.12$VCV[,"traitMale:traitMale.units"])
within_slope2019 <- mean(mod.13$VCV[,"traitMale:traitFemale.units"]/
                           mod.13$VCV[,"traitMale:traitMale.units"])



#Within-nest plot 2018
Within_nest2018<-ggplot(d2018, aes(x = male_cw2, y = female_cw2, group = NestID)) +
  geom_point(alpha = 1) + geom_abline(intercept =0 , slope = within_slope2018,size=1)+
  labs(x = "",
       y = "") +  scale_x_continuous(limits = c(-2,2)) + 
  scale_y_continuous( limits = c(-2,2)) +
  theme_classic() +  ggtitle("") + theme(plot.title = element_text(hjust=0.5)) + 
  theme(axis.title.y =element_text(size=11))



Within_nest2018

#within nest plot 2019

coef(lm(female_cw2 ~ male_cw2, data = d2019))
Within_nest2019<-ggplot(d2019, aes(x = male_cw2, y = female_cw2, group = NestID)) +
  geom_point(alpha = 1) + geom_abline(intercept =0 , slope = within_slope2019,size=1)+
  labs(x = "",
       y = "") +  scale_x_continuous(limits = c(-2,2)) + 
  scale_y_continuous( limits = c(-2,2)) +
  theme_classic() +  ggtitle("") + theme(plot.title = element_text(hjust=0.5)) + 
  theme(axis.title.y =element_text(size=11))



Within_nest2019

require(gridExtra)
wnp<-grid.arrange(Within_nest2018, Within_nest2019, ncol=2)

#Between-nest plot 2018


Between_nest2018<-ggplot(d2018, aes(x = male_cb2, y = female_cb2, group = NestID)) +
  geom_point(alpha = 1) + geom_abline(intercept = 0, slope =between_slope2018, size=1)+
  labs(x = " ",
       y = " ") +  scale_x_continuous(limits = c(-2.5,2)) + 
  scale_y_continuous( limits = c(-2,2)) +
  theme_classic() + ggtitle("") + theme(plot.title = element_text(hjust=0.5)) + 
  theme(axis.title.y =element_text(size=11))


Between_nest2018


#Between-nest plot 2019

Between_nest2019<-ggplot(d2019, aes(x = male_cb2, y = female_cb2, group = NestID)) +
  geom_point(alpha = 1) + geom_abline(intercept = 0, slope =between_slope2019, size=1)+
  labs(x = " ",
       y = " ") +  scale_x_continuous(limits = c(-2,2)) + 
  scale_y_continuous( limits = c(-2,2)) +
  theme_classic() +  ggtitle("") + theme(plot.title = element_text(hjust=0.5)) + 
  theme(axis.title.y =element_text(size=11))



Between_nest2019



bnp<-grid.arrange(Between_nest2018, Between_nest2019, ncol=2)

#combine plots
require(gridExtra)
require(cowplot)

g4<-grid.arrange(Between_nest2018, top= "Between-nest, 2018")
g5<-grid.arrange(Between_nest2019, top="Between-nest, 2019")
g6<-grid.arrange(Within_nest2018, top="Within-nest, 2018")
g7<-grid.arrange(Within_nest2019, top="Within-nest, 2019")


#Plot----

final_plot<-grid.arrange(g4,g6,g5,g7, left= "Female nest defense (log x + 1)", 
                         bottom= "Male nest defense (log x + 1)")
plot(final_plot)
#save plot
ggsave("Within_Between_centered.jpg", plot = final_plot, device= "jpeg", units="mm", dpi= "retina")



require(tidyr)

df_mcmc_cors <- data_frame(Traits = c("Among pair-correlation (2018)",    
                                      "Within pair-correlation (2018)",
                                      "Among pair-correlation (2019)",
                                      "Within pair-correlation (2019)"
                                        ),
                           Estimate = c(posterior.mode(pair.correlation_3),
                                        posterior.mode(residual.correlation_3),
                                        posterior.mode(pair.correlation_4),
                                        posterior.mode(residual.correlation4)),
                           Lower = c(HPDinterval(pair.correlation_3)[,"lower"],
                                     HPDinterval(residual.correlation_3)[,"lower"],
                                     HPDinterval(pair.correlation_4)[,"lower"],
                                     HPDinterval(residual.correlation4)[,"lower"]
                                     ),
                           Upper = c(HPDinterval(pair.correlation_3)[,"upper"],
                                     HPDinterval(residual.correlation_3)[,"upper"],
                                     HPDinterval(pair.correlation_4)[,"upper"],
                                     HPDinterval(residual.correlation4)[,"upper"]))




ggplot(df_mcmc_cors, aes(x = Traits, y = Estimate)) +
  geom_pointrange(aes(ymin = Lower,
                      ymax = Upper)) +
  geom_hline(yintercept = 0,
             linetype = "dotted",alpha = 0.3) +
  labs(x = "",
       y = "Correlation (Estimate +/- 95% CIs)") +
  ylim(-1,1) +
  coord_flip() +
  theme_classic()






df_mcmc_cors1 <- data_frame(Traits = c("Female (2018)",    
                                      "Male (2018)",
                                      "Female (2019)",
                                      "Male (2019)"
),
Estimate = c(posterior.mode(Female_sel_18),
             posterior.mode(Male_sel_18),
             posterior.mode(Female_sel_19),
             posterior.mode(Male_sel_19)),
Lower = c(HPDinterval(Female_sel_18)[,"lower"],
          HPDinterval(Male_sel_18)[,"lower"],
          HPDinterval(Female_sel_19)[,"lower"],
          HPDinterval(Male_sel_19)[,"lower"]
),
Upper = c(HPDinterval(Female_sel_18)[,"upper"],
          HPDinterval(Male_sel_18)[,"upper"],
          HPDinterval(Female_sel_19)[,"upper"],
          HPDinterval(Male_sel_19)[,"upper"]))




ggplot(df_mcmc_cors1, aes(x = Traits, y = Estimate)) +
  geom_pointrange(aes(ymin = Lower,
                      ymax = Upper)) +
  geom_hline(yintercept = 0,
             linetype = "dotted",alpha = 0.3) +
  labs(x = "",
       y = "Correlation (Estimate +/- 95% CIs)") +
  ylim(-1,1) +
  coord_flip() +
  theme_classic()
