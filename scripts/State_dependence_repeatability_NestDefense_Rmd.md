R Script: Does fluctuating selection maintain variation in nest defense
behavior in Arctic peregrine falcons (Falco peregrinus tundrius)?
================

``` r
#required packages
require(readr)
require(lme4)
require(MCMCglmm)
require(dplyr)
require(ggplot2)
require(rlang)
require(tidyverse)
require(rlang)
require(ggeffects)
require(arm)
require(bayestestR)

# Set wd
setwd("C:/Users/nickg/OneDrive/Desktop/R projects/krmp_nest-defense")
#load data
data<-read.csv("C:/Users/nickg/OneDrive/Desktop/R projects/krmp_nest-defense/data/nest_defense_test_final.csv")
names(data)
```

#### **Make Year,ID, and ID_Series factors. Inspect raw data and transform.**

``` r
data$Year2<-as.factor(data$Year2)
data$ID_Series<-as.factor(data$ID_Series)
data$ID<-as.factor(data$ID)

summary(data$MinDisRaw)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##    0.10    5.00   15.00   27.72   25.00  600.00

``` r
data$MinDis100<-ifelse(data$MinDisRaw>100,100,data$MinDisRaw)
data$LogMin<-log(data$MinDis100+1)
```

#### **Minimum distance is capped at 100 meters (i.e., \>100 is not a response)**

``` r
data2<-subset(data, data$MinDisRaw<=100)
hist(data2$MinDisRaw)
```

![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/subset-1.png)<!-- -->

``` r
data2$logMin<-log(data2$MinDisRaw+1)
hist(data2$logMin)
```

![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/subset-2.png)<!-- -->

### **1. Sources of variation in nest defense**

-   *Univariate analyses of each of the three measures to i) evaluate if
    all three behaviors are representations of nest defense, and ii)
    choose best measure to use in main text.*

#### **Univariate model– Log-transformed minimum distance**

``` r
#Minimum distance---- 
require(lme4)
m1<-lmer(logMin~ Sex + NestStage + Year2+ (1|ID)+(1|Observer) ,
         data=data2) 

#summary
summary(m1)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: logMin ~ Sex + NestStage + Year2 + (1 | ID) + (1 | Observer)
    ##    Data: data2
    ## 
    ## REML criterion at convergence: 883.8
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.46859 -0.64241 -0.00094  0.59131  2.53522 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  ID       (Intercept) 0.29899  0.5468  
    ##  Observer (Intercept) 0.03097  0.1760  
    ##  Residual             0.52508  0.7246  
    ## Number of obs: 351, groups:  ID, 103; Observer, 8
    ## 
    ## Fixed effects:
    ##                        Estimate Std. Error t value
    ## (Intercept)             2.97079    0.15338  19.369
    ## SexMale                -0.30939    0.14235  -2.173
    ## NestStage2Incubating   -0.24865    0.11971  -2.077
    ## NestStage3Provisioning -0.56368    0.12764  -4.416
    ## Year21                  0.08779    0.12461   0.704
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) SexMal NstS2I NstS3P
    ## SexMale     -0.437                     
    ## NstStg2Incb -0.453  0.038              
    ## NstStg3Prvs -0.464  0.040  0.622       
    ## Year21      -0.409 -0.014  0.008  0.166

``` r
#model check
plot(m1)
```

![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/mod1-1.png)<!-- -->

``` r
hist(resid(m1)) # residuals
```

![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/mod1-2.png)<!-- -->

``` r
lattice::qqmath(m1) #normality of errors
```

![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/mod1-3.png)<!-- -->

#### **Simulate posterior distribution–Log-transformed minimum distance**

``` r
require(MCMCglmm)
require(arm)
require(lme4)
sm1<-sim(m1)
sm1<-as.mcmc(sm1)
smfixef=sm1@fixef

smranef=sm1@ranef
smfixef=as.mcmc(smfixef)
posterior.mode(smfixef)
```

    ##            (Intercept)                SexMale   NestStage2Incubating 
    ##              3.1144340             -0.3092407             -0.1785240 
    ## NestStage3Provisioning                 Year21 
    ##             -0.5276880              0.1121864

``` r
HPDinterval(smfixef)
```

    ##                             lower       upper
    ## (Intercept)             2.7122457  3.31641560
    ## SexMale                -0.5954352 -0.04671351
    ## NestStage2Incubating   -0.4652405  0.08820849
    ## NestStage3Provisioning -0.8017719 -0.29522147
    ## Year21                 -0.1636502  0.30812800
    ## attr(,"Probability")
    ## [1] 0.95

``` r
bID<-sm1@ranef$ID
bvar<-as.vector(apply(bID, 1, var)) ##between individual variance posterior distribution
bvar<-as.mcmc(bvar)
posterior.mode(bvar )## mode of the distribution
```

    ##      var1 
    ## 0.2976508

``` r
HPDinterval(bvar)
```

    ##          lower     upper
    ## var1 0.2420548 0.3741901
    ## attr(,"Probability")
    ## [1] 0.95

``` r
bObs<-sm1@ranef$Observer
ObsVar<-as.vector(apply(bObs, 1, var)) ##between individual variance posterior distribution
ObsVar<-as.mcmc(ObsVar)
posterior.mode(ObsVar )## mode of the distribution
```

    ##       var1 
    ## 0.02469483

``` r
HPDinterval(ObsVar)
```

    ##            lower      upper
    ## var1 0.004991606 0.05285004
    ## attr(,"Probability")
    ## [1] 0.95

``` r
rvar<-sm1@sigma^2
rvar<-as.mcmc(rvar)
posterior.mode(rvar)
```

    ##      var1 
    ## 0.5335908

``` r
HPDinterval(rvar)
```

    ##          lower     upper
    ## var1 0.4574328 0.6043529
    ## attr(,"Probability")
    ## [1] 0.95

#### **Covariance matrix fixed effects**

``` r
#smfixef
```

#### **Covariance matrix random effects**

``` r
#smranef
```

#### **Repeatability–Log-transformed minimum distance**

``` r
r <- bvar/(bvar + rvar) ###individual
posterior.mode(r)
```

    ##      var1 
    ## 0.3396046

``` r
HPDinterval(r)
```

    ##          lower     upper
    ## var1 0.3117125 0.4177184
    ## attr(,"Probability")
    ## [1] 0.95

#### **Univariate model– Log-transformed truncated minimum distance to the observer (i.e., \>100m is not response)**

``` r
###min distance (truncated to max of 100)----
data$MinDis100<-ifelse(data$MinDisRaw>100,100,data$MinDisRaw)
data$logMinDis100<-log(data$MinDis100+1)
hist(data$logMinDis100)
```

![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/mod2-1.png)<!-- -->

``` r
m1.1<-lmer(logMinDis100~ Sex + NestStage + Year2+ (1|ID)+(1|Observer) ,
         data=data) 
#summary
summary(m1.1)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: logMinDis100 ~ Sex + NestStage + Year2 + (1 | ID) + (1 | Observer)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 986.4
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.7813 -0.6582 -0.0389  0.6082  2.3707 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  ID       (Intercept) 0.33362  0.5776  
    ##  Observer (Intercept) 0.07101  0.2665  
    ##  Residual             0.61580  0.7847  
    ## Number of obs: 369, groups:  ID, 108; Observer, 8
    ## 
    ## Fixed effects:
    ##                        Estimate Std. Error t value
    ## (Intercept)              3.0078     0.1764  17.054
    ## SexMale                 -0.1685     0.1478  -1.140
    ## NestStage2Incubating    -0.3290     0.1245  -2.642
    ## NestStage3Provisioning  -0.5716     0.1349  -4.236
    ## Year21                   0.1734     0.1367   1.268
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) SexMal NstS2I NstS3P
    ## SexMale     -0.400                     
    ## NstStg2Incb -0.405  0.036              
    ## NstStg3Prvs -0.411  0.030  0.604       
    ## Year21      -0.386 -0.022  0.033  0.203

``` r
#model check
plot(m1.1)
```

![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/mod2-2.png)<!-- -->

``` r
hist(resid(m1.1)) # residuals
```

![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/mod2-3.png)<!-- -->

``` r
lattice::qqmath(m1.1) #normality of errors
```

![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/mod2-4.png)<!-- -->

#### **Simulate posterior distribution–Log-transformed truncated minimum distance to the observer (i.e., \>100m is not response)**

``` r
sm1<-sim(m1.1)
smfixef=sm1@fixef
smranef=sm1@ranef
smfixef=as.mcmc(smfixef)
posterior.mode(smfixef)
```

    ##            (Intercept)                SexMale   NestStage2Incubating 
    ##              2.8312227             -0.1116254             -0.3389741 
    ## NestStage3Provisioning                 Year21 
    ##             -0.6188973              0.2587132

``` r
HPDinterval(smfixef)
```

    ##                              lower      upper
    ## (Intercept)             2.77379558  3.3799109
    ## SexMale                -0.41527195  0.1033397
    ## NestStage2Incubating   -0.54896062 -0.1309863
    ## NestStage3Provisioning -0.77900053 -0.3107714
    ## Year21                 -0.05763239  0.4653870
    ## attr(,"Probability")
    ## [1] 0.95

``` r
bID<-sm1@ranef$ID
bvar<-as.vector(apply(bID, 1, var)) ##between individual variance posterior distribution
bvar<-as.mcmc(bvar)
posterior.mode(bvar )## mode of the distribution
```

    ##      var1 
    ## 0.3502112

``` r
HPDinterval(bvar)
```

    ##          lower     upper
    ## var1 0.2710317 0.4180961
    ## attr(,"Probability")
    ## [1] 0.95

``` r
bObs<-sm1@ranef$Observer
ObsVar<-as.vector(apply(bObs, 1, var)) ##between individual variance posterior distribution
ObsVar<-as.mcmc(ObsVar)
posterior.mode(ObsVar )## mode of the distribution
```

    ##       var1 
    ## 0.07088557

``` r
HPDinterval(ObsVar)
```

    ##           lower     upper
    ## var1 0.01767623 0.1139684
    ## attr(,"Probability")
    ## [1] 0.95

``` r
rvar<-sm1@sigma^2
rvar<-as.mcmc(rvar)
posterior.mode(rvar)
```

    ##      var1 
    ## 0.6173605

``` r
HPDinterval(rvar)
```

    ##          lower     upper
    ## var1 0.5369171 0.6895798
    ## attr(,"Probability")
    ## [1] 0.95

#### **Covariance matrix fixed effects**

``` r
#smfixef
```

#### **Covariance matrix random effects**

``` r
#smranef
```

#### **Repeatability–Log-transformed truncated minimum distance to the observer (i.e., \>100m is not response)**

``` r
r <- bvar/(bvar + rvar)        ###individual
posterior.mode(r)
```

    ##      var1 
    ## 0.3416593

``` r
HPDinterval(r)
```

    ##          lower     upper
    ## var1 0.2940487 0.3993569
    ## attr(,"Probability")
    ## [1] 0.95

#### **Univariate model– Number of Stoops**

``` r
###Number of stoops----
#names(data)
m2<-glmer(Dives~ Sex + NestStage + Year2+ (1|ID)+ (1|Observer) ,
         data=data, family="poisson") 
#summary
summary(m2)
```

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: poisson  ( log )
    ## Formula: Dives ~ Sex + NestStage + Year2 + (1 | ID) + (1 | Observer)
    ##    Data: data
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##   2510.4   2537.7  -1248.2   2496.4      362 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.5774 -1.3142 -0.5969  0.8893  6.7682 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  ID       (Intercept) 1.4901   1.2207  
    ##  Observer (Intercept) 0.2758   0.5252  
    ## Number of obs: 369, groups:  ID, 108; Observer, 8
    ## 
    ## Fixed effects:
    ##                        Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)             0.58284    0.27920   2.088 0.036839 *  
    ## SexMale                 0.55899    0.25335   2.206 0.027356 *  
    ## NestStage2Incubating    0.19798    0.08685   2.280 0.022636 *  
    ## NestStage3Provisioning  0.20617    0.09424   2.188 0.028683 *  
    ## Year21                 -0.33174    0.10022  -3.310 0.000932 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) SexMal NstS2I NstS3P
    ## SexMale     -0.464                     
    ## NstStg2Incb -0.177  0.019              
    ## NstStg3Prvs -0.167  0.015  0.734       
    ## Year21      -0.149 -0.005  0.081  0.226

#### **Simulate posterior distribution–Number of Stoops**

``` r
sm2<-sim(m2)
smfixef=sm2@fixef
smranef=sm2@ranef
smfixef=as.mcmc(smfixef)
posterior.mode(smfixef)
```

    ##            (Intercept)                SexMale   NestStage2Incubating 
    ##              0.7273439              0.5824616              0.2293012 
    ## NestStage3Provisioning                 Year21 
    ##              0.2412294             -0.3081508

``` r
HPDinterval(smfixef)
```

    ##                              lower       upper
    ## (Intercept)             0.06645423  1.11216924
    ## SexMale                 0.01334285  0.94270842
    ## NestStage2Incubating    0.01823280  0.36034794
    ## NestStage3Provisioning  0.04599919  0.36537906
    ## Year21                 -0.52311728 -0.06826131
    ## attr(,"Probability")
    ## [1] 0.95

``` r
bID<-sm2@ranef$ID
bvar<-as.vector(apply(bID, 1, var)) ##between individual variance posterior distribution
bvar<-as.mcmc(bvar)
posterior.mode(bvar )## mode of the distribution
```

    ##     var1 
    ## 1.419279

``` r
HPDinterval(bvar)
```

    ##        lower    upper
    ## var1 1.17172 1.587777
    ## attr(,"Probability")
    ## [1] 0.95

``` r
bObs<-sm2@ranef$Observer
ObsVar<-as.vector(apply(bObs, 1, var)) ##between individual variance posterior distribution
ObsVar<-as.mcmc(ObsVar)
posterior.mode(ObsVar )## mode of the distribution
```

    ##      var1 
    ## 0.2908955

``` r
HPDinterval(ObsVar)
```

    ##          lower     upper
    ## var1 0.2050345 0.3739773
    ## attr(,"Probability")
    ## [1] 0.95

#### **Covariance matrix fixed effects**

``` r
#smfixef
```

#### **Covariance matrix random effects**

``` r
#smranef
```

#### **Repeatability– Number of Stoops**

``` r
r <- bvar/(bvar + 1)        ###individual
posterior.mode(r)
```

    ##      var1 
    ## 0.5865412

``` r
HPDinterval(r)
```

    ##          lower     upper
    ## var1 0.5395356 0.6135679
    ## attr(,"Probability")
    ## [1] 0.95

#### **Univariate model– Log-transformed FID**

``` r
###FID----
data3<-subset(data, data$FID>0)
hist(data3$FID)
```

![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/mod4-1.png)<!-- -->

``` r
data3$logFID<-log(data3$FID+1)
hist(data3$logFID)
```

![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/mod4-2.png)<!-- -->

``` r
m3<-lmer(logFID~ Sex + NestStage + Year2+ (1|ID)+(1|Observer) ,
         data=data3) 
#summary
summary(m3)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: logFID ~ Sex + NestStage + Year2 + (1 | ID) + (1 | Observer)
    ##    Data: data3
    ## 
    ## REML criterion at convergence: 938.5
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.9408 -0.5841  0.1320  0.6995  2.2551 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  ID       (Intercept) 0.2042   0.4519  
    ##  Observer (Intercept) 0.0000   0.0000  
    ##  Residual             0.9121   0.9550  
    ## Number of obs: 321, groups:  ID, 105; Observer, 8
    ## 
    ## Fixed effects:
    ##                        Estimate Std. Error t value
    ## (Intercept)             4.38743    0.15729  27.894
    ## SexMale                 0.39660    0.15040   2.637
    ## NestStage2Incubating   -0.47823    0.15153  -3.156
    ## NestStage3Provisioning -0.05798    0.14886  -0.390
    ## Year21                 -0.27550    0.13033  -2.114
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) SexMal NstS2I NstS3P
    ## SexMale     -0.432                     
    ## NstStg2Incb -0.560  0.054              
    ## NstStg3Prvs -0.627  0.073  0.613       
    ## Year21      -0.457 -0.026 -0.030  0.109
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see ?isSingular

``` r
#model check
plot(m3)
```

![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/mod4-3.png)<!-- -->

``` r
hist(resid(m3)) # residuals
```

![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/mod4-4.png)<!-- -->

``` r
lattice::qqmath(m3) #normality of errors
```

![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/mod4-5.png)<!-- -->

#### **Simulate posterior distribution– Log-transformed FID**

``` r
sm3<-sim(m3)
smfixef=sm3@fixef
smranef=sm3@ranef
smfixef=as.mcmc(smfixef)
posterior.mode(smfixef)
```

    ##            (Intercept)                SexMale   NestStage2Incubating 
    ##             4.31996663             0.37158445            -0.46208229 
    ## NestStage3Provisioning                 Year21 
    ##            -0.02500051            -0.30296857

``` r
HPDinterval(smfixef)
```

    ##                              lower       upper
    ## (Intercept)             4.11005173  4.76984195
    ## SexMale                 0.08591674  0.67269418
    ## NestStage2Incubating   -0.76797214 -0.17256696
    ## NestStage3Provisioning -0.39249686  0.27260840
    ## Year21                 -0.53981778 -0.02328108
    ## attr(,"Probability")
    ## [1] 0.95

``` r
bID<-sm3@ranef$ID
bvar<-as.vector(apply(bID, 1, var)) ##between individual variance posterior distribution
bvar<-as.mcmc(bvar)
posterior.mode(bvar )## mode of the distribution
```

    ##      var1 
    ## 0.1848171

``` r
HPDinterval(bvar)
```

    ##          lower     upper
    ## var1 0.1595443 0.2697814
    ## attr(,"Probability")
    ## [1] 0.95

``` r
bObs<-sm3@ranef$Observer
ObsVar<-as.vector(apply(bObs, 1, var)) ##between individual variance posterior distribution
ObsVar<-as.mcmc(ObsVar)
posterior.mode(ObsVar )## mode of the distribution
```

    ##          var1 
    ## -0.0002103502

``` r
HPDinterval(ObsVar)
```

    ##      lower upper
    ## var1     0     0
    ## attr(,"Probability")
    ## [1] 0.95

``` r
rvar<-sm1@sigma^2
rvar<-as.mcmc(rvar)
posterior.mode(rvar)
```

    ##      var1 
    ## 0.6173605

``` r
HPDinterval(rvar)
```

    ##          lower     upper
    ## var1 0.5369171 0.6895798
    ## attr(,"Probability")
    ## [1] 0.95

#### **Covariance matrix fixed effects**

``` r
#smfixef
```

#### **Covariance matrix random effects**

``` r
#smranef
```

#### **Repeatability– Log-transformed FID**

``` r
r <- bvar/(bvar + rvar)        ###individual
posterior.mode(r)
```

    ##      var1 
    ## 0.2587882

``` r
HPDinterval(r)
```

    ##          lower     upper
    ## var1 0.2010614 0.3084137
    ## attr(,"Probability")
    ## [1] 0.95

### **2. Multivariate model with 9 nest defense traits (across contexts)**

-   *Since we had variable effects of nest stage with each behavior
    measured, we constructed three multivariate models to understand the
    correlation of each measure of nest defense with nest stage.*

#### **Multivariate: FID across breeding contexts**

``` r
prior3var=list(R=list(V=diag(3),nu=3.002),G=list(G1=list(V=diag(3), nu=3.002))) #3 trait prior w/one random effect

m3FID<-MCMCglmm(cbind(logLFID, logIFID,logPFID)~(trait-1),
             random=~us(trait):ID,
             rcov=~idh(trait):units,
             family=c("gaussian","gaussian", "gaussian"),
             prior=prior3var,
             nitt=1300000,thin=1000,burnin=300000, 
             data=data_3c,verbose=F)

summary(m3FID)
```

    ## 
    ##  Iterations = 300001:1299001
    ##  Thinning interval  = 1000
    ##  Sample size  = 1000 
    ## 
    ##  DIC: 908.0646 
    ## 
    ##  G-structure:  ~us(trait):ID
    ## 
    ##                              post.mean l-95% CI u-95% CI eff.samp
    ## traitlogLFID:traitlogLFID.ID   0.56260  0.24385   0.8846   1000.0
    ## traitlogIFID:traitlogLFID.ID  -0.01404 -0.31305   0.2414    891.9
    ## traitlogPFID:traitlogLFID.ID  -0.22391 -0.58239   0.0731   1000.0
    ## traitlogLFID:traitlogIFID.ID  -0.01404 -0.31305   0.2414    891.9
    ## traitlogIFID:traitlogIFID.ID   0.47231  0.22446   0.7707   1000.0
    ## traitlogPFID:traitlogIFID.ID   0.17719 -0.05108   0.4471   1000.0
    ## traitlogLFID:traitlogPFID.ID  -0.22391 -0.58239   0.0731   1000.0
    ## traitlogIFID:traitlogPFID.ID   0.17719 -0.05108   0.4471   1000.0
    ## traitlogPFID:traitlogPFID.ID   0.65405  0.29576   1.0413   1000.0
    ## 
    ##  R-structure:  ~idh(trait):units
    ## 
    ##                    post.mean l-95% CI u-95% CI eff.samp
    ## traitlogLFID.units    0.5531   0.3065   0.8274     1158
    ## traitlogIFID.units    0.9743   0.7110   1.2466     1000
    ## traitlogPFID.units    0.7282   0.5291   0.9797     1000
    ## 
    ##  Location effects: cbind(logLFID, logIFID, logPFID) ~ (trait - 1) 
    ## 
    ##              post.mean l-95% CI u-95% CI eff.samp  pMCMC    
    ## traitlogLFID     4.434    4.175    4.703     1336 <0.001 ***
    ## traitlogIFID     3.924    3.653    4.203     1093 <0.001 ***
    ## traitlogPFID     4.281    4.053    4.568     1000 <0.001 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(m3FID)
```

![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/mvFID-1.png)<!-- -->![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/mvFID-2.png)<!-- -->![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/mvFID-3.png)<!-- -->![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/mvFID-4.png)<!-- -->![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/mvFID-5.png)<!-- -->

#### **Covariance matrix fixed effects**

``` r
#m3FID$Sol
```

#### **Covariance matrix random effects**

``` r
#m3FID$VCV
```

##### **Among-individual covariance-FID**

``` r
c1 <- posterior.cor(m3FID$VCV[,1:9])
round(apply(c1,2,mean),2)
```

    ##  var1  var2  var3  var4  var5  var6  var7  var8  var9 
    ##  1.00 -0.02 -0.36 -0.02  1.00  0.32 -0.36  0.32  1.00

``` r
round(apply(c1,2, quantile, c(0.025, 0.975)),2)
```

    ##       var1  var2  var3  var4 var5  var6  var7  var8 var9
    ## 2.5%     1 -0.53 -0.72 -0.53    1 -0.11 -0.72 -0.11    1
    ## 97.5%    1  0.46  0.17  0.46    1  0.65  0.17  0.65    1

``` r
#c1
```

#### **Multivariate: minimum approach distance across breeding contexts**

``` r
m3Min<-MCMCglmm(cbind(logLMinDis, logIMinDis,logPMinDis)~(trait-1),
                random=~us(trait):ID,
                rcov=~idh(trait):units,
                family=c("gaussian","gaussian", "gaussian"),
                prior=prior3var,
                nitt=1300000,thin=1000,burnin=300000, 
                data=data_3c,verbose=FALSE)

summary(m3Min)
```

    ## 
    ##  Iterations = 300001:1299001
    ##  Thinning interval  = 1000
    ##  Sample size  = 1000 
    ## 
    ##  DIC: 844.0413 
    ## 
    ##  G-structure:  ~us(trait):ID
    ## 
    ##                                    post.mean l-95% CI u-95% CI eff.samp
    ## traitlogLMinDis:traitlogLMinDis.ID    0.5488  0.22435   0.9151   1000.0
    ## traitlogIMinDis:traitlogLMinDis.ID    0.1602 -0.08392   0.4031   1000.0
    ## traitlogPMinDis:traitlogLMinDis.ID    0.1422 -0.08378   0.4490   1000.0
    ## traitlogLMinDis:traitlogIMinDis.ID    0.1602 -0.08392   0.4031   1000.0
    ## traitlogIMinDis:traitlogIMinDis.ID    0.5714  0.29596   0.8774   1000.0
    ## traitlogPMinDis:traitlogIMinDis.ID    0.2657  0.03835   0.4928    831.8
    ## traitlogLMinDis:traitlogPMinDis.ID    0.1422 -0.08378   0.4490   1000.0
    ## traitlogIMinDis:traitlogPMinDis.ID    0.2657  0.03835   0.4928    831.8
    ## traitlogPMinDis:traitlogPMinDis.ID    0.6200  0.31976   0.9586   1126.0
    ## 
    ##  R-structure:  ~idh(trait):units
    ## 
    ##                       post.mean l-95% CI u-95% CI eff.samp
    ## traitlogLMinDis.units    0.6763   0.3859   1.0146     1308
    ## traitlogIMinDis.units    0.6090   0.4375   0.8014     1000
    ## traitlogPMinDis.units    0.5777   0.4187   0.7638     1000
    ## 
    ##  Location effects: cbind(logLMinDis, logIMinDis, logPMinDis) ~ (trait - 1) 
    ## 
    ##                 post.mean l-95% CI u-95% CI eff.samp  pMCMC    
    ## traitlogLMinDis     3.002    2.739    3.226   1000.0 <0.001 ***
    ## traitlogIMinDis     2.722    2.485    2.987    865.2 <0.001 ***
    ## traitlogPMinDis     2.412    2.184    2.683   1000.0 <0.001 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(m3Min)
```

![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/mvMin-1.png)<!-- -->![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/mvMin-2.png)<!-- -->![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/mvMin-3.png)<!-- -->![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/mvMin-4.png)<!-- -->![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/mvMin-5.png)<!-- -->

#### **Covariance matrix fixed effects**

``` r
#m3Min$Sol
```

#### **Covariance matrix random effects**

``` r
#m3Min$VCV
```

##### **Among-individual covariance–minimum approach distance across breeding contexts**

``` r
c1 <- posterior.cor(m3Min$VCV[,1:9])
round(apply(c1,2,mean),2)
```

    ## var1 var2 var3 var4 var5 var6 var7 var8 var9 
    ## 1.00 0.29 0.24 0.29 1.00 0.45 0.24 0.45 1.00

``` r
round(apply(c1,2, quantile, c(0.025, 0.975)),2)
```

    ##       var1  var2  var3  var4 var5 var6  var7 var8 var9
    ## 2.5%     1 -0.15 -0.17 -0.15    1 0.11 -0.17 0.11    1
    ## 97.5%    1  0.62  0.62  0.62    1 0.71  0.62 0.71    1

``` r
#c1
```

#### **Multivariate: number of dives across breeding contexts**

``` r
m3Dives<-MCMCglmm(cbind(Ldives, Idives,Pdives)~(trait-1),
                random=~us(trait):ID,
                rcov=~idh(trait):units,
                family=c("poisson","poisson", "poisson"), 
                prior=prior3var,
                nitt=1300000,thin=1000,burnin=300000, 
                data=data_3c,verbose=FALSE)

summary(m3Dives)
```

    ## 
    ##  Iterations = 300001:1299001
    ##  Thinning interval  = 1000
    ##  Sample size  = 1000 
    ## 
    ##  DIC: 1273.432 
    ## 
    ##  G-structure:  ~us(trait):ID
    ## 
    ##                            post.mean l-95% CI u-95% CI eff.samp
    ## traitLdives:traitLdives.ID    1.9246  0.46021    3.517     1110
    ## traitIdives:traitLdives.ID    0.9051  0.08086    1.804     1000
    ## traitPdives:traitLdives.ID    1.6446  0.30928    3.184     1127
    ## traitLdives:traitIdives.ID    0.9051  0.08086    1.804     1000
    ## traitIdives:traitIdives.ID    1.4379  0.68340    2.555     1000
    ## traitPdives:traitIdives.ID    1.6573  0.58315    2.819     1000
    ## traitLdives:traitPdives.ID    1.6446  0.30928    3.184     1127
    ## traitIdives:traitPdives.ID    1.6573  0.58315    2.819     1000
    ## traitPdives:traitPdives.ID    3.3190  1.25416    5.856     1000
    ## 
    ##  R-structure:  ~idh(trait):units
    ## 
    ##                   post.mean l-95% CI u-95% CI eff.samp
    ## traitLdives.units     1.827   0.6231    3.170     1122
    ## traitIdives.units     1.245   0.6447    1.954     1400
    ## traitPdives.units     2.175   1.1986    3.465     1000
    ## 
    ##  Location effects: cbind(Ldives, Idives, Pdives) ~ (trait - 1) 
    ## 
    ##             post.mean l-95% CI u-95% CI eff.samp  pMCMC    
    ## traitLdives   0.39474 -0.18779  0.92822     1000  0.202    
    ## traitIdives   0.78970  0.37382  1.22161     1000 <0.001 ***
    ## traitPdives   0.01581 -0.60554  0.67711     1000  0.926    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(m3Dives)
```

![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/mvStoop-1.png)<!-- -->![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/mvStoop-2.png)<!-- -->![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/mvStoop-3.png)<!-- -->![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/mvStoop-4.png)<!-- -->![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/mvStoop-5.png)<!-- -->

#### **Covariance matrix fixed effects**

``` r
#m3Dives$Sol
```

#### **Covariance matrix random effects**

``` r
#m3Dives$VCV
```

#### **Among-individual covariance–number of dives across breeding contexts**

``` r
c1 <- posterior.cor(m3Dives$VCV[,1:9])
round(apply(c1,2,mean),2)
```

    ## var1 var2 var3 var4 var5 var6 var7 var8 var9 
    ## 1.00 0.55 0.66 0.55 1.00 0.77 0.66 0.77 1.00

``` r
round(apply(c1,2, quantile, c(0.025, 0.975)),2)
```

    ##       var1 var2 var3 var4 var5 var6 var7 var8 var9
    ## 2.5%     1 0.05  0.2 0.05    1 0.50  0.2 0.50    1
    ## 97.5%    1 0.84  0.9 0.84    1 0.91  0.9 0.91    1

``` r
#c1
```

### **3. Bivariate model with minimum distance and number of dives**

-   *Evaluate covariance between each measure of nest defense*

``` r
prior_2var= list(R = list(V = diag(2), nu = 1.002),
                 G = list(G1 = list(V = diag(2), nu = 2,
                                    alpha.mu = rep(0,2),
                                    alpha.V = diag(25^2,2)))) # 3-trait prior w/one random effect
#names(data)

model.2var <- MCMCglmm(
  cbind(LogMin, Dives) ~ trait - 1,
  random = ~ us(trait):ID,
  rcov = ~ us(trait):units,
  family = c("gaussian", "poisson"),
  data = data,
  prior = prior_2var,
  verbose = FALSE,
  nitt = 1300000, thin = 1000, burnin = 300000
  )


summary(model.2var)
```

    ## 
    ##  Iterations = 300001:1299001
    ##  Thinning interval  = 1000
    ##  Sample size  = 1000 
    ## 
    ##  DIC: 2337.386 
    ## 
    ##  G-structure:  ~us(trait):ID
    ## 
    ##                            post.mean l-95% CI u-95% CI eff.samp
    ## traitLogMin:traitLogMin.ID    0.3883   0.2267   0.5595     1000
    ## traitDives:traitLogMin.ID    -0.5300  -0.8180  -0.2607     1096
    ## traitLogMin:traitDives.ID    -0.5300  -0.8180  -0.2607     1096
    ## traitDives:traitDives.ID      1.2814   0.6739   2.0727     1228
    ## 
    ##  R-structure:  ~us(trait):units
    ## 
    ##                               post.mean l-95% CI u-95% CI eff.samp
    ## traitLogMin:traitLogMin.units    0.7252   0.6086   0.8523     1000
    ## traitDives:traitLogMin.units    -0.6280  -0.8275  -0.4542     1000
    ## traitLogMin:traitDives.units    -0.6280  -0.8275  -0.4542     1000
    ## traitDives:traitDives.units      1.9098   1.3499   2.4889     1000
    ## 
    ##  Location effects: cbind(LogMin, Dives) ~ trait - 1 
    ## 
    ##             post.mean l-95% CI u-95% CI eff.samp  pMCMC    
    ## traitLogMin    2.7329   2.5801   2.8992    894.6 <0.001 ***
    ## traitDives     0.4682   0.1827   0.7977   1000.0   0.01 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(model.2var)
```

![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/3bivariate-1.png)<!-- -->![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/3bivariate-2.png)<!-- -->![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/3bivariate-3.png)<!-- -->

#### **Covariance matrix fixed effects**

``` r
#model.2var$Sol
```

#### **Covariance matrix random effects**

``` r
#model.2var$VCV
```

#### **Among-individual and within-individual covariance–minimum distance and number of dives**

``` r
c1<- posterior.cor(model.2var$VCV[,1:4])
round(apply(c1,2,mean),2)
```

    ##  var1  var2  var3  var4 
    ##  1.00 -0.75 -0.75  1.00

``` r
round(apply(c1,2, quantile, c(0.025, 0.975)),2)
```

    ##       var1  var2  var3 var4
    ## 2.5%     1 -0.90 -0.90    1
    ## 97.5%    1 -0.54 -0.54    1

``` r
#c1
#within individual covariance
c2<- posterior.cor(model.2var$VCV[,5:8])
round(apply(c2,2,mean),2)
```

    ##  var1  var2  var3  var4 
    ##  1.00 -0.53 -0.53  1.00

``` r
round(apply(c2,2, quantile, c(0.025, 0.975)),2)
```

    ##       var1  var2  var3 var4
    ## 2.5%     1 -0.63 -0.63    1
    ## 97.5%    1 -0.43 -0.43    1

``` r
#c2
```

### **4. Main effects–Univariate model to estimate short-versus long-term repeatability with interaction**

-   *Interaction is included since we predicted that females would be
    more defensive with increasing nest stages.*

``` r
cor.test(data$Stage, data$NumberOfChicks)#correlation b/n nest stage and number of chicks. r=0.78 thus drop Number of Chicks from model 
```

    ## 
    ##  Pearson's product-moment correlation
    ## 
    ## data:  data$Stage and data$NumberOfChicks
    ## t = 21.628, df = 367, p-value < 2.2e-16
    ## alternative hypothesis: true correlation is not equal to 0
    ## 95 percent confidence interval:
    ##  0.6999745 0.7902693
    ## sample estimates:
    ##       cor 
    ## 0.7485722

``` r
m1<-lmer(LogMin~ (Sex+NestStage)^2 + Year2+ (1|ID) + (1|ID_Series) ,
         data=data) 
#summary
summary(m1)
```

    ## Linear mixed model fit by REML ['lmerMod']
    ## Formula: LogMin ~ (Sex + NestStage)^2 + Year2 + (1 | ID) + (1 | ID_Series)
    ##    Data: data
    ## 
    ## REML criterion at convergence: 1003.1
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.94297 -0.61202 -0.02309  0.52714  2.76760 
    ## 
    ## Random effects:
    ##  Groups    Name        Variance Std.Dev.
    ##  ID_Series (Intercept) 0.08504  0.2916  
    ##  ID        (Intercept) 0.27454  0.5240  
    ##  Residual              0.65003  0.8062  
    ## Number of obs: 369, groups:  ID_Series, 123; ID, 108
    ## 
    ## Fixed effects:
    ##                                Estimate Std. Error t value
    ## (Intercept)                     2.94076    0.16208  18.144
    ## SexMale                        -0.09501    0.21235  -0.447
    ## NestStage2Incubating           -0.28349    0.16770  -1.690
    ## NestStage3Provisioning         -0.55002    0.15752  -3.492
    ## Year21                          0.29835    0.12778   2.335
    ## SexMale:NestStage2Incubating   -0.10133    0.25100  -0.404
    ## SexMale:NestStage3Provisioning -0.13211    0.24497  -0.539
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) SexMal NstS2I NstS3P Year21 SM:NS2
    ## SexMale     -0.628                                   
    ## NstStg2Incb -0.586  0.457                            
    ## NstStg3Prvs -0.602  0.444  0.647                     
    ## Year21      -0.411 -0.016 -0.029  0.045              
    ## SxMl:NstS2I  0.396 -0.642 -0.668 -0.433  0.008       
    ## SxMl:NstS3P  0.369 -0.615 -0.417 -0.641  0.013  0.604

``` r
#model check
plot(m1)
```

![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/unvME-1.png)<!-- -->

``` r
hist(resid(m1)) # residuals
```

![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/unvME-2.png)<!-- -->

``` r
lattice::qqmath(m1) #normality of errors
```

![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/unvME-3.png)<!-- -->

##### **Simulate posterior distribution– Univariate model with no interaction**

``` r
#simulated parameters
require(dplyr)
library(arm)
require(MCMCglmm)
require(ggplot2)
smod<-sim(m1,1000)
describe_posterior(smod, ci = 0.95)
```

    ## Summary of Posterior Distribution
    ## 
    ## Parameter                      | Median |         95% CI |     pd |          ROPE | % in ROPE
    ## ---------------------------------------------------------------------------------------------
    ## (Intercept)                    |   2.95 | [ 2.61,  3.26] |   100% | [-0.10, 0.10] |        0%
    ## SexMale                        |  -0.10 | [-0.51,  0.29] | 66.40% | [-0.10, 0.10] |    33.16%
    ## NestStage2Incubating           |  -0.29 | [-0.61,  0.04] | 96.10% | [-0.10, 0.10] |    10.95%
    ## NestStage3Provisioning         |  -0.55 | [-0.86, -0.27] |   100% | [-0.10, 0.10] |        0%
    ## Year21                         |   0.30 | [ 0.05,  0.55] | 99.20% | [-0.10, 0.10] |     3.16%
    ## SexMale:NestStage2Incubating   |  -0.08 | [-0.58,  0.37] | 63.00% | [-0.10, 0.10] |    31.79%
    ## SexMale:NestStage3Provisioning |  -0.13 | [-0.61,  0.36] | 69.00% | [-0.10, 0.10] |    27.47%

``` r
posterior.mode(as.mcmc(smod@fixef))
```

    ##                    (Intercept)                        SexMale 
    ##                      2.9084015                     -0.1823287 
    ##           NestStage2Incubating         NestStage3Provisioning 
    ##                     -0.3199375                     -0.5319123 
    ##                         Year21   SexMale:NestStage2Incubating 
    ##                      0.2978666                     -0.2621571 
    ## SexMale:NestStage3Provisioning 
    ##                     -0.1425156

``` r
HPDinterval(as.mcmc(smod@fixef))
```

    ##                                      lower       upper
    ## (Intercept)                     2.61773030  3.26475776
    ## SexMale                        -0.46964788  0.31599217
    ## NestStage2Incubating           -0.62189654  0.02599652
    ## NestStage3Provisioning         -0.85215405 -0.26034640
    ## Year21                          0.05897453  0.54787348
    ## SexMale:NestStage2Incubating   -0.51007121  0.41229142
    ## SexMale:NestStage3Provisioning -0.58685428  0.38060616
    ## attr(,"Probability")
    ## [1] 0.95

``` r
##Between individual variance
bID<-smod@ranef$ID[,,1]
bvar<-as.vector(apply(bID, 1, var)) ##ID variance posterior distribution

require(MCMCglmm)
bvar4<-as.mcmc(bvar)
posterior.mode(bvar4)## mode of the distribution
```

    ##      var1 
    ## 0.2927175

``` r
HPDinterval(bvar4)
```

    ##          lower    upper
    ## var1 0.2404513 0.410325
    ## attr(,"Probability")
    ## [1] 0.95

``` r
require(bayestestR)
bvar11<-as.data.frame(bvar) #posterior plot
describe_posterior(bvar11)
```

    ## Summary of Posterior Distribution
    ## 
    ## Parameter | Median |       95% CI |   pd |          ROPE | % in ROPE
    ## --------------------------------------------------------------------
    ## bvar      |   0.31 | [0.24, 0.41] | 100% | [-0.10, 0.10] |        0%

``` r
bv<-describe_posterior(bvar11)#data frame for CrI and median

ggplot(bvar11, aes(x = bvar)) +
  geom_density(fill = "orange") +
  geom_vline(xintercept = median(bv$Median), color="red", size=1)+
   geom_vline(xintercept = (bv$CI_low), color="black", size=1, linetype="longdash")+
  geom_vline(xintercept = (bv$CI_high), color="black", size=1, linetype="longdash")+
  labs(y="Density", x="Between-Individual Variance (ID)")+
  ggtitle(label ="Density Plot Between-Indvidual Variance (ID)",
          subtitle = "Posterior distribution plotted")+
  theme(plot.title = element_text(face = "bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5))
```

![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/unMESim-1.png)<!-- -->

``` r
##Between individual variance, ID_Series
bID2<-smod@ranef$ID_Series[,,1]
bvar2<-as.vector(apply(bID2, 1, var)) ##ID_Series variance posterior distribution
require(MCMCglmm)
bvar22<-as.mcmc(bvar2)
posterior.mode(bvar22)## mode of the distribution
```

    ##      var1 
    ## 0.0728038

``` r
HPDinterval(bvar22)
```

    ##           lower      upper
    ## var1 0.05475697 0.09562029
    ## attr(,"Probability")
    ## [1] 0.95

``` r
bvar112<-as.data.frame(bvar2) #posterior plot
describe_posterior(bvar112)
```

    ## Summary of Posterior Distribution
    ## 
    ## Parameter | Median |       95% CI |   pd |          ROPE | % in ROPE
    ## --------------------------------------------------------------------
    ## bvar2     |   0.07 | [0.06, 0.10] | 100% | [-0.10, 0.10] |      100%

``` r
bv2<-describe_posterior(bvar112)#data frame for CrI and median

ggplot(bvar112, aes(x = bvar2)) +
  geom_density(fill = "orange") +
  geom_vline(xintercept = median(bv2$Median), color="red", size=1)+
   geom_vline(xintercept = (bv2$CI_low), color="black", size=1, linetype="longdash")+
  geom_vline(xintercept = (bv2$CI_high), color="black", size=1, linetype="longdash")+
  labs(y="Density", x="Between-Individual Variance (ID_Series)")+
  ggtitle(label ="Density Plot Between-Indvidual Variance (ID_Series)",
          subtitle = "Posterior distribution plotted")+
  theme(plot.title = element_text(face = "bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5))
```

![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/unMESim-2.png)<!-- -->

``` r
###residual variance
rvar<-smod@sigma^2
rvar<-as.mcmc(rvar)
posterior.mode(rvar)
```

    ##      var1 
    ## 0.6587313

``` r
HPDinterval(rvar)
```

    ##          lower     upper
    ## var1 0.5593159 0.7500881
    ## attr(,"Probability")
    ## [1] 0.95

``` r
rvarr<-smod@sigma^2
describe_posterior(rvarr)
```

    ## Summary of Posterior Distribution
    ## 
    ## Parameter | Median |       95% CI |   pd |          ROPE | % in ROPE
    ## --------------------------------------------------------------------
    ## Posterior |   0.65 | [0.56, 0.76] | 100% | [-0.10, 0.10] |        0%

``` r
rv<-describe_posterior(rvarr)#data frame for CrI and median
rvar11<-as.data.frame(rvarr) #posterior plot
ggplot(rvar11, aes(x = rvarr)) +
  geom_density(fill = "orange") +
  geom_vline(xintercept = median(rv$Median), color="red", size=1)+
   geom_vline(xintercept = (rv$CI_low), color="black", size=1, linetype="longdash")+
  geom_vline(xintercept = (rv$CI_high), color="black", size=1, linetype="longdash")+
  labs(y="Density", x="Residual Variance")+
  ggtitle(label ="Density Plot Residual Variance",
          subtitle = "Posterior distribution plotted")+
  theme(plot.title = element_text(face = "bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5))
```

![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/unMESim-3.png)<!-- -->

#### **Covariance matrix fixed effects**

``` r
#smod@fixef
```

#### **Covariance matrix random effects**

``` r
#smod@ranef
```

#### **Contrasts differences among-sexes**

``` r
#compute difference among-sexes

p222<-smod@fixef[,1:2]#gather posterior for intercept(female) and Male 
p222<-as.data.frame(p222)#make a dataframe to manipulate easier
colnames(p222)<-c("Female", "SexMale")#change column names
p22<-p222%>%
  summarise(Male= Female+ SexMale)#add posteriors for Intercept to each value of Male

con<-as.data.frame(c(p222, p22))#combine estimates into one dataframe

difference<-con%>%
  summarise(diff= Female - Male)#compute difference

contrast<-as.data.frame(c(con, difference))#dataframe for contrasts

require(bayestestR)
describe_posterior(contrast)# pd is 87% for difference among sexes 
```

    ## Summary of Posterior Distribution
    ## 
    ## Parameter | Median |        95% CI |     pd |          ROPE | % in ROPE
    ## -----------------------------------------------------------------------
    ## Female    |   2.95 | [ 2.61, 3.26] |   100% | [-0.10, 0.10] |        0%
    ## SexMale   |  -0.10 | [-0.51, 0.29] | 66.40% | [-0.10, 0.10] |    33.16%
    ## Male      |   2.85 | [ 2.52, 3.18] |   100% | [-0.10, 0.10] |        0%
    ## diff      |   0.10 | [-0.29, 0.51] | 66.40% | [-0.10, 0.10] |    33.16%

``` r
contable<-describe_posterior(contrast)#estimates for difference b/n male and female nest defense

contab<-contable%>%filter(Parameter=="diff")



ggplot(contrast, aes(x=diff)) + geom_density(fill="orange") +
  geom_vline(xintercept = median(contrast$diff), color="red", size=1)+
   geom_vline(xintercept = (contab$CI_low), color="black", size=1, linetype="longdash")+
  geom_vline(xintercept = (contab$CI_high), color="black", size=1, linetype="longdash")+
  labs(x= "Difference", y="Density")+
  ggtitle(label ="Difference between male and female defense",
          subtitle = "Posterior Distribution, Median, and Highest Density Interval")+
  theme(plot.title = element_text(face = "bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5))
```

![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/contrast-1.png)<!-- -->

#### **Contrast differences among-nest stages**

``` r
p2223<-smod@fixef[,3:4]#gather posterior for incubation and provisioning nest stages 
p2223<-as.data.frame(p2223)#make a dataframe to manipulate easier
p2224<-smod@fixef[,1]#extract intercept (estimate for egg-laying)
p2224<-as.data.frame(p2224)
colnames(p2224)<-"Intercept"

nestcon<-as.data.frame(c(p2224, p2223))

#compute estimates for each nest stage
neststage<-nestcon%>%
  summarise(Incubation   = Intercept  + NestStage2Incubating,
            Provisioning = Intercept + NestStage3Provisioning,
            EggLaying    = Intercept)#add posteriors for Intercept to each value of Male

#compute differences between nest stage groups
contrasts<-neststage%>%
  summarise(diff_ie = Incubation-EggLaying,
            diff_ip = Incubation-Provisioning,
            diff_ep = EggLaying-Provisioning)
require(bayestestR)
describe_posterior(contrasts)#median, 95% CrI and PD for each difference
```

    ## Summary of Posterior Distribution
    ## 
    ## Parameter | Median |        95% CI |     pd |          ROPE | % in ROPE
    ## -----------------------------------------------------------------------
    ## diff_ie   |  -0.29 | [-0.61, 0.04] | 96.10% | [-0.10, 0.10] |    10.95%
    ## diff_ip   |   0.27 | [ 0.00, 0.53] | 97.50% | [-0.10, 0.10] |     9.16%
    ## diff_ep   |   0.55 | [ 0.27, 0.86] |   100% | [-0.10, 0.10] |        0%

``` r
contable<-describe_posterior(contrasts)#estimates for difference b/n male and female nest defense
diff_ie<-contable%>%filter(Parameter=="diff_ie")
diff_ip<-contable%>%filter(Parameter=="diff_ip")
diff_ep<-contable%>%filter(Parameter=="diff_ep")

# incubation vs egg-laying
p3<-ggplot(contrasts, aes(x=diff_ie)) + geom_density(fill="orange") +
  geom_vline(xintercept = median(contrasts$diff_ie), color="red", size=1)+
   geom_vline(xintercept = (diff_ie$CI_low), color="black", size=1, linetype="longdash")+
  geom_vline(xintercept = (diff_ie$CI_high), color="black", size=1, linetype="longdash")+
  labs(x= "Difference", y="Density")+
  ggtitle(label ="Difference between incubation and egg laying",
          subtitle = "Posterior Distribution, Median, and Highest Density Interval")+
  theme(plot.title = element_text(face = "bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5))
p3
```

![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/contrast2-1.png)<!-- -->

``` r
# incubation vs provisioning
p2<-ggplot(contrasts, aes(x=diff_ip)) + geom_density(fill="orange") +
  geom_vline(xintercept = median(contrasts$diff_ip), color="red", size=1)+
   geom_vline(xintercept = (diff_ip$CI_low), color="black", size=1, linetype="longdash")+
  geom_vline(xintercept = (diff_ip$CI_high), color="black", size=1, linetype="longdash")+
  labs(x= "Difference", y="Density")+
  ggtitle(label ="Difference between incubation and provisioning",
          subtitle = "Posterior Distribution, Median, and Highest Density Interval")+
  theme(plot.title = element_text(face = "bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5))
p2
```

![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/contrast2-2.png)<!-- -->

``` r
# egg laying vs provisioning
p1<-ggplot(contrasts, aes(x=diff_ep)) + geom_density(fill="orange") +
  geom_vline(xintercept = median(contrasts$diff_ep), color="red", size=1)+
   geom_vline(xintercept = (diff_ep$CI_low), color="black", size=1, linetype="longdash")+
  geom_vline(xintercept = (diff_ep$CI_high), color="black", size=1, linetype="longdash")+
  labs(x= "Difference", y="Density")+
  ggtitle(label ="Difference between incubation and provisioning",
          subtitle = "Posterior Distribution, Median, and Highest Density Interval")+
  theme(plot.title = element_text(face = "bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5))
p1
```

![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/contrast2-3.png)<!-- -->

### **Contrast difference among years**

``` r
p55<-smod@fixef[,5]#gather posterior for year
p55<-as.data.frame(p222)#make a dataframe to manipulate easier
colnames(p55)<-c("Year2019")#change column names
p2224<-smod@fixef[,1]#extract intercep
p2224<-as.data.frame(p2224)
colnames(p2224)<-"Intercept"

yearcon<-as.data.frame(c(p2224, p55))

p2299<-yearcon%>%
  summarise(Y2019 = Intercept+ Year2019)#add posteriors for Intercept(2018) to each value                                         #of Male

ycon<-as.data.frame(c(yearcon, p2299)) 


difference_y<-ycon%>%
  summarise(diff_year= Intercept - Y2019)#compute difference


require(bayestestR)
describe_posterior(difference_y)# pd is 87% for difference among sexes 
```

    ## Summary of Posterior Distribution
    ## 
    ## Parameter | Median |         95% CI |   pd |          ROPE | % in ROPE
    ## ----------------------------------------------------------------------
    ## diff_year |  -2.95 | [-3.26, -2.61] | 100% | [-0.10, 0.10] |        0%

``` r
contable<-describe_posterior(difference_y)#estimates for difference b/n male and female nest defense




ggplot(difference_y, aes(x=diff_year)) + geom_density(fill="orange") +
  geom_vline(xintercept = median(difference_y$diff_year), color="red", size=1)+
   geom_vline(xintercept = (contable$CI_low), color="black", size=1, linetype="longdash")+
  geom_vline(xintercept = (contable$CI_high), color="black", size=1, linetype="longdash")+
  labs(x= "Difference", y="Density")+
  ggtitle(label ="Difference between 2018 and 2019",
          subtitle = "Posterior Distribution, Median, and Highest Density Interval")+
  theme(plot.title = element_text(face = "bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5))
```

![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/contrast5-1.png)<!-- -->

##### **Long-term and short-term repeatability–Univariate model with no interaction**

``` r
### Long-term repeatability
r<-bvar4/(rvar+bvar4 +bvar22)
r<-as.mcmc(r)
posterior.mode(r)
```

    ##      var1 
    ## 0.2990986

``` r
HPDinterval(r)  ##repeatability 
```

    ##          lower    upper
    ## var1 0.2506514 0.359143
    ## attr(,"Probability")
    ## [1] 0.95

``` r
rpt<-bvar4/(rvar+bvar4 +bvar22)
rpt<-as.data.frame(rpt)
rptt<-describe_posterior(rpt)#data frame for CrI and median
rpt11<-as.data.frame(rpt) #posterior plot
ggplot(rpt11, aes(x = var1)) +
  geom_density(fill = "orange") +
  geom_vline(xintercept = median(rptt$Median), color="red", size=1)+
   geom_vline(xintercept = (rptt$CI_low), color="black", size=1, linetype="longdash")+
  geom_vline(xintercept = (rptt$CI_high), color="black", size=1, linetype="longdash")+
  labs(y="Density", x="Long-term Repeatability")+
  ggtitle(label ="Density Plot Long-term Repeatability",
          subtitle = "Posterior distribution plotted")+
  theme(plot.title = element_text(face = "bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5))
```

![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/lvsrpt-1.png)<!-- -->

``` r
###Short-term
r1<-(bvar4 +bvar22  )/(rvar+bvar4 +bvar22)
r1<-as.mcmc(r1)
posterior.mode(r1)
```

    ##      var1 
    ## 0.3715559

``` r
HPDinterval(r1)  ##repeatability
```

    ##         lower     upper
    ## var1 0.323377 0.4224197
    ## attr(,"Probability")
    ## [1] 0.95

``` r
rptr<-(bvar4 +bvar22  )/(rvar+bvar4 +bvar22)
rptr<-as.data.frame(rptr)
describe_posterior(rptr)
```

    ## Summary of Posterior Distribution
    ## 
    ## Parameter | Median |       95% CI |   pd |          ROPE | % in ROPE
    ## --------------------------------------------------------------------
    ## var1      |   0.37 | [0.32, 0.42] | 100% | [-0.10, 0.10] |        0%

``` r
rpttt<-describe_posterior(rptr)#data frame for CrI and median
rpt111<-as.data.frame(rptr) #posterior plot
ggplot(rpt111, aes(x = var1)) +
  geom_density(fill = "orange") +
  geom_vline(xintercept = median(rpttt$Median), color="red", size=1)+
   geom_vline(xintercept = (rpttt$CI_low), color="black", size=1, linetype="longdash")+
  geom_vline(xintercept = (rpttt$CI_high), color="black", size=1, linetype="longdash")+
  labs(y="Density", x="Short-term Repeatability")+
  ggtitle(label ="Density Plot Short-term Repeatability",
          subtitle = "Posterior distribution plotted")+
  theme(plot.title = element_text(face = "bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5))
```

![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/lvsrpt-2.png)<!-- -->

##### **Plot–Univariate model with no interaction**

``` r
#load packages for predictive plot
require(ggplot2)
require(ggeffects)
require(cowplot)

d18<-subset(data, data$Year=="2018")
d19<-subset(data, data$Year=="2019")

#plot
p18<-ggplot(d18, aes(x = NestStage, y = MinDisRaw, factor=NestStage, fill=Sex)) + geom_boxplot() +
  scale_y_continuous("Minimum approach distnace to observer (m)",limits = c(0, 100)) +
  scale_x_discrete("Nest Stage", labels=c('Egg-laying', 'Incubating', 'Provisioning'))+ 
  guides(fill="none")+
  guides(color="none")+ theme_classic(base_size = 12)  +
  ggtitle(label ="2018",
          subtitle = "")+
  theme(plot.title = element_text(face = "bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5))



p19<-ggplot(d19, aes(x = NestStage, y = MinDisRaw, factor=NestStage, fill=Sex)) + geom_boxplot() +
  scale_y_continuous("Minimum approach distnace to observer (m)",limits = c(0, 100)) +
  scale_x_discrete("Nest Stage", labels=c('Egg-laying', 'Incubating', 'Provisioning'))+ 
  guides(color="none")+ theme_classic(base_size = 12)  +
  ggtitle(label ="2019",
          subtitle = "")+
  theme(plot.title = element_text(face = "bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5))

plot_grid(p18,p19)
```

![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/plot1-1.png)<!-- -->

``` r
#posterior of difference plots for groups in main effects (supp material)
cowplot::plot_grid(p3, p2, p1, nrow = 2,labels = c("A.", "B.", "C."))
```

![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/plot1-2.png)<!-- -->

#### **Density plots of behavioral types**

``` r
## density plots

post1<-smod@ranef$ID
post1<-as.data.frame(post1)

colnames(post1)<-(c("01H", "03M", "119Female", 
                         "119Male", "11A", "12Female",
                         "12Male", "145Male","14H","155Female",
                         "155Male", "156Female", "156Male", 
                         "157Female", "15C ", "160Male",
                         "168Female", "168Male","175Male",
                         "179Female", "179Male","184Male", 
                         "186Female", "18K", "193Male", 
                         "196Male", "19K", "1Male",
                         "201Female", "20H", "214Male", 
                         "215Female", "215Male", "217Female",       
                         "217Male", "219Female","219Male",     
                         "223Female", "223Male", "227Male", 
                         "228Female", "228Male", "231Female",       
                         "234Male", "236Female", "23Female",
                         "23Male",  "241Female", "241Male",
                         "244Female", "244Male", "24K",     
                         "251Female", "253Male", "254Female" ,      
                         "254Male", "256Female", "25K",
                         "304Female", "304Male", "30Female",
                         "30Male", "318Female", "318Male", 
                         "319Male", "32K", "38A",       
                         "39E", "40Male", "41H",
                         "43A", "43H", "44H",
                         "45K", "46K", "47E", "47Female",
                         "47Male","49E", "49H", 
                         "4Female", "4Male", "50H",
                         "53Female", "53Male", "55K",       
                         "57Male", "58Male" , "59Male",
                         "75B", "78Male", "80Female", "83K",
                         "85C", "86C", "87C",
                         "87K", "88C",  "89C",
                         "89E", "90C", "91C",
                         "92H", "93E", "97Female",
                         "97Male", "98Female", "98Male"))


#convert data to long format for plotting
require(tidyr)
plot_data<-post1%>% 
  pivot_longer(c("01H", "03M", "119Female", 
                         "119Male", "11A", "12Female",
                         "12Male", "145Male","14H","155Female",
                         "155Male", "156Female", "156Male", 
                         "157Female", "15C ", "160Male",
                         "168Female", "168Male","175Male",
                         "179Female", "179Male","184Male", 
                         "186Female", "18K", "193Male", 
                         "196Male", "19K", "1Male",
                         "201Female", "20H", "214Male", 
                         "215Female", "215Male", "217Female",       
                         "217Male", "219Female","219Male",     
                         "223Female", "223Male", "227Male", 
                         "228Female", "228Male", "231Female",       
                         "234Male", "236Female", "23Female",
                         "23Male",  "241Female", "241Male",
                         "244Female", "244Male", "24K",     
                         "251Female", "253Male", "254Female" ,      
                         "254Male", "256Female", "25K",
                         "304Female", "304Male", "30Female",
                         "30Male", "318Female", "318Male", 
                         "319Male", "32K", "38A",       
                         "39E", "40Male", "41H",
                         "43A", "43H", "44H",
                         "45K", "46K", "47E", "47Female",
                         "47Male","49E", "49H", 
                         "4Female", "4Male", "50H",
                         "53Female", "53Male", "55K",       
                         "57Male", "58Male" , "59Male",
                         "75B", "78Male", "80Female", "83K",
                         "85C", "86C", "87C",
                         "87K", "88C",  "89C",
                         "89E", "90C", "91C",
                         "92H", "93E", "97Female",
                         "97Male", "98Female", "98Male"), names_to="ID", values_to="value")


# add population level mean residual variance
plot_data$value <-
plot_data$value +
fixef(m1, pars = "Intercept")[1]




BT<-plot_data%>%
  dplyr:: group_by(ID)%>%
  dplyr:: mutate(meanBT=mean(value))%>%
  dplyr:: ungroup()


data_r<-data%>%
  dplyr::select(ID,Sex)

ref<-data_r[!duplicated(data_r),]


bTT<-merge(BT, ref, all.x = T, no.dups = T)

BTT<-bTT%>%replace_na(list(Sex='Male'))



bt<-ggplot()+
  ggridges::geom_density_ridges(data=BTT,
                                aes(x=value,
                                    y=reorder(as.factor(ID), meanBT),
                                    height=..density..,
                                    fill=ID, scale=3), alpha=0.6) +
  geom_point(data = BTT[!duplicated(BTT$ID),],
             aes(x = meanBT,
                 y = as.factor(ID),
                 col = Sex),
             size = 1.8)+
  labs(y="Tag ID", x="Behavioral Types (Log-transformed Nest Defense)")+
  ggtitle(label ="Among-individual differences in nest defense",
          subtitle = "Posterior distribution for each individual is plotted")+
  theme(plot.title = element_text(face = "bold", hjust=0.5), plot.subtitle=element_text(hjust=0.7), 
        text = element_text(size=18),axis.ticks.y=element_blank(),axis.text.y=element_blank()) #remove y axis ticks)
bt + theme(aspect.ratio=12/5)#+guides(fill="none")+ scale_fill_discrete(guide="none")+ theme(legend.position="none")
```

![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/denplot-1.png)<!-- -->

### **5. Models of (dis-)assortative mating + relative fitness**

-   *Multivariate models that estimate assortative mating (among-pair
    correlation), response to labile environment (within-pair
    correlation), and selection gradients for each combination of sex
    and year.*

#### **2018 3 trait model**

``` r
mod.12 <- MCMCglmm(cbind(Male, Female, rfit) ~ (trait-1),  
                   random = ~us(trait):NestID ,
                   rcov = ~us(trait):units, 
                   family = c("gaussian", "gaussian", "gaussian"),
                   data=data2018, 
                   prior = prior_E_B_fit_1px, 
                   verbose = FALSE,
                   nitt = 1300000, thin = 1000, burnin = 300000
                   )
summary(mod.12)
```

    ## 
    ##  Iterations = 300001:1299001
    ##  Thinning interval  = 1000
    ##  Sample size  = 1000 
    ## 
    ##  DIC: 117.6162 
    ## 
    ##  G-structure:  ~us(trait):NestID
    ## 
    ##                                post.mean   l-95% CI u-95% CI eff.samp
    ## traitMale:traitMale.NestID       0.58115  7.899e-02   1.1814   1000.0
    ## traitFemale:traitMale.NestID     0.04334 -1.002e-01   0.2062   1000.0
    ## traitrfit:traitMale.NestID       0.02227 -3.646e-01   0.3898   1000.0
    ## traitMale:traitFemale.NestID     0.04334 -1.002e-01   0.2062   1000.0
    ## traitFemale:traitFemale.NestID   0.09135  1.679e-06   0.2297   1000.0
    ## traitrfit:traitFemale.NestID     0.10413 -4.976e-02   0.3479   1000.0
    ## traitMale:traitrfit.NestID       0.02227 -3.646e-01   0.3898   1000.0
    ## traitFemale:traitrfit.NestID     0.10413 -4.976e-02   0.3479   1000.0
    ## traitrfit:traitrfit.NestID       0.89323  4.267e-01   1.4749    274.4
    ## 
    ##  R-structure:  ~us(trait):units
    ## 
    ##                                post.mean  l-95% CI u-95% CI eff.samp
    ## traitMale:traitMale.units      6.333e-01  0.382869 0.900644     1000
    ## traitFemale:traitMale.units    4.346e-02 -0.091399 0.190046     1000
    ## traitrfit:traitMale.units      2.172e-04 -0.005862 0.006547     1000
    ## traitMale:traitFemale.units    4.346e-02 -0.091399 0.190046     1000
    ## traitFemale:traitFemale.units  3.814e-01  0.235881 0.516824     1000
    ## traitrfit:traitFemale.units   -2.598e-05 -0.004621 0.004024     1000
    ## traitMale:traitrfit.units      2.172e-04 -0.005862 0.006547     1000
    ## traitFemale:traitrfit.units   -2.598e-05 -0.004621 0.004024     1000
    ## traitrfit:traitrfit.units      1.000e-04  0.000100 0.000100        0
    ## 
    ##  Location effects: cbind(Male, Female, rfit) ~ (trait - 1) 
    ## 
    ##             post.mean l-95% CI u-95% CI eff.samp  pMCMC    
    ## traitMale     -2.3734  -2.6958  -2.0158     1000 <0.001 ***
    ## traitFemale   -2.5874  -2.7833  -2.3916     1000 <0.001 ***
    ## traitrfit      0.7933   0.4570   1.2014     1000 <0.001 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(mod.12)
```

![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/2018trait-1.png)<!-- -->![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/2018trait-2.png)<!-- -->![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/2018trait-3.png)<!-- -->![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/2018trait-4.png)<!-- -->![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/2018trait-5.png)<!-- -->![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/2018trait-6.png)<!-- -->

``` r
#auto-correlation
autocorr.diag(mod.12$VCV)
```

    ##           traitMale:traitMale.NestID traitFemale:traitMale.NestID
    ## Lag 0                     1.00000000                 1.0000000000
    ## Lag 1000                 -0.03474776                 0.0002459054
    ## Lag 5000                  0.06280472                 0.0517662517
    ## Lag 10000                 0.03713087                -0.0005260866
    ## Lag 50000                -0.02500029                 0.0039811555
    ##           traitrfit:traitMale.NestID traitMale:traitFemale.NestID
    ## Lag 0                     1.00000000                 1.0000000000
    ## Lag 1000                 -0.03562206                 0.0002459054
    ## Lag 5000                 -0.03596761                 0.0517662517
    ## Lag 10000                 0.01581216                -0.0005260866
    ## Lag 50000                -0.04594951                 0.0039811555
    ##           traitFemale:traitFemale.NestID traitrfit:traitFemale.NestID
    ## Lag 0                        1.000000000                  1.000000000
    ## Lag 1000                    -0.003170393                 -0.027769944
    ## Lag 5000                    -0.018145788                  0.002070384
    ## Lag 10000                    0.011986780                  0.008254789
    ## Lag 50000                   -0.069373584                  0.013591875
    ##           traitMale:traitrfit.NestID traitFemale:traitrfit.NestID
    ## Lag 0                     1.00000000                  1.000000000
    ## Lag 1000                 -0.03562206                 -0.027769944
    ## Lag 5000                 -0.03596761                  0.002070384
    ## Lag 10000                 0.01581216                  0.008254789
    ## Lag 50000                -0.04594951                  0.013591875
    ##           traitrfit:traitrfit.NestID traitMale:traitMale.units
    ## Lag 0                     1.00000000                1.00000000
    ## Lag 1000                  0.05059957               -0.01314756
    ## Lag 5000                  0.05487032                0.02901474
    ## Lag 10000                 0.02652232               -0.02664002
    ## Lag 50000                -0.03519444               -0.03771583
    ##           traitFemale:traitMale.units traitrfit:traitMale.units
    ## Lag 0                     1.000000000              1.0000000000
    ## Lag 1000                 -0.011922038              0.0217897344
    ## Lag 5000                  0.014053275              0.0172967372
    ## Lag 10000                 0.019021634             -0.0163914620
    ## Lag 50000                 0.006902238              0.0007934218
    ##           traitMale:traitFemale.units traitFemale:traitFemale.units
    ## Lag 0                     1.000000000                    1.00000000
    ## Lag 1000                 -0.011922038                   -0.02947878
    ## Lag 5000                  0.014053275                    0.05677563
    ## Lag 10000                 0.019021634                    0.01963331
    ## Lag 50000                 0.006902238                   -0.01650059
    ##           traitrfit:traitFemale.units traitMale:traitrfit.units
    ## Lag 0                     1.000000000              1.0000000000
    ## Lag 1000                  0.012407842              0.0217897344
    ## Lag 5000                 -0.007527296              0.0172967372
    ## Lag 10000                 0.088732711             -0.0163914620
    ## Lag 50000                 0.077918584              0.0007934218
    ##           traitFemale:traitrfit.units traitrfit:traitrfit.units
    ## Lag 0                     1.000000000                       NaN
    ## Lag 1000                  0.012407842                       NaN
    ## Lag 5000                 -0.007527296                       NaN
    ## Lag 10000                 0.088732711                       NaN
    ## Lag 50000                 0.077918584                       NaN

``` r
autocorr(mod.12$Sol)
```

    ## , , traitMale
    ## 
    ##             traitMale traitFemale   traitrfit
    ## Lag 0      1.00000000  0.12093152  0.05059305
    ## Lag 1000  -0.01721503 -0.06261559 -0.03077089
    ## Lag 5000   0.04272718  0.02540332  0.07704342
    ## Lag 10000  0.03689627  0.02204279  0.01479684
    ## Lag 50000  0.01185146 -0.02132288 -0.04785055
    ## 
    ## , , traitFemale
    ## 
    ##             traitMale  traitFemale    traitrfit
    ## Lag 0      0.12093152  1.000000000  0.230906944
    ## Lag 1000  -0.05706587 -0.005818565 -0.009555432
    ## Lag 5000   0.01465919 -0.044456165 -0.002806319
    ## Lag 10000 -0.06269602 -0.026310087  0.005717193
    ## Lag 50000 -0.01402497  0.009492818 -0.005777581
    ## 
    ## , , traitrfit
    ## 
    ##             traitMale   traitFemale    traitrfit
    ## Lag 0      0.05059305  0.2309069441  1.000000000
    ## Lag 1000   0.06061723 -0.0417729844  0.013987357
    ## Lag 5000  -0.06085983  0.0005389801 -0.069460317
    ## Lag 10000  0.01365451  0.0217137721  0.015463460
    ## Lag 50000 -0.01709034  0.0030669438 -0.003066144

#### **Covariance matrix random effects**

``` r
#mod.12$VCV
```

#### **Among and within-pair correlations and selection gradients– 2018 3 trait model**

``` r
# posteriors
posteriors_3<-as.mcmc(mod.12$VCV)
posterior.mode(posteriors_3)
```

    ##     traitMale:traitMale.NestID   traitFemale:traitMale.NestID 
    ##                   4.986134e-01                   5.552063e-02 
    ##     traitrfit:traitMale.NestID   traitMale:traitFemale.NestID 
    ##                   3.721502e-02                   5.552063e-02 
    ## traitFemale:traitFemale.NestID   traitrfit:traitFemale.NestID 
    ##                   2.194865e-03                   6.334161e-02 
    ##     traitMale:traitrfit.NestID   traitFemale:traitrfit.NestID 
    ##                   3.721502e-02                   6.334161e-02 
    ##     traitrfit:traitrfit.NestID      traitMale:traitMale.units 
    ##                   7.213663e-01                   5.839139e-01 
    ##    traitFemale:traitMale.units      traitrfit:traitMale.units 
    ##                   3.468687e-02                   4.761173e-04 
    ##    traitMale:traitFemale.units  traitFemale:traitFemale.units 
    ##                   3.468687e-02                   3.700588e-01 
    ##    traitrfit:traitFemale.units      traitMale:traitrfit.units 
    ##                  -8.352214e-05                   4.761173e-04 
    ##    traitFemale:traitrfit.units      traitrfit:traitrfit.units 
    ##                  -8.352214e-05                   9.998673e-05

``` r
# among-pair correlations
pair.correlation_3<-posteriors_3[,"traitFemale:traitMale.NestID"]/
  sqrt(posteriors_3[,"traitFemale:traitFemale.NestID"]*
         posteriors_3[,"traitMale:traitMale.NestID"])


posterior.mode(pair.correlation_3)
```

    ##      var1 
    ## 0.3760725

``` r
HPDinterval(pair.correlation_3)
```

    ##           lower     upper
    ## var1 -0.4598017 0.8081129
    ## attr(,"Probability")
    ## [1] 0.95

``` r
pair<-pair.correlation_3
pair<-as.data.frame(pair)
describe_posterior(pair)
```

    ## Summary of Posterior Distribution
    ## 
    ## Parameter | Median |        95% CI |     pd |          ROPE | % in ROPE
    ## -----------------------------------------------------------------------
    ## var1      |   0.25 | [-0.49, 0.80] | 73.80% | [-0.10, 0.10] |    16.63%

``` r
rptt<-describe_posterior(pair)#data frame for CrI and median
rpt11<-as.data.frame(pair) #posterior plot
ggplot(rpt11, aes(x = var1)) +
  geom_density(fill = "orange") +
  geom_vline(xintercept = median(rptt$Median), color="red", size=1)+
   geom_vline(xintercept = (rptt$CI_low), color="black", size=1, linetype="longdash")+
  geom_vline(xintercept = (rptt$CI_high), color="black", size=1, linetype="longdash")+
  labs(y="Density", x="Among-Pair correlation")+
  ggtitle(label ="Density Plot Among-Pair Correlation 2018",
          subtitle = "Posterior distribution plotted")+
  theme(plot.title = element_text(face = "bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5))
```

![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/cor18-1.png)<!-- -->

``` r
# within-pair correlations (residuals)
residual.correlation_3<-posteriors_3[,"traitFemale:traitMale.units"]/
  sqrt(posteriors_3[,"traitFemale:traitFemale.units"]*
         posteriors_3[,"traitMale:traitMale.units"])

posterior.mode(residual.correlation_3)
```

    ##       var1 
    ## 0.06089368

``` r
HPDinterval(residual.correlation_3)
```

    ##           lower    upper
    ## var1 -0.1838632 0.343404
    ## attr(,"Probability")
    ## [1] 0.95

``` r
wpair<-as.data.frame(residual.correlation_3)

describe_posterior(wpair)
```

    ## Summary of Posterior Distribution
    ## 
    ## Parameter | Median |        95% CI |     pd |          ROPE | % in ROPE
    ## -----------------------------------------------------------------------
    ## var1      |   0.09 | [-0.18, 0.34] | 72.40% | [-0.10, 0.10] |    45.89%

``` r
rptt<-describe_posterior(wpair)#data frame for CrI and median
rpt11<-as.data.frame(wpair) #posterior plot
ggplot(rpt11, aes(x = var1)) +
  geom_density(fill = "orange") +
  geom_vline(xintercept = median(rptt$Median), color="red", size=1)+
   geom_vline(xintercept = (rptt$CI_low), color="black", size=1, linetype="longdash")+
  geom_vline(xintercept = (rptt$CI_high), color="black", size=1, linetype="longdash")+
  labs(y="Density", x="Within-Pair correlation")+
  ggtitle(label ="Density Plot Within-Pair Correlation 2018",
          subtitle = "Posterior distribution plotted")+
  theme(plot.title = element_text(face = "bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5))
```

![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/cor18-2.png)<!-- -->

``` r
#male selection gradient
Male_sel_18<- posteriors_3[,"traitrfit:traitMale.NestID"]/
  (sqrt(posteriors_3[,"traitrfit:traitrfit.NestID"])*
     sqrt(posteriors_3[,"traitMale:traitMale.NestID"]))
posterior.mode(Male_sel_18)
```

    ##       var1 
    ## 0.06628769

``` r
HPDinterval(Male_sel_18)
```

    ##           lower     upper
    ## var1 -0.4312316 0.5082828
    ## attr(,"Probability")
    ## [1] 0.95

``` r
male18<-as.data.frame(Male_sel_18)
describe_posterior(male18)
```

    ## Summary of Posterior Distribution
    ## 
    ## Parameter | Median |        95% CI |     pd |          ROPE | % in ROPE
    ## -----------------------------------------------------------------------
    ## var1      |   0.05 | [-0.46, 0.49] | 56.40% | [-0.10, 0.10] |    32.42%

``` r
rptt<-describe_posterior(male18)#data frame for CrI and median
rpt11<-as.data.frame(male18) #posterior plot
ggplot(rpt11, aes(x = var1)) +
  geom_density(fill = "orange") +
  geom_vline(xintercept = median(rptt$Median), color="red", size=1)+
   geom_vline(xintercept = (rptt$CI_low), color="black", size=1, linetype="longdash")+
  geom_vline(xintercept = (rptt$CI_high), color="black", size=1, linetype="longdash")+
  labs(y="Density", x="Male Selection Gradient")+
  ggtitle(label ="Density Plot Male Selection Gradient 2018",
          subtitle = "Posterior distribution plotted")+
  theme(plot.title = element_text(face = "bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5))
```

![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/cor18-3.png)<!-- -->

``` r
#female selection gradient
Female_sel_18<- posteriors_3[,"traitrfit:traitFemale.NestID"]/
  (sqrt(posteriors_3[,"traitrfit:traitrfit.NestID"])*
     sqrt(posteriors_3[,"traitFemale:traitFemale.NestID"]))
posterior.mode(Female_sel_18)
```

    ##      var1 
    ## 0.3492295

``` r
HPDinterval(Female_sel_18)
```

    ##           lower     upper
    ## var1 -0.2073388 0.8736082
    ## attr(,"Probability")
    ## [1] 0.95

``` r
female18<-as.data.frame(Female_sel_18)
describe_posterior(female18)
```

    ## Summary of Posterior Distribution
    ## 
    ## Parameter | Median |        95% CI |     pd |          ROPE | % in ROPE
    ## -----------------------------------------------------------------------
    ## var1      |   0.41 | [-0.28, 0.83] | 89.40% | [-0.10, 0.10] |     9.89%

``` r
rptt<-describe_posterior(female18)#data frame for CrI and median
rpt11<-as.data.frame(female18) #posterior plot
ggplot(rpt11, aes(x = var1)) +
  geom_density(fill = "orange") +
  geom_vline(xintercept = median(rptt$Median), color="red", size=1)+
   geom_vline(xintercept = (rptt$CI_low), color="black", size=1, linetype="longdash")+
  geom_vline(xintercept = (rptt$CI_high), color="black", size=1, linetype="longdash")+
  labs(y="Density", x="Female Selection Gradient")+
  ggtitle(label ="Density Plot Female Selection Gradient 2018",
          subtitle = "Posterior distribution plotted")+
  theme(plot.title = element_text(face = "bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5))
```

![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/cor18-4.png)<!-- -->

#### **2019 3 trait model**

``` r
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
                   verbose = FALSE,
                   nitt = 1300000, thin = 1000, burnin = 300000
                   )
summary(mod.13)
```

    ## 
    ##  Iterations = 300001:1299001
    ##  Thinning interval  = 1000
    ##  Sample size  = 1000 
    ## 
    ##  DIC: 117.0247 
    ## 
    ##  G-structure:  ~us(trait):NestID
    ## 
    ##                                post.mean   l-95% CI u-95% CI eff.samp
    ## traitMale:traitMale.NestID        0.7258  1.462e-01  1.45757     1000
    ## traitFemale:traitMale.NestID     -0.2501 -6.154e-01  0.08516     1000
    ## traitrfit:traitMale.NestID        0.4915 -2.095e-01  1.29249     1000
    ## traitMale:traitFemale.NestID     -0.2501 -6.154e-01  0.08516     1000
    ## traitFemale:traitFemale.NestID    0.4556  6.702e-05  0.95110     1000
    ## traitrfit:traitFemale.NestID     -0.2576 -9.656e-01  0.33481     1000
    ## traitMale:traitrfit.NestID        0.4915 -2.095e-01  1.29249     1000
    ## traitFemale:traitrfit.NestID     -0.2576 -9.656e-01  0.33481     1000
    ## traitrfit:traitrfit.NestID        3.1909  1.752e+00  5.28423     1000
    ## 
    ##  R-structure:  ~us(trait):units
    ## 
    ##                               post.mean  l-95% CI u-95% CI eff.samp
    ## traitMale:traitMale.units     9.241e-01  0.542790 1.357233   1000.0
    ## traitFemale:traitMale.units   2.804e-01 -0.011201 0.585761   1000.0
    ## traitrfit:traitMale.units     1.410e-04 -0.007680 0.008278   1000.0
    ## traitMale:traitFemale.units   2.804e-01 -0.011201 0.585761   1000.0
    ## traitFemale:traitFemale.units 7.906e-01  0.460014 1.165358   1000.0
    ## traitrfit:traitFemale.units   6.154e-05 -0.006836 0.007909    838.9
    ## traitMale:traitrfit.units     1.410e-04 -0.007680 0.008278   1000.0
    ## traitFemale:traitrfit.units   6.154e-05 -0.006836 0.007909    838.9
    ## traitrfit:traitrfit.units     1.000e-04  0.000100 0.000100      0.0
    ## 
    ##  Location effects: cbind(Male, Female, rfit) ~ (trait - 1) 
    ## 
    ##             post.mean l-95% CI u-95% CI eff.samp  pMCMC    
    ## traitMale    -2.40534 -2.85542 -1.99259     1000 <0.001 ***
    ## traitFemale  -2.71981 -3.04057 -2.33248     1000 <0.001 ***
    ## traitrfit     0.67707  0.02085  1.37082     1000  0.048 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
plot(mod.13)
```

![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/2019Assort-1.png)<!-- -->![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/2019Assort-2.png)<!-- -->![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/2019Assort-3.png)<!-- -->![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/2019Assort-4.png)<!-- -->![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/2019Assort-5.png)<!-- -->![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/2019Assort-6.png)<!-- -->

``` r
#auto-correlation
autocorr.diag(mod.13$VCV)
```

    ##           traitMale:traitMale.NestID traitFemale:traitMale.NestID
    ## Lag 0                    1.000000000                  1.000000000
    ## Lag 1000                 0.003255865                  0.024796231
    ## Lag 5000                -0.007100073                  0.003761449
    ## Lag 10000                0.035055231                  0.024905239
    ## Lag 50000               -0.014182584                  0.001419709
    ##           traitrfit:traitMale.NestID traitMale:traitFemale.NestID
    ## Lag 0                     1.00000000                  1.000000000
    ## Lag 1000                  0.03395259                  0.024796231
    ## Lag 5000                 -0.03271509                  0.003761449
    ## Lag 10000                -0.02171407                  0.024905239
    ## Lag 50000                -0.04308871                  0.001419709
    ##           traitFemale:traitFemale.NestID traitrfit:traitFemale.NestID
    ## Lag 0                        1.000000000                  1.000000000
    ## Lag 1000                     0.014346860                 -0.011944350
    ## Lag 5000                    -0.019969292                  0.058341868
    ## Lag 10000                   -0.017577443                 -0.003304697
    ## Lag 50000                    0.004256966                 -0.008498481
    ##           traitMale:traitrfit.NestID traitFemale:traitrfit.NestID
    ## Lag 0                     1.00000000                  1.000000000
    ## Lag 1000                  0.03395259                 -0.011944350
    ## Lag 5000                 -0.03271509                  0.058341868
    ## Lag 10000                -0.02171407                 -0.003304697
    ## Lag 50000                -0.04308871                 -0.008498481
    ##           traitrfit:traitrfit.NestID traitMale:traitMale.units
    ## Lag 0                     1.00000000                1.00000000
    ## Lag 1000                  0.03747785                0.01322421
    ## Lag 5000                 -0.02387547                0.03314490
    ## Lag 10000                 0.04819820               -0.02513536
    ## Lag 50000                -0.03829091                0.04807731
    ##           traitFemale:traitMale.units traitrfit:traitMale.units
    ## Lag 0                     1.000000000               1.000000000
    ## Lag 1000                  0.020751264               0.005897439
    ## Lag 5000                 -0.043851046               0.039095023
    ## Lag 10000                -0.005730889               0.019301383
    ## Lag 50000                 0.006813945               0.017702441
    ##           traitMale:traitFemale.units traitFemale:traitFemale.units
    ## Lag 0                     1.000000000                    1.00000000
    ## Lag 1000                  0.020751264                   -0.01145308
    ## Lag 5000                 -0.043851046                   -0.00398356
    ## Lag 10000                -0.005730889                    0.03498865
    ## Lag 50000                 0.006813945                   -0.06207765
    ##           traitrfit:traitFemale.units traitMale:traitrfit.units
    ## Lag 0                     1.000000000               1.000000000
    ## Lag 1000                  0.087093881               0.005897439
    ## Lag 5000                 -0.043673541               0.039095023
    ## Lag 10000                 0.009549914               0.019301383
    ## Lag 50000                 0.010286635               0.017702441
    ##           traitFemale:traitrfit.units traitrfit:traitrfit.units
    ## Lag 0                     1.000000000                       NaN
    ## Lag 1000                  0.087093881                       NaN
    ## Lag 5000                 -0.043673541                       NaN
    ## Lag 10000                 0.009549914                       NaN
    ## Lag 50000                 0.010286635                       NaN

``` r
autocorr(mod.13$Sol)
```

    ## , , traitMale
    ## 
    ##              traitMale traitFemale     traitrfit
    ## Lag 0      1.000000000 -0.12187153  0.2297051130
    ## Lag 1000  -0.007543261  0.02879602  0.0194862663
    ## Lag 5000  -0.002679404 -0.01297298  0.0007714628
    ## Lag 10000  0.015753466 -0.03690020  0.0220050717
    ## Lag 50000 -0.018787913 -0.05965266 -0.0089199668
    ## 
    ## , , traitFemale
    ## 
    ##             traitMale traitFemale   traitrfit
    ## Lag 0     -0.12187153  1.00000000 -0.18521259
    ## Lag 1000   0.03091951 -0.01514647  0.06622646
    ## Lag 5000  -0.04985377 -0.02538271 -0.04006823
    ## Lag 10000 -0.00903294 -0.00287185 -0.03291006
    ## Lag 50000 -0.02423899 -0.01081976 -0.01197076
    ## 
    ## , , traitrfit
    ## 
    ##              traitMale  traitFemale    traitrfit
    ## Lag 0      0.229705113 -0.185212589  1.000000000
    ## Lag 1000  -0.041748842 -0.004825510 -0.028384715
    ## Lag 5000   0.025034713 -0.038845202  0.031227344
    ## Lag 10000  0.005843479 -0.007827011 -0.025991231
    ## Lag 50000 -0.006263483  0.006544551 -0.003870122

``` r
# posteriors
posteriors1<-as.mcmc(mod.13$VCV)
posterior.mode(posteriors1)
```

    ##     traitMale:traitMale.NestID   traitFemale:traitMale.NestID 
    ##                   6.828565e-01                  -2.546004e-01 
    ##     traitrfit:traitMale.NestID   traitMale:traitFemale.NestID 
    ##                   4.510442e-01                  -2.546004e-01 
    ## traitFemale:traitFemale.NestID   traitrfit:traitFemale.NestID 
    ##                   3.917944e-01                  -1.982412e-01 
    ##     traitMale:traitrfit.NestID   traitFemale:traitrfit.NestID 
    ##                   4.510442e-01                  -1.982412e-01 
    ##     traitrfit:traitrfit.NestID      traitMale:traitMale.units 
    ##                   3.007943e+00                   7.784028e-01 
    ##    traitFemale:traitMale.units      traitrfit:traitMale.units 
    ##                   2.430113e-01                  -2.658113e-03 
    ##    traitMale:traitFemale.units  traitFemale:traitFemale.units 
    ##                   2.430113e-01                   6.851933e-01 
    ##    traitrfit:traitFemale.units      traitMale:traitrfit.units 
    ##                   3.353292e-05                  -2.658113e-03 
    ##    traitFemale:traitrfit.units      traitrfit:traitrfit.units 
    ##                   3.353292e-05                   9.998673e-05

#### **Covariance matrix random effects**

``` r
#mod.13$VCV
```

#### **Among and within-pair correlations and selection gradients– 2019 3 trait model**

``` r
# among-pair correlations
pair.correlation_4<-posteriors1[,"traitFemale:traitMale.NestID"]/
  sqrt(posteriors1[,"traitFemale:traitFemale.NestID"]*
         posteriors1[,"traitMale:traitMale.NestID"])

posterior.mode(pair.correlation_4)
```

    ##       var1 
    ## -0.5378701

``` r
HPDinterval(pair.correlation_4)
```

    ##           lower      upper
    ## var1 -0.9391147 0.07140852
    ## attr(,"Probability")
    ## [1] 0.95

``` r
pair<-as.data.frame(pair.correlation_4)
describe_posterior(pair)
```

    ## Summary of Posterior Distribution
    ## 
    ## Parameter | Median |        95% CI |     pd |          ROPE | % in ROPE
    ## -----------------------------------------------------------------------
    ## var1      |  -0.52 | [-0.92, 0.15] | 94.00% | [-0.10, 0.10] |     8.00%

``` r
rptt<-describe_posterior(pair)#data frame for CrI and median
rpt11<-as.data.frame(pair) #posterior plot
ggplot(rpt11, aes(x = var1)) +
  geom_density(fill = "orange") +
  geom_vline(xintercept = median(rptt$Median), color="red", size=1)+
   geom_vline(xintercept = (rptt$CI_low), color="black", size=1, linetype="longdash")+
  geom_vline(xintercept = (rptt$CI_high), color="black", size=1, linetype="longdash")+
  labs(y="Density", x="Among-Pair correlation")+
  ggtitle(label ="Density Plot Among-Pair Correlation 2019",
          subtitle = "Posterior distribution plotted")+
  theme(plot.title = element_text(face = "bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5))
```

![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/cor19-1.png)<!-- -->

``` r
# within-pair correlations (residuals)
residual.correlation4<-posteriors1[,"traitFemale:traitMale.units"]/
  sqrt(posteriors1[,"traitFemale:traitFemale.units"]*
         posteriors1[,"traitMale:traitMale.units"])

posterior.mode(residual.correlation4)
```

    ##     var1 
    ## 0.415727

``` r
HPDinterval(residual.correlation4)
```

    ##           lower     upper
    ## var1 0.05122956 0.6062937
    ## attr(,"Probability")
    ## [1] 0.95

``` r
wpair<-as.data.frame(residual.correlation4)
describe_posterior(wpair)
```

    ## Summary of Posterior Distribution
    ## 
    ## Parameter | Median |       95% CI |     pd |          ROPE | % in ROPE
    ## ----------------------------------------------------------------------
    ## var1      |   0.34 | [0.00, 0.58] | 97.50% | [-0.10, 0.10] |     4.11%

``` r
rptt<-describe_posterior(wpair)#data frame for CrI and median
rpt11<-as.data.frame(wpair) #posterior plot
ggplot(rpt11, aes(x = var1)) +
  geom_density(fill = "orange") +
  geom_vline(xintercept = median(rptt$Median), color="red", size=1)+
   geom_vline(xintercept = (rptt$CI_low), color="black", size=1, linetype="longdash")+
  geom_vline(xintercept = (rptt$CI_high), color="black", size=1, linetype="longdash")+
  labs(y="Density", x="Within-Pair correlation")+
  ggtitle(label ="Density Plot Within-Pair Correlation 2019",
          subtitle = "Posterior distribution plotted")+
  theme(plot.title = element_text(face = "bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5))
```

![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/cor19-2.png)<!-- -->

``` r
#male selection gradient
Male_sel_19<- posteriors1[,"traitrfit:traitMale.NestID"]/
  (sqrt(posteriors1[,"traitrfit:traitrfit.NestID"])*
     sqrt(posteriors1[,"traitMale:traitMale.NestID"]))
posterior.mode(Male_sel_19)
```

    ##      var1 
    ## 0.3814246

``` r
HPDinterval(Male_sel_19)
```

    ##            lower     upper
    ## var1 -0.06968465 0.7641123
    ## attr(,"Probability")
    ## [1] 0.95

``` r
male19<-as.data.frame(Male_sel_19)
describe_posterior(male19)
```

    ## Summary of Posterior Distribution
    ## 
    ## Parameter | Median |        95% CI |     pd |          ROPE | % in ROPE
    ## -----------------------------------------------------------------------
    ## var1      |   0.35 | [-0.15, 0.70] | 92.30% | [-0.10, 0.10] |    11.68%

``` r
rptt<-describe_posterior(male19)#data frame for CrI and median
rpt11<-as.data.frame(male19) #posterior plot
ggplot(rpt11, aes(x = var1)) +
  geom_density(fill = "orange") +
  geom_vline(xintercept = median(rptt$Median), color="red", size=1)+
   geom_vline(xintercept = (rptt$CI_low), color="black", size=1, linetype="longdash")+
  geom_vline(xintercept = (rptt$CI_high), color="black", size=1, linetype="longdash")+
  labs(y="Density", x="Male Selection Gradient")+
  ggtitle(label ="Density Plot Male Selection Gradient 2019",
          subtitle = "Posterior distribution plotted")+
  theme(plot.title = element_text(face = "bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5))
```

![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/cor19-3.png)<!-- -->

``` r
#female selection gradient
Female_sel_19<- posteriors1[,"traitrfit:traitFemale.NestID"]/
  (sqrt(posteriors1[,"traitrfit:traitrfit.NestID"])*
     sqrt(posteriors1[,"traitFemale:traitFemale.NestID"]))
posterior.mode(Female_sel_19)
```

    ##       var1 
    ## -0.2934552

``` r
HPDinterval(Female_sel_19)
```

    ##          lower     upper
    ## var1 -0.714596 0.2789848
    ## attr(,"Probability")
    ## [1] 0.95

``` r
female19<-as.data.frame(Female_sel_19)
describe_posterior(female19)
```

    ## Summary of Posterior Distribution
    ## 
    ## Parameter | Median |        95% CI |     pd |          ROPE | % in ROPE
    ## -----------------------------------------------------------------------
    ## var1      |  -0.24 | [-0.70, 0.32] | 80.90% | [-0.10, 0.10] |    19.37%

``` r
rptt<-describe_posterior(female19)#data frame for CrI and median
rpt11<-as.data.frame(female19) #posterior plot
ggplot(rpt11, aes(x = var1)) +
  geom_density(fill = "orange") +
  geom_vline(xintercept = median(rptt$Median), color="red", size=1)+
   geom_vline(xintercept = (rptt$CI_low), color="black", size=1, linetype="longdash")+
  geom_vline(xintercept = (rptt$CI_high), color="black", size=1, linetype="longdash")+
  labs(y="Density", x="Female Selection Gradient")+
  ggtitle(label ="Density Plot Female Selection Gradient 2019",
          subtitle = "Posterior distribution plotted")+
  theme(plot.title = element_text(face = "bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5))
```

![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/cor19-4.png)<!-- -->

#### **Comparing 2018 and 2019**

``` r
##is male selection different across years?
male_diff<-Male_sel_18-Male_sel_19
posterior.mode(male_diff)
```

    ##       var1 
    ## -0.3733127

``` r
HPDinterval(male_diff)
```

    ##           lower     upper
    ## var1 -0.9111791 0.3752959
    ## attr(,"Probability")
    ## [1] 0.95

``` r
male_pos<-ifelse(male_diff>0,1,0)
male_neg<-ifelse(male_diff>0,0,1)
sum(male_pos)
```

    ## [1] 170

``` r
sum(male_neg)
```

    ## [1] 830

``` r
##is female section different across years?
female_diff<-Female_sel_18-Female_sel_19
posterior.mode(female_diff)
```

    ##      var1 
    ## 0.6382944

``` r
HPDinterval(female_diff)
```

    ##           lower    upper
    ## var1 -0.1698064 1.320827
    ## attr(,"Probability")
    ## [1] 0.95

``` r
female_pos<-ifelse(female_diff>0,1,0)
female_neg<-ifelse(female_diff>0,0,1)
sum(female_pos)
```

    ## [1] 941

``` r
sum(female_neg)
```

    ## [1] 59

``` r
#do among pair correlations differ across years?
among_diff<-pair.correlation_3-pair.correlation_4
posterior.mode(among_diff)
```

    ##      var1 
    ## 0.7204608

``` r
HPDinterval(among_diff)
```

    ##           lower    upper
    ## var1 -0.1103407 1.621404
    ## attr(,"Probability")
    ## [1] 0.95

``` r
among_pos<-ifelse(among_diff>0,1,0)
among_neg<-ifelse(among_diff>0,0,1)
sum(among_pos)
```

    ## [1] 933

``` r
sum(among_neg)
```

    ## [1] 67

``` r
#do within pair correlations differ across years?
within_diff<-residual.correlation4-residual.correlation_3
posterior.mode(within_diff)
```

    ##      var1 
    ## 0.1800889

``` r
HPDinterval(within_diff)
```

    ##           lower     upper
    ## var1 -0.1285764 0.6750763
    ## attr(,"Probability")
    ## [1] 0.95

``` r
within_pos<-ifelse(within_diff>0,1,0)
within_neg<-ifelse(within_diff>0,0,1)
sum(within_pos)
```

    ## [1] 889

``` r
sum(within_neg)
```

    ## [1] 111

#### **2018/2019 3 trait model**

``` r
#pooled years----

require(readr)
require(tidyr)
require(MCMCglmm)
#data
data2018<-read.csv("C:/Users/nickg/OneDrive/Desktop/R projects/krmp_nest-defense/data/data_2018.csv")
data2019<-read.csv("C:/Users/nickg/OneDrive/Desktop/R projects/krmp_nest-defense/data/data_2019.csv")



require(dplyr)
data2018<-data2018%>%mutate(Site=SiteID)
require(tidyr)
data2018<-data2018%>% 
  unite(SiteID_Series,c("SiteID", "Year"))
data2018<-data2018%>%mutate(Year="2018")



data2019<-data2019%>%mutate(Site=SiteID)
data2019<-data2019%>% 
  unite(SiteID_Series,c("SiteID", "Year"))
data2019<-data2019%>%mutate(Year="2019")


data_2<-full_join(data2018, data2019) #join datasets
data_2$Male<-log(data_2$male_raw+1)*-1 #log transform
data_2$Female<-log(data_2$female_raw+1)*-1 #log transform
#prior for three traits- from Houslay tutorial
prior_E_B_fit_1px = list(R = list(V = diag(c(1,1,0.0001),3,3), nu = 1.002, fix = 3),
                         G = list(G1 = list(V = diag(3), nu = 3,
                                            alpha.mu = rep(0,3),
                                            alpha.V = diag(25^2,3,3))))
#set residual variance for fitness near 0
prior_E_B_fit_1px$R$V[3,3]<-0.0001

#2018/2019 3 trait model ----
mod.122 <- MCMCglmm(cbind(Male, Female, rfit) ~ (trait-1),  
                   random = ~us(trait):SiteID_Series ,
                   rcov = ~us(trait):units, 
                   family = c("gaussian", "gaussian", "gaussian"),
                   data=data_2, 
                   prior = prior_E_B_fit_1px, 
                   verbose = FALSE,
                   nitt = 1300000, thin = 1000, burnin = 300000
)

plot(mod.122)
```

![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/pooled-1.png)<!-- -->![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/pooled-2.png)<!-- -->![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/pooled-3.png)<!-- -->![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/pooled-4.png)<!-- -->![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/pooled-5.png)<!-- -->![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/pooled-6.png)<!-- -->

``` r
summary(mod.122)
```

    ## 
    ##  Iterations = 300001:1299001
    ##  Thinning interval  = 1000
    ##  Sample size  = 1000 
    ## 
    ##  DIC: 234.7829 
    ## 
    ##  G-structure:  ~us(trait):SiteID_Series
    ## 
    ##                                       post.mean l-95% CI u-95% CI eff.samp
    ## traitMale:traitMale.SiteID_Series       0.59125  0.25028  0.99703     1149
    ## traitFemale:traitMale.SiteID_Series    -0.11205 -0.28679  0.05139     1000
    ## traitrfit:traitMale.SiteID_Series       0.25048 -0.10984  0.62311     1000
    ## traitMale:traitFemale.SiteID_Series    -0.11205 -0.28679  0.05139     1000
    ## traitFemale:traitFemale.SiteID_Series   0.25587  0.07626  0.47732     1000
    ## traitrfit:traitFemale.SiteID_Series    -0.06285 -0.34123  0.18078     1000
    ## traitMale:traitrfit.SiteID_Series       0.25048 -0.10984  0.62311     1000
    ## traitFemale:traitrfit.SiteID_Series    -0.06285 -0.34123  0.18078     1000
    ## traitrfit:traitrfit.SiteID_Series       1.98212  1.29205  2.77215     1285
    ## 
    ##  R-structure:  ~us(trait):units
    ## 
    ##                                post.mean  l-95% CI u-95% CI eff.samp
    ## traitMale:traitMale.units      7.225e-01  0.521689 0.958526   1000.0
    ## traitFemale:traitMale.units    1.536e-01  0.017006 0.285767   1000.0
    ## traitrfit:traitMale.units      2.669e-05 -0.006887 0.006828    864.1
    ## traitMale:traitFemale.units    1.536e-01  0.017006 0.285767   1000.0
    ## traitFemale:traitFemale.units  5.286e-01  0.383953 0.718333   1098.6
    ## traitrfit:traitFemale.units   -6.150e-05 -0.005495 0.005614   1000.0
    ## traitMale:traitrfit.units      2.669e-05 -0.006887 0.006828    864.1
    ## traitFemale:traitrfit.units   -6.150e-05 -0.005495 0.005614   1000.0
    ## traitrfit:traitrfit.units      1.000e-04  0.000100 0.000100      0.0
    ## 
    ##  Location effects: cbind(Male, Female, rfit) ~ (trait - 1) 
    ## 
    ##             post.mean l-95% CI u-95% CI eff.samp  pMCMC    
    ## traitMale     -2.3746  -2.6351  -2.1199     1000 <0.001 ***
    ## traitFemale   -2.6514  -2.8521  -2.4687     1000 <0.001 ***
    ## traitrfit      0.7299   0.3875   1.1433     1000 <0.001 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#auto-correlation
autocorr.diag(mod.122$VCV)
```

    ##           traitMale:traitMale.SiteID_Series traitFemale:traitMale.SiteID_Series
    ## Lag 0                            1.00000000                          1.00000000
    ## Lag 1000                        -0.06976720                         -0.04103222
    ## Lag 5000                        -0.01432652                         -0.03100777
    ## Lag 10000                       -0.01023715                         -0.02423127
    ## Lag 50000                        0.04013006                         -0.02989298
    ##           traitrfit:traitMale.SiteID_Series traitMale:traitFemale.SiteID_Series
    ## Lag 0                            1.00000000                          1.00000000
    ## Lag 1000                         0.04153844                         -0.04103222
    ## Lag 5000                         0.01861047                         -0.03100777
    ## Lag 10000                        0.02128101                         -0.02423127
    ## Lag 50000                       -0.02428383                         -0.02989298
    ##           traitFemale:traitFemale.SiteID_Series
    ## Lag 0                              1.0000000000
    ## Lag 1000                           0.0143721732
    ## Lag 5000                           0.0377051243
    ## Lag 10000                          0.0050527855
    ## Lag 50000                         -0.0004002886
    ##           traitrfit:traitFemale.SiteID_Series traitMale:traitrfit.SiteID_Series
    ## Lag 0                             1.000000000                        1.00000000
    ## Lag 1000                         -0.004485045                        0.04153844
    ## Lag 5000                          0.040445127                        0.01861047
    ## Lag 10000                        -0.006764720                        0.02128101
    ## Lag 50000                         0.012051408                       -0.02428383
    ##           traitFemale:traitrfit.SiteID_Series traitrfit:traitrfit.SiteID_Series
    ## Lag 0                             1.000000000                       1.000000000
    ## Lag 1000                         -0.004485045                       0.007670160
    ## Lag 5000                          0.040445127                      -0.014683806
    ## Lag 10000                        -0.006764720                      -0.004022242
    ## Lag 50000                         0.012051408                       0.004866768
    ##           traitMale:traitMale.units traitFemale:traitMale.units
    ## Lag 0                    1.00000000                  1.00000000
    ## Lag 1000                 0.01057613                  0.01962769
    ## Lag 5000                 0.02058026                  0.05427485
    ## Lag 10000                0.01693657                  0.01321521
    ## Lag 50000               -0.03776340                 -0.01553351
    ##           traitrfit:traitMale.units traitMale:traitFemale.units
    ## Lag 0                  1.0000000000                  1.00000000
    ## Lag 1000               0.0434939006                  0.01962769
    ## Lag 5000               0.0518288382                  0.05427485
    ## Lag 10000              0.0002471502                  0.01321521
    ## Lag 50000             -0.0259692949                 -0.01553351
    ##           traitFemale:traitFemale.units traitrfit:traitFemale.units
    ## Lag 0                        1.00000000                 1.000000000
    ## Lag 1000                    -0.04749889                 0.005000181
    ## Lag 5000                    -0.05281134                 0.087934358
    ## Lag 10000                   -0.02547927                -0.024096830
    ## Lag 50000                   -0.02468861                -0.009687367
    ##           traitMale:traitrfit.units traitFemale:traitrfit.units
    ## Lag 0                  1.0000000000                 1.000000000
    ## Lag 1000               0.0434939006                 0.005000181
    ## Lag 5000               0.0518288382                 0.087934358
    ## Lag 10000              0.0002471502                -0.024096830
    ## Lag 50000             -0.0259692949                -0.009687367
    ##           traitrfit:traitrfit.units
    ## Lag 0                           NaN
    ## Lag 1000                        NaN
    ## Lag 5000                        NaN
    ## Lag 10000                       NaN
    ## Lag 50000                       NaN

``` r
autocorr(mod.122$Sol)
```

    ## , , traitMale
    ## 
    ##              traitMale   traitFemale    traitrfit
    ## Lag 0      1.000000000 -0.0398548321  0.178282505
    ## Lag 1000  -0.023816162  0.0007517514 -0.029088385
    ## Lag 5000  -0.002640596 -0.0272845219  0.004307238
    ## Lag 10000 -0.015177305 -0.0159052860  0.014510003
    ## Lag 50000 -0.020516345 -0.0266200366  0.019657648
    ## 
    ## , , traitFemale
    ## 
    ##              traitMale traitFemale   traitrfit
    ## Lag 0     -0.039854832  1.00000000 -0.02021594
    ## Lag 1000   0.009956558 -0.02090325  0.03122404
    ## Lag 5000  -0.028197495 -0.01587403 -0.04809554
    ## Lag 10000 -0.017603056 -0.01718580 -0.02164424
    ## Lag 50000  0.013298876 -0.08013802 -0.01048096
    ## 
    ## , , traitrfit
    ## 
    ##             traitMale  traitFemale    traitrfit
    ## Lag 0      0.17828251 -0.020215939  1.000000000
    ## Lag 1000   0.01064287 -0.050650581 -0.000100731
    ## Lag 5000   0.01945395 -0.009441877  0.006639041
    ## Lag 10000 -0.04598225 -0.006769760 -0.001129813
    ## Lag 50000 -0.02755478 -0.009700478 -0.036600284

#### **Covariance matrix random effects**

``` r
#mod.122$VCV
```

##### **Among and within-pair correlations and selection gradients– 2018/2019 3 trait model**

``` r
# posteriors
require(bayestestR)
require(ggplot2)
posteriors_3<-as.mcmc(mod.122$VCV)


# among-pair correlations
pair.correlation_33<-posteriors_3[,"traitFemale:traitMale.SiteID_Series"]/
  sqrt(posteriors_3[,"traitFemale:traitFemale.SiteID_Series"]*
         posteriors_3[,"traitMale:traitMale.SiteID_Series"])

posterior.mode(pair.correlation_33)
```

    ##      var1 
    ## -0.399878

``` r
HPDinterval(pair.correlation_33)
```

    ##           lower     upper
    ## var1 -0.7511533 0.1186297
    ## attr(,"Probability")
    ## [1] 0.95

``` r
pair<-as.data.frame(pair.correlation_33)
describe_posterior(pair)
```

    ## Summary of Posterior Distribution
    ## 
    ## Parameter | Median |        95% CI |     pd |          ROPE | % in ROPE
    ## -----------------------------------------------------------------------
    ## var1      |  -0.32 | [-0.73, 0.15] | 90.00% | [-0.10, 0.10] |    14.53%

``` r
rptt<-describe_posterior(pair)#data frame for CrI and median
rpt11<-as.data.frame(pair) #posterior plot
ggplot(rpt11, aes(x = var1)) +
  geom_density(fill = "orange") +
  geom_vline(xintercept = median(rptt$Median), color="red", size=1)+
   geom_vline(xintercept = (rptt$CI_low), color="black", size=1, linetype="longdash")+
  geom_vline(xintercept = (rptt$CI_high), color="black", size=1, linetype="longdash")+
  labs(y="Density", x="Among-Pair correlation")+
  ggtitle(label ="Density Plot Among-Pair Correlation 2018/2019",
          subtitle = "Posterior distribution plotted")+
  theme(plot.title = element_text(face = "bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5))
```

![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/plotmod122-1.png)<!-- -->

``` r
# within-pair correlations (residuals)
residual.correlation_33<-posteriors_3[,"traitFemale:traitMale.units"]/
  sqrt(posteriors_3[,"traitFemale:traitFemale.units"]*
         posteriors_3[,"traitMale:traitMale.units"])

posterior.mode(residual.correlation_33)
```

    ##      var1 
    ## 0.2651647

``` r
HPDinterval(residual.correlation_33)
```

    ##           lower     upper
    ## var1 0.03675965 0.4307064
    ## attr(,"Probability")
    ## [1] 0.95

``` r
wpair<-as.data.frame(residual.correlation_33)
describe_posterior(wpair)
```

    ## Summary of Posterior Distribution
    ## 
    ## Parameter | Median |       95% CI |     pd |          ROPE | % in ROPE
    ## ----------------------------------------------------------------------
    ## var1      |   0.26 | [0.04, 0.44] | 99.10% | [-0.10, 0.10] |     5.37%

``` r
rptt<-describe_posterior(wpair)#data frame for CrI and median
rpt11<-as.data.frame(wpair) #posterior plot
ggplot(rpt11, aes(x = var1)) +
  geom_density(fill = "orange") +
  geom_vline(xintercept = median(rptt$Median), color="red", size=1)+
   geom_vline(xintercept = (rptt$CI_low), color="black", size=1, linetype="longdash")+
  geom_vline(xintercept = (rptt$CI_high), color="black", size=1, linetype="longdash")+
  labs(y="Density", x="Within-Pair correlation")+
  ggtitle(label ="Density Plot Within-Pair Correlation 2018/2019",
          subtitle = "Posterior distribution plotted")+
  theme(plot.title = element_text(face = "bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5))
```

![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/plotmod122-2.png)<!-- -->

``` r
#male selection gradient
Male_sel<- posteriors_3[,"traitrfit:traitMale.SiteID_Series"]/
  (sqrt(posteriors_3[,"traitrfit:traitrfit.SiteID_Series"])*
     sqrt(posteriors_3[,"traitMale:traitMale.SiteID_Series"]))
posterior.mode(Male_sel)
```

    ##      var1 
    ## 0.1707595

``` r
HPDinterval(Male_sel)
```

    ##            lower     upper
    ## var1 -0.05233826 0.5507422
    ## attr(,"Probability")
    ## [1] 0.95

``` r
male<-as.data.frame(Male_sel)
describe_posterior(male)
```

    ## Summary of Posterior Distribution
    ## 
    ## Parameter | Median |        95% CI |     pd |          ROPE | % in ROPE
    ## -----------------------------------------------------------------------
    ## var1      |   0.24 | [-0.08, 0.53] | 92.40% | [-0.10, 0.10] |    17.79%

``` r
rptt<-describe_posterior(male)#data frame for CrI and median
rpt11<-as.data.frame(male) #posterior plot
ggplot(rpt11, aes(x = var1)) +
  geom_density(fill = "orange") +
  geom_vline(xintercept = median(rptt$Median), color="red", size=1)+
   geom_vline(xintercept = (rptt$CI_low), color="black", size=1, linetype="longdash")+
  geom_vline(xintercept = (rptt$CI_high), color="black", size=1, linetype="longdash")+
  labs(y="Density", x="Male Selection Gradient")+
  ggtitle(label ="Density Plot Male Selection Gradient 2018/2019",
          subtitle = "Posterior distribution plotted")+
  theme(plot.title = element_text(face = "bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5))
```

![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/plotmod122-3.png)<!-- -->

``` r
#female selection gradient
Female_sel<- posteriors_3[,"traitrfit:traitFemale.SiteID_Series"]/
  (sqrt(posteriors_3[,"traitrfit:traitrfit.SiteID_Series"])*
     sqrt(posteriors_3[,"traitFemale:traitFemale.SiteID_Series"]))
posterior.mode(Female_sel)
```

    ##        var1 
    ## -0.08397669

``` r
HPDinterval(Female_sel)
```

    ##           lower     upper
    ## var1 -0.4781202 0.2471888
    ## attr(,"Probability")
    ## [1] 0.95

``` r
female<-as.data.frame(Female_sel)
describe_posterior(female)
```

    ## Summary of Posterior Distribution
    ## 
    ## Parameter | Median |        95% CI |     pd |          ROPE | % in ROPE
    ## -----------------------------------------------------------------------
    ## var1      |  -0.09 | [-0.47, 0.26] | 67.30% | [-0.10, 0.10] |    39.58%

``` r
rptt<-describe_posterior(female)#data frame for CrI and median
rpt11<-as.data.frame(female) #posterior plot
ggplot(rpt11, aes(x = var1)) +
  geom_density(fill = "orange") +
  geom_vline(xintercept = median(rptt$Median), color="red", size=1)+
   geom_vline(xintercept = (rptt$CI_low), color="black", size=1, linetype="longdash")+
  geom_vline(xintercept = (rptt$CI_high), color="black", size=1, linetype="longdash")+
  labs(y="Density", x="Female Selection Gradient")+
  ggtitle(label ="Density Plot Female Selection Gradient 2018/2019",
          subtitle = "Posterior distribution plotted")+
  theme(plot.title = element_text(face = "bold", hjust=0.5), plot.subtitle=element_text(hjust=0.5))
```

![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/plotmod122-4.png)<!-- -->

##### **Correlation Plot**

``` r
### create data frame and within-subject center for plot

require(tidyr)

df_mcmc_cors <- data_frame(Traits = c("Among-pair (2018)",    
                                      "Within-pair (2018)",
                                      "Among-pair (2019)",
                                      "Within-pair (2019)",
                                      "Among-pair (Pooled)",
                                      "Within-pair (Pooled)"
                                        ),
                           Estimate = c(posterior.mode(pair.correlation_3),
                                        posterior.mode(residual.correlation_3),
                                        posterior.mode(pair.correlation_4),
                                        posterior.mode(residual.correlation4),
                                        posterior.mode(pair.correlation_33),
                                        posterior.mode(residual.correlation_33)),
                           Lower = c(HPDinterval(pair.correlation_3)[,"lower"],
                                     HPDinterval(residual.correlation_3)[,"lower"],
                                     HPDinterval(pair.correlation_4)[,"lower"],
                                     HPDinterval(residual.correlation4)[,"lower"],
                                     HPDinterval(pair.correlation_33)[,"lower"],
                                     HPDinterval(residual.correlation_33)[,"lower"]
                                     ),
                           Upper = c(HPDinterval(pair.correlation_3)[,"upper"],
                                     HPDinterval(residual.correlation_3)[,"upper"],
                                     HPDinterval(pair.correlation_4)[,"upper"],
                                     HPDinterval(residual.correlation4)[,"upper"],
                                     HPDinterval(pair.correlation_33)[,"upper"],
                                     HPDinterval(residual.correlation_33)[,"upper"]
                                     ))




correlation1<-ggplot(df_mcmc_cors, aes(x = Traits, y = Estimate)) +
  geom_pointrange(aes(ymin = Lower,
                      ymax = Upper)) +
  geom_hline(yintercept = 0,
             linetype = "dotted",alpha = 0.35, size=0.3) +
  labs(x = "",
       y = "Correlation (Estimate +/- 95% CIs)") +
  ylim(-1,1) +
  coord_flip() + theme_classic()





df_mcmc_cors1 <- data_frame(Traits = c("Female (2018)",    
                                      "Male (2018)",
                                      "Female (2019)",
                                      "Male (2019)",
                                      "Female (Pooled)",
                                      "Male (Pooled)"
),
Estimate = c(posterior.mode(Female_sel_18),
             posterior.mode(Male_sel_18),
             posterior.mode(Female_sel_19),
             posterior.mode(Male_sel_19),
             posterior.mode(Female_sel),
             posterior.mode(Male_sel)
             ),
Lower = c(HPDinterval(Female_sel_18)[,"lower"],
          HPDinterval(Male_sel_18)[,"lower"],
          HPDinterval(Female_sel_19)[,"lower"],
          HPDinterval(Male_sel_19)[,"lower"],
          HPDinterval(Female_sel)[,"lower"],
          HPDinterval(Male_sel)[,"lower"]
),
Upper = c(HPDinterval(Female_sel_18)[,"upper"],
          HPDinterval(Male_sel_18)[,"upper"],
          HPDinterval(Female_sel_19)[,"upper"],
          HPDinterval(Male_sel_19)[,"upper"],
          HPDinterval(Female_sel)[,"upper"],
          HPDinterval(Male_sel)[,"upper"]
          ))




correlation2<-ggplot(df_mcmc_cors1, aes(x = Traits, y = Estimate)) +
  geom_pointrange(aes(ymin = Lower,
                      ymax = Upper)) +
  geom_hline(yintercept = 0,
             linetype = "dotted",alpha = 0.35, size=0.3) +
  labs(x = "",
       y = "Correlation (Estimate +/- 95% CIs)") +
  ylim(-1,1) +
  coord_flip() + theme_classic()

require(cowplot)
plot_grid(correlation1,correlation2, labels = c("A.", "B."), nrow = 1)
```

![](State_dependence_repeatability_NestDefense_Rmd_files/figure-gfm/corrplot-1.png)<!-- -->
