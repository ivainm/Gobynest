#Analysing data from nest pictures and male phenotypes


#setwd("C://Users//ivainm//Working_Document//DYNAMAR//Manuscripts//Nests//data")
setwd("C://Users//ivmar9032//Documents//DYNAMAR local//Manuscripts//Nests//data//Gobynest")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# packages ----
library(ggplot2)
library(ggExtra)
library(emmeans)
#library(reshape)
#library(lsmeans)
library(effects)
#library(multcomp)
library(ggpubr)
library(lubridate)
library(car)
library(lme4)
#library(brms)
#library(gamlss)
library(betareg)
library(plyr)
library(dplyr)
library(corrplot)
#library(postHoc)
#library(rjags)
#library(zoib)
library(MuMIn)
library(glmmTMB)

#install.packages(c("StanHeaders", "rstan"))
#install.packages("cmdstanr")
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(usethis)
edit_git_config()
edit_git_ignore()
use_git()

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# 1. Nest data import and format ----
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

## 1.1  nest data import----------------------------------------------------------
nestdat<-read.delim("fullnest2.txt", header=T)
nestdat$area<- factor(nestdat$area, levels=c('Kristineberg', 'Arendal', 'Austevoll', 'Hitra', 'Helligvaer', 'Ringstad'))
nestdat$location<-as.factor(nestdat$location)
table(nestdat$Nest_size)
nestdat<-nestdat[which(nestdat$Nest_size=="l" |nestdat$Nest_size=="s" ) ,]
nestdat$Nest_size<-as.factor(nestdat$Nest_size)
nestdat$male<-as.factor(nestdat$male)
nestdat$Eyed_analysis<-as.factor(nestdat$Eyed_analysis)
nestdat$Live_eggs_analysis<-as.factor(nestdat$Live_eggs_analysis)
nestdat$Timeseason<-as.factor(nestdat$Timepoint)
levels(nestdat$Timeseason)<-c("Peak season","Late season")
levels(nestdat$Nest_size)<-c("Large nest","Small nest")
str(nestdat)

#Add latitude to each population
Coord<-read.delim("Coord.txt", header=T)
#Add environmental data
Env<-read.delim("env_info.txt", header=T, sep=",")

#add latitude info
nestdat<-nestdat %>%
  mutate(latitude=case_when(
    area=="Kristineberg" ~ Coord[Coord$Area=="Kristineberg",]$Lat,
    area=="Arendal" ~ Coord[Coord$Area=="Arendal",]$Lat,
    area=="Austevoll" ~ Coord[Coord$Area=="Austevoll",]$Lat,
    area=="Hitra" ~ Coord[Coord$Area=="Hitra",]$Lat,
    area=="Helligvaer" ~ Coord[Coord$Area=="Helligvaer",]$Lat,
    area=="Ringstad" ~ Coord[Coord$Area=="Ringstad",]$Lat))

#add degreedays fixed info
nestdat<-nestdat %>%
  mutate(DDF=case_when(
    area=="Kristineberg" ~ Env[Env$Area=="Kristineberg",]$Ddaysfixed ,
    area=="Arendal" ~ Env[Env$Area=="Arendal",]$Ddaysfixed ,
    area=="Austevoll" ~ Env[Env$Area=="Austevoll",]$Ddaysfixed ,
    area=="Hitra" ~ Env[Env$Area=="Hitra",]$Ddaysfixed ,
    area=="Helligvaer" ~ Env[Env$Area=="Helligvaer",]$Ddaysfixed ,
    area=="Ringstad" ~ Env[Env$Area=="Ringstad",]$Ddaysfixed ))

ggplot(nestdat,aes(x=area,y=DDF))+geom_point()
### calculation of new variables ---------------------------------------
#The only consistent thing between nests is the internal surface. 
#Unfortunately, the grid patterns used on the lining sheet were not always the same size
#Thus, for each nest we need to calculate the total nest surface, in units of grid squares, which is:
# square_height*square_length
#This total surface will be used to standardize other surfaces
nestdat$nestsurf<-nestdat$square_height*nestdat$square_length

#Total brood area relative to nest size (including empty spaces and overgrown spaces where presumably there were eggs before)
nestdat$TotBroodFrac<-nestdat$TotBrood_surface/(nestdat$X1sq_surface*nestdat$nestsurf)
#Surface of  Large nests in cm^2: 17.5*2*3.14159 | of small nest: 8*1.6*3.14159
nestdat$surfcm<-NA
nestdat[nestdat$Nest_size=="Large nest",]$surfcm<-17.5*2*3.14159
nestdat[nestdat$Nest_size=="Small nest",]$surfcm<-8*1.6*3.14159
#Total brood area in cm^2
nestdat$TotBroodCm<-nestdat$TotBroodFrac*nestdat$surfcm

# The cannibalised area and overgrown area in fraction of nest and cm^2
nestdat$CannibFrac<-nestdat$Cann_surface/(nestdat$X1sq_surface*nestdat$nestsurf)
nestdat$CannibCm<-nestdat$CannibFrac*nestdat$surfcm

nestdat$OvergrownFrac<-nestdat$Grown_surface/(nestdat$X1sq_surface*nestdat$nestsurf)
nestdat$OvergrownCm<-nestdat$OvergrownFrac*nestdat$surfcm

# The area still covered with eggs within the brood area, in fraction of nest surface and cm^2
nestdat$BroodFrac<-nestdat$TotBroodFrac-(nestdat$CannibFrac+nestdat$OvergrownFrac)
nestdat$BroodCm<-nestdat$TotBroodCm-(nestdat$CannibCm+nestdat$OvergrownCm)

#Cannibalised and overgrwon area, in proportion of the total brood surface
nestdat$PropCann<-nestdat$CannibCm/nestdat$TotBroodCm
nestdat$PropGrown<-nestdat$OvergrownCm/nestdat$TotBroodCm

## Egg number in total brood area
nestdat$Egg_Nr<-nestdat$Egg1+nestdat$Egg2+nestdat$Egg3+nestdat$Egg4+nestdat$Egg5
nestdat$Totegg<-(nestdat$TotBrood_surface/(5*nestdat$Egg_surface))*nestdat$Egg_Nr
# Egg number left after cannibalism and overgrowth
nestdat$Eggs<-((nestdat$TotBrood_surface-(nestdat$Cann_surface+nestdat$Grown_surface))/(5*nestdat$Egg_surface))*nestdat$Egg_Nr


head(nestdat)
str(nestdat)

#Subset for when we want to work only with nest that have a male measured
nestmale<-nestdat[!is.na(nestdat$length),]

# Estimate condition as residuals from a length~weightregression over all pops
lmcondall<-lm(log(weight)~log(length), data=nestdat)
summary(nestdat$weight)
#add residuals to the dataset
hist(residuals (lmcondall))
nestdat$con_all<-NA
nestdat[rownames(nestdat) %in% rownames(nestmale),]$con_all<-residuals (lmcondall)


## 1.2 Egg size data import-------------------------------------

eggiz<-read.delim("eggsize.txt", header=T)
eggiz$Area<-as.factor(eggiz$Area)
eggiz$Area<- factor(eggiz$Area, levels=c('Kristineberg', 'Arendal','Austevoll','Hitra','Helligvaer',  'Ringstad'))
eggiz$Nest_size_pic<-as.factor(eggiz$Nest_size_pic)
levels(eggiz$Nest_size_pic)<-c("Large","Small")
eggiz$eyed<-as.factor(eggiz$eyed)
levels(eggiz$eyed)<-c("No","Yes")
eggiz$Timeseason<-as.factor(eggiz$Timeseason)
eggiz$Timeseason<- factor(eggiz$Timeseason, levels=c('Peak season', 'Late season'))
eggiz$meansize<-(eggiz[,7]+eggiz[,8]+eggiz[,9]+eggiz[,10]+eggiz[,11]+eggiz[,12]+eggiz[,13]+eggiz[,14]+eggiz[,15]+eggiz[,16])/10
eggiz$meansizeMm<-eggiz$meansize/eggiz$scale*10

#add latitude info
eggiz<-eggiz %>%
  mutate(latitude=case_when(
    Area=="Kristineberg" ~ Coord[Coord$Area=="Kristineberg",]$Lat,
    Area=="Arendal" ~ Coord[Coord$Area=="Arendal",]$Lat,
    Area=="Austevoll" ~ Coord[Coord$Area=="Austevoll",]$Lat,
    Area=="Hitra" ~ Coord[Coord$Area=="Hitra",]$Lat,
    Area=="Helligvaer" ~ Coord[Coord$Area=="Helligvaer",]$Lat,
    Area=="Ringstad" ~ Coord[Coord$Area=="Ringstad",]$Lat))

#add degreedays fixed info
eggiz<-eggiz %>%
  mutate(DDF=case_when(
    Area=="Kristineberg" ~ Env[Env$Area=="Kristineberg",]$Ddaysfixed ,
    Area=="Arendal" ~ Env[Env$Area=="Arendal",]$Ddaysfixed ,
    Area=="Austevoll" ~ Env[Env$Area=="Austevoll",]$Ddaysfixed ,
    Area=="Hitra" ~ Env[Env$Area=="Hitra",]$Ddaysfixed ,
    Area=="Helligvaer" ~ Env[Env$Area=="Helligvaer",]$Ddaysfixed ,
    Area=="Ringstad" ~ Env[Env$Area=="Ringstad",]$Ddaysfixed ))


str(eggiz)

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# 2. Exploration of the data prior to path analysis (and producing manuscript figures) ----
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#Color scheme for populations
popcols<-c("red3",  "brown1", "darkorchid4",  "mediumpurple1", "blue3",  "steelblue1", "grey15")



## 2.1  Exploring for outliers, cleaning up the data in general ----
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


###-----------Brood cannibalism and overgrowth---------------------
#__--__--__--__--__--__--__--__--__--__--__--__--__--__--__--__--__--__--

#These are percentage of surfaces, close to zero. What is a good distribution for such data?

#First of all, which data should we use? Live eggs confirmed during analysis of course!
nestlive<-nestdat[nestdat$Live_eggs_analysis=="y",]

# For cannibalism, we had agreed to only look at nests that are not eyed (not about to hatch)
ggplot(data=nestlive[which(nestlive$Eyed_analysis=="n"),], aes(y=area,x=PropCann*100))+
  geom_point(aes(fill=area), shape=21, size=2,color="grey90",alpha=0.5, position = position_jitterdodge(jitter.width = 2))+
  geom_boxplot(aes(fill=area),color="black", alpha=0.5, outlier.shape=NA)+
  scale_fill_manual("Population", values=popcols)+
  #scale_x_continuous()+
  scale_x_sqrt(breaks=c(0,2,5,10,20,30,40,60))+
  facet_grid(cols=vars(Timeseason))+
  labs(title="Cannibalism", y="Population", x="Cannibalised (%)")
#The distributions don't look too awful. There is a few high values outliers I want to check. Basically the few points above 15%, as it also seems biologically surprising
nestlive[which(nestlive$Eyed_analysis=="n" & nestlive$PropCann>0.15),c("area","location","Timepoint","nest_ID","comment_analysis","Eyed_analysis","PropCann")]
#After visual inspection, the following two nests are removed
nestlive<-nestlive[which(!(nestlive$area=="Arendal" & nestlive$location=="Langrumpa" & nestlive$Timepoint==2 & nestlive$nest_ID==11)),]
nestlive<-nestlive[which(!(nestlive$area=="Arendal" & nestlive$location=="Langrumpa" & nestlive$Timepoint==2 & nestlive$nest_ID==12)),]
# Make sure cannibalism is not considered outside of its range when doing multivariate models
nestanalysis<-nestdat
nestanalysis[!(rownames(nestanalysis) %in% rownames(nestlive[which(nestlive$Eyed_analysis=="n"),])),]$PropCann<-NA


#Now looking at overgrowth
ggplot(data=nestlive[which(nestlive$Eyed_analysis=="n"),], aes(y=area,x=PropGrown*100))+
  geom_point(aes(fill=area), shape=21, size=2,color="grey90",alpha=0.5, position = position_jitterdodge(jitter.width = 2))+
  geom_boxplot(aes(fill=area),color="black", alpha=0.5, outlier.shape=NA)+
  scale_fill_manual("Population", values=popcols)+
  #scale_x_continuous()+
  scale_x_sqrt(breaks=c(0,2,5,10,20,30,40,60))+
  facet_grid(cols=vars(Timeseason))+
  labs(title="Overgrowth", y="Population", x="Overgrwon (%)")
# I see one super weird point in ringstad T3. It looks legit actually
nestlive[which(nestlive$area=="Ringstad" & nestlive$PropGrown>0.02),c("area","location","Timepoint","nest_ID","comment_analysis","Eyed_analysis","PropCann")]
# Make sure overgrowth is not considered outside of its range when doing multivariate models
length(nestlive[which(nestlive$Eyed_analysis=="n"),]$PropGrown)
nestlive[which(nestlive$Eyed_analysis=="n"),]$PropGrown
nestanalysis[!(rownames(nestanalysis) %in% rownames(nestlive[which(nestlive$Eyed_analysis=="n"),])),]$PropGrown<-NA
length(nestanalysis[!is.na(nestanalysis$PropGrown),]$PropGrown)

###-----------Male length and conditioon data---------------------
#__--__--__--__--__--__--__--__--__--__--__--__--__--__--__--__--__--__--
ggplot(data=nestanalysis, aes(y=area,x=length))+
  geom_point(aes(fill=area), shape=21, size=2,color="grey90",alpha=0.5, position = position_jitterdodge(jitter.width = 2))+
  geom_boxplot(aes(fill=area),color="black", alpha=0.5, outlier.shape=NA)+
  scale_fill_manual("Population", values=popcols)+
  #scale_x_continuous()+
  facet_grid(cols=vars(Timeseason))+
  labs(title="Male length (mm)", y="Population", x="Male length (mm)")
#In Helligvaer there is a bimodal distribution, but we know that these are the two year classes. Nothing odd here

ggplot(data=nestanalysis, aes(y=area,x=con_all))+
  geom_point(aes(fill=area), shape=21, size=2,color="grey90",alpha=0.5, position = position_jitterdodge(jitter.width = 2))+
  geom_boxplot(aes(fill=area),color="black", alpha=0.5, outlier.shape=NA)+
  scale_fill_manual("Population", values=popcols)+
  #scale_x_continuous()+
  facet_grid(cols=vars(Timeseason))+
  labs(title="Male length (mm)", y="Population", x="Male length (mm)")
#Nothing very suspicious here

###-----------Eggs (reproductive success)---------------------
#__--__--__--__--__--__--__--__--__--__--__--__--__--__--__--__--__--__--
#can we include the few nests with "eyed" eggs that have an estimate of egg number? are they outliers?
ggplot(data=nestlive[which(nestlive$Eyed_analysis=="n"),], aes(y=area,x=Eggs))+
  geom_point(aes(fill=area), shape=21, size=2,color="grey90",alpha=0.5, position = position_jitterdodge(jitter.width = 2))+
  geom_point(data=nestlive[which(nestlive$Eyed_analysis=="y"),],aes(y=area,x=Eggs, fill=area),color="black", shape=21, size=5,alpha=0.5, position = position_jitterdodge(jitter.width = 2))+
  geom_boxplot(aes(fill=area),color="black", alpha=0.5, outlier.shape=NA)+
  scale_fill_manual("Population", values=popcols)+
  #scale_x_continuous()+
  facet_grid(cols=vars(Timeseason), rows=vars(Nest_size))+
  labs(title="Estimated RS", y="Population", x="Eggs")
#They look fine. No clear outliers we can include everyone.

#Okay the data is sorter, wer can use:
str(nestanalysis)



## 2.2  Risks of collinearity, VIF, corrplots, PCA ----
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
nestcor<-nestanalysis[,c("area","Timeseason","surfcm","length","DDF","PropCann","PropGrown","Eggs","con_all")]
str(nestcor)
#Peak season
corrplot(cor(nestcor[nestcor$Timeseason=="Peak season",c(3:9)], use="complete.obs"), 
         type="upper",col = COL2('RdBu', 10),addCoef.col = 'black',tl.col="grey40", diag=F, outline=F)
#Late season
corrplot(cor(nestcor[nestcor$Timeseason=="Late season",c(3:9)], use="complete.obs"), 
         type="upper",col = COL2('RdBu', 10),addCoef.col = 'black',tl.col="grey40", diag=F, outline=F)




## 2.3  get a feeling of relationships between variables ----
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&




#or not

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# 3. DETERMINANTS of REPRODUCTIVE SUCCESS ACROSS LATITUDE ----
#Path analysis combining within and across population effects
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# We want to understand if and how the environment is driving reproductive success in the nests
# We expect some latitudinal effects through variation on climate,
#but also some within population effects (e.g. male-male competition)
# We propose a global framework to investigate these two angles together
#For some variables, it will be important to disentangle the across population trends from the within population effects
#this will be interesting for example for male body size, where the mean population size could have an effect on RS across pops
# while the relative size within pop can also affect RS
#
# To disentengle these two effects in a single model, we can decompose the trait male body length into
# (i) population mean
# (ii) individual deviation from population mean
# In theory this could also be done at the level of the locality but we may lack sample size in some locations.
#Question: can we add a random population effect too? This should be needed to get real within population effects, since latitude etc don't fully account for among population differences

# The variables concerned by within population effects are: male traits (length and condition), and nest status (cannibalism and overgrowth)
#For these variables I will calculate the population mean and integrate it in the dataset, then calculate the deviation of each individual from the population mean.
#I also want to replace Latitude by a more biologically meaningful variable representing the environment, I will use DDF from the MEPS paper.
# Before decomposing the male traits and nest variables, I should probably scale them.

## DATA FORMATTING----------------------------------------------------------------------------------
#__---__---__---__---__---__---__---__---__---__---__---__---__---__---__---__---__---__---__---__---

# A. Create a scaled version of the variables that will be decomposed
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
nestanalysis$length_sc<-scale(nestanalysis$length)
nestanalysis$con_sc<-scale(nestanalysis$con_all)
# It does not make sense to scale the Cannibalism and Overgrowth as they are true proportions. I think
hist(nestanalysis$con_sc)

# B. Decompose the variables that have both across and within population effects
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#First I need to format the data to add the mean per population of the variables of interest. Let's start with length 
poplengthz<-as.data.frame(tapply(nestanalysis$length, list(area=nestanalysis$area, Timepoint=nestanalysis$Timepoint), mean, na.rm=T))
nestanalysis$poplength<-NA
for (i in levels(nestanalysis$area)) {for (j in (2:3)) {nestanalysis[nestanalysis$area==i & nestanalysis$Timepoint==j,]$poplength<-poplengthz[i,as.character(j)]}}
ggplot(nestanalysis,aes(x=area:factor(Timepoint),y=poplength))+geom_point(position=position_jitter()) #sanity check
#then I had one column for the within population deviation from the mean
nestanalysis$withlength<-nestanalysis$length-nestanalysis$poplength
ggplot(nestanalysis,aes(x=area:factor(Timepoint),y=withlength))+geom_point(position=position_jitter()) #sanity check

#Scaled version of male bodylength
poplengthz_sc<-as.data.frame(tapply(nestanalysis$length_sc, list(area=nestanalysis$area, Timepoint=nestanalysis$Timepoint), mean, na.rm=T))
nestanalysis$poplength_sc<-NA
for (i in levels(nestanalysis$area)) {for (j in (2:3)) {nestanalysis[nestanalysis$area==i & nestanalysis$Timepoint==j,]$poplength_sc<-poplengthz_sc[i,as.character(j)]}}
ggplot(nestanalysis,aes(x=area:factor(Timepoint),y=poplength_sc))+geom_point(position=position_jitter()) #sanity check
#then I had one column for the within population deviation from the mean
nestanalysis$withlength_sc<-nestanalysis$length_sc-nestanalysis$poplength_sc
ggplot(nestanalysis,aes(x=area:factor(Timepoint),y=withlength_sc))+geom_point(position=position_jitter()) #sanity check

#Does it matter if we scale before or after decomposition?
ggplot(nestanalysis,aes(x=scale(withlength), y=withlength_sc))+geom_point()
hist(scale(nestanalysis$withlength))
#Nope

#now male condition
popconz<-as.data.frame(tapply(nestanalysis$con_all, list(area=nestanalysis$area, Timepoint=nestanalysis$Timepoint), mean, na.rm=T))
popconz
nestanalysis$popcon<-NA
for (i in levels(nestanalysis$area)) {for (j in (2:3)) {nestanalysis[nestanalysis$area==i & nestanalysis$Timepoint==j,]$popcon<-popconz[i,as.character(j)]}}
ggplot(nestanalysis,aes(x=area:factor(Timepoint),y=popcon))+geom_point(position=position_jitter()) #sanity check
#then I had one column for the within population deviation from the mean
nestanalysis$withcon<-nestanalysis$con_all-nestanalysis$popcon
ggplot(nestanalysis,aes(x=area:factor(Timepoint),y=withcon))+geom_point(aes(color=Nest_size),position=position_jitter()) #sanity check

#cannibalism
popCanz<-as.data.frame(tapply(nestanalysis$PropCann, list(area=nestanalysis$area, Timepoint=nestanalysis$Timepoint), mean, na.rm=T))
popCanz
nestanalysis$popCan<-NA
for (i in levels(nestanalysis$area)) {for (j in (2:3)) {nestanalysis[nestanalysis$area==i & nestanalysis$Timepoint==j,]$popCan<-popCanz[i,as.character(j)]}}
ggplot(nestanalysis,aes(x=area:factor(Timepoint),y=popCan))+geom_point(position=position_jitter()) #sanity check
#then I had one column for the within population deviation from the mean
nestanalysis$withCan<-nestanalysis$PropCann-nestanalysis$popCan
ggplot(nestanalysis,aes(x=area:factor(Timepoint),y=withCan))+geom_point(aes(color=Nest_size),position=position_jitter()) #sanity check

#overgrowth
popGroz<-as.data.frame(tapply(nestanalysis$PropGrown, list(area=nestanalysis$area, Timepoint=nestanalysis$Timepoint), mean, na.rm=T))
popGroz
nestanalysis$popGro<-NA
for (i in levels(nestanalysis$area)) {for (j in (2:3)) {nestanalysis[nestanalysis$area==i & nestanalysis$Timepoint==j,]$popGro<-popGroz[i,as.character(j)]}}
ggplot(nestanalysis,aes(x=area:factor(Timepoint),y=popGro))+geom_point(position=position_jitter()) #sanity check
#then I had one column for the within population deviation from the mean
nestanalysis$withGro<-nestanalysis$PropGrown-nestanalysis$popGro
ggplot(nestanalysis,aes(x=area:factor(Timepoint),y=withGro))+geom_point(aes(color=Nest_size),position=position_jitter()) #sanity check

#Overgrwoth with offset 
nestanalysis$PropGoff<-(nestanalysis$PropGrown*196+0.5)/197
popGrofz<-as.data.frame(tapply(nestanalysis$PropGoff, list(area=nestanalysis$area, Timepoint=nestanalysis$Timepoint), mean, na.rm=T))
popGrofz
nestanalysis$popGroff<-NA
for (i in levels(nestanalysis$area)) {for (j in (2:3)) {nestanalysis[nestanalysis$area==i & nestanalysis$Timepoint==j,]$popGroff<-popGrofz[i,as.character(j)]}}
ggplot(nestanalysis,aes(x=area:factor(Timepoint),y=popGroff))+geom_point(position=position_jitter()) #sanity check
#then I had one column for the within population deviation from the mean
nestanalysis$withGroff<-nestanalysis$PropGoff-nestanalysis$popGroff
ggplot(nestanalysis,aes(x=area:factor(Timepoint),y=withGroff))+geom_point(aes(color=Nest_size),position=position_jitter()) #sanity check


## STATISTICAL MODELS-----------------------
#__---__---__---__---__---__---__---__---__---__---__---__---__---__---__---__---__---__---__---__---
#These 6 models are needed to describe the path diagram.
#We run them once for peak season and once for late season

##PEAK SEASON----
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

###Model 1: Reproductive success explained by EVRYTHING ----------------- or almost everything (collinearuty issues) ----
#scale variables to get meaningful effect sizes
hist(nestanalysis$Eggs)
momod1.2<-lmer(scale(Eggs)~scale(DDF)+scale(poplength)+scale(popcon)+popCan+popGro+relevel(Nest_size, ref="Small nest")+scale(withlength)+scale(withcon)+withCan+withGro+(1|location), data=nestanalysis[nestanalysis$Timeseason=="Peak season",])
vif(momod1.2, type="predictor")
summary(momod1.2)
Anova(momod1.2, type=2)
#dropping pop level condition because of collinearity with DDF
momod1.2.2<-lmer(scale(Eggs)~scale(DDF)+scale(poplength)+popCan+popGro+Nest_size+scale(withlength)+scale(withcon)+withCan+withGro+(1|location), data=nestanalysis[nestanalysis$Timeseason=="Peak season",])
vif(momod1.2.2, type="predictor")
#now we have an acceptable model (outcome is unchanged since neither DDF or condition affected RS
summary(momod1.2.2)
Anova(momod1.2.2, type=2)

#visualise significant effects to make sure I interpret effect sizes properly
plot(effect("withCan",momod1.2.2))
plot(effect("withlength",momod1.2.2))
plot(effect("Nest_size",momod1.2.2))
#Adjusted R square
r.squaredGLMM(momod1.2.2) #very similar R squared value to the model without random effect here
#I find the effect size of cannibalism on RS a bit strong. Could be due to outliers driving the stuff
ggplot(nestanalysis, aes(x=withCan, y= Eggs))+geom_point()+facet_wrap(~Nest_size:Timeseason)
#Ok looks legit, remember this is within pop. What is driving cannibalism within pop? not known


###Model 2: Cannibalism is the response---------------------------------------------
#Problem: how to handle these percentage (true proportion ) data?
ggplot (nestmale, aes(x=PropCann))+geom_histogram(binwidth=0.005)+labs(x="Proportion surface canniblaised")
#let\s have a closer look at the data first> do we have a lot of true zeroes_ what about the very high ouliers like 60% cannibalism ? check
#looking at the plots, I see a few true zeroes in cannibalism, and a lot in overgrowth, so we should probably include zero inflation
#there might be a few extreme upwards outliers to investigate, in cannibalism and overgrowth

#betareg... maybe not ideal because I have zero inflation in my data? I actually don't have a lot of true zeros... there is always a few eggs missing
#problem with betareg: no random effect

#data transformation that is sort of not sneered upon by statisticians
length(nestanalysis[!is.na(nestanalysis$PropCann),]$PropGrown)
nestanalysis$PropCaoff<-(nestanalysis$PropCann*193+0.5)/194
beta2 <- betareg(PropCaoff~scale(DDF)+scale(poplength)+scale(popcon)+Nest_size+scale(withlength)+scale(withcon), data = nestanalysis[nestanalysis$Timeseason=="Peak season",])
#for cannibalism, there are very few zeroes and when can question whether they are true zeroes so I am happy with the offset betareg except that it lacks random structure
vif(beta2, type="predictor") #as usual remove condition at the pop level as it is collinear to DDF
beta2.2 <- betareg(PropCaoff~scale(DDF)+scale(poplength)+Nest_size+scale(withlength)+scale(withcon), data = nestanalysis[nestanalysis$Timeseason=="Peak season",])
vif(beta2.2, type="predictor") #NOICE
Anova(beta2.2, type=2)
summary(beta2.2)

#What about glmmTMB? still with offset data but includes the random effect
Tmbeta3<-glmmTMB(PropCaoff~scale(DDF)+scale(poplength)+Nest_size+scale(withlength)+scale(withcon)+(1|location), 
                 data = nestanalysis[nestanalysis$Timeseason=="Peak season",],
                 family = beta_family(link = "logit") )
summary(Tmbeta3)
Anova(Tmbeta3)


#use a linear model to calculate effect size to report on the path diagram
momod2<-lmer((PropCann)~scale(DDF)+scale(poplength)+Nest_size+scale(withlength)+scale(withcon)+(1|area), data=nestanalysis[nestanalysis$Timeseason=="Peak season",])
#Boundary fit warning: random effects are very small but the model converged. Use area instead of location!
summary(momod2)
plot(effect("DDF",momod2))



###Model 3: Overgrowth is the response----------------------------------------------
#here there are a lot of true zeroes in the north, we NEED zero inflation model I think + still a problem with no random structure
#In Bert's paper they mention a common transformation to deal with the zeroes in betareg. Let's have a look
hist(nestanalysis$PropGrown, breaks=150)
length(nestanalysis[!is.na(nestanalysis$PropGrown),]$PropGrown)
ggplot(nestanalysis, aes(x=PropGrown, y=(PropGrown*196+0.5)/197))+geom_point()

###old way to deal with beta distributed stuff
#it's for sure better than offsetting the point, at least the order is respected
#see what gives with betareg
beta3.2 <- betareg(PropGoff~scale(DDF)+scale(poplength)+scale(popcon)+relevel(Nest_size, ref="Small nest")+scale(withlength)+scale(withcon), data = nestanalysis[nestanalysis$Timeseason=="Peak season",])
vif(beta3.2, type="predictor")
#again collinearity issue with DDF and condition. remove condition
beta3.4 <- betareg(PropGoff~scale(DDF)+scale(poplength)+relevel(Nest_size, ref="Small nest")+scale(withlength)+scale(withcon), data = nestanalysis[nestanalysis$Timeseason=="Peak season",])
vif(beta3.4, type="predictor") #BEAUTIFUL
summary(beta3.4)  #SHIIIIIIIIT it looks so different with the collinearity out of the way
Anova(beta3.4, type=2)

#There are a lot of zeroes in the growth data so we use a hurdle model (There is not such thing as a zero inflated betregression)
#binomial model to account for zeroes
#1. create the binomial variable
nestanalysis$GrownBin<-NA
nestanalysis[which(nestanalysis$PropGrown>0),]$GrownBin<-"y"
nestanalysis[which(nestanalysis$PropGrown==0),]$GrownBin<-"n"
table(nestanalysis$GrownBin)
#visual check
ggplot(nestanalysis[!is.na(nestanalysis$GrownBin),],aes(x=area, fill=GrownBin))+geom_bar()
#2. run the binomial model
hurdlbin2<-glmer(factor(GrownBin)~scale(DDF)+scale(poplength)+relevel(Nest_size, ref="Small nest")+scale(withlength)+scale(withcon)+(1|location), 
                 data = nestanalysis[nestanalysis$Timeseason=="Peak season",],
                 family=binomial(link = "logit"))
summary(hurdlbin2)
Anova(hurdlbin2)

#betareg for strictly positive values
Tmbeta1.2<-glmmTMB(PropGoff~scale(DDF)+scale(poplength)+relevel(Nest_size, ref="Small nest")+scale(withlength)+scale(withcon)+(1|location), 
                   data = nestanalysis[nestanalysis$Timeseason=="Peak season" & nestanalysis$GrownBin=="y",],
                   family = beta_family(link = "logit") )
summary(Tmbeta1.2)
Anova(Tmbeta1.2)

# DDF and nest size affect how much growth there is. Not Huge difference from the simple offset model. But cool :)

#lm for effect sizes
momod3<-lmer(PropGrown~scale(DDF)+scale(poplength)+Nest_size+scale(withlength)+scale(withcon)+(1|location), data=nestanalysis[nestanalysis$Timeseason=="Peak season",])
summary(momod3)
Anova(momod3)
plot(effect("DDF",momod3))
plot(effect("Nest_size",momod3))
plot(effect("scale(poplength)",momod3))





###model 4: Nest size is the response ----
momod4<-glmer(family=binomial(link = "logit"),relevel(Nest_size, ref="Small nest")~scale(poplength)+scale(popcon)+scale(withlength)+scale(withcon)+(1|area), data=nestanalysis[nestanalysis$Timeseason=="Peak season",])
#singular fit? estimation of random effect failed. Compare to non mixed model?
vif(momod4, type="predictor") #no problem
summary(momod4)
Anova(momod4)
r.squaredGLMM(momod4)

momod4.2<-glm(family=binomial(link = "logit"),relevel(Nest_size, ref="Small nest")~scale(poplength)+scale(popcon)+scale(withlength)+scale(withcon), data=nestanalysis[nestanalysis$Timeseason=="Peak season",])
#singular fit? estimation of random effect failed. Compare to non mixed model?
summary(momod4.2)
Anova(momod4.2)
r.squaredGLMM(momod4.2)
#how to explain the effect sizes with binomial model?
plot(effect("scale(withlength)",momod4.2))
plot(momod4.2$fitted.values~(nestanalysis[nestanalysis$Timeseason=="Peak season" & !is.na(nestanalysis$withlength),]$withlength))

plot(momod4.2$fitted.values/(1-momod4.2$fitted.values)~nestanalysis[nestanalysis$Timeseason=="Peak season" & !is.na(nestanalysis$withlength),]$withlength)

plot(momod4.2$fitted.values~(nestanalysis[nestanalysis$Timeseason=="Peak season" & !is.na(nestanalysis$withlength),]$withcon))


###model 5: DDF on Length ----
momod5<-lmer(scale(length)~scale(DDF)+(1|location),data=nestanalysis[nestanalysis$Timeseason=="Peak season",] )
summary(momod5)
Anova(momod5, type=2)
r.squaredGLMM(momod5)

###model 6: DDF on condition ----
momod6<-lmer(scale(con_all)~scale(DDF)+(1|location),data=nestanalysis[nestanalysis$Timeseason=="Peak season",] )
summary(momod6)
Anova(momod6, type=2)
r.squaredGLMM(momod6)


##Late SEASON----
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

###Model 1: Reproductive success explained by EVRYTHING ----------------- or almost everything (collinearuty issues)----
###scale variables to get meaningful effect sizes
momod1.2<-lmer(scale(Eggs)~scale(DDF)+scale(poplength)+scale(popcon)+popCan+popGro+relevel(Nest_size, ref="Small nest")+scale(withlength)+scale(withcon)+withCan+withGro+(1|location), data=nestanalysis[nestanalysis$Timeseason=="Late season",])
vif(momod1.2, type="predictor")
#This time we have a problem with population level growth we might have to drop it. Collinear with what? probably DDF or nest size
momod1.2.2<-lmer(scale(Eggs)~scale(DDF)+scale(poplength)+scale(popcon)+popCan+Nest_size+scale(withlength)+scale(withcon)+withCan+withGro+(1|location), data=nestanalysis[nestanalysis$Timeseason=="Late season",])
vif(momod1.2.2, type="predictor")
#acceptable
summary(momod1.2.2)
Anova(momod1.2.2, type=2)

#visualise significant effects to make sure I interpret effect sizes properly
plot(effect("withCan",momod1.2.2))
plot(effect("withlength",momod1.2.2))
plot(effect("Nest_size",momod1.2.2))
#Adjusted R square
r.squaredGLMM(momod1.2.2) #very similar R squared value to the model without random effect here
#I find the effect size of cannibalism on RS a bit strong. Could be due to outliers driving the stuff
ggplot(nestanalysis, aes(x=withCan, y= Eggs))+geom_point()+facet_wrap(~Nest_size:Timeseason)
#Ok looks legit, remember this is within pop. What is driving cannibalism within pop? not known


###Model 2: Cannibalism is the response---------------------------------------------
#Problem: how to handle these percentage (true proportion ) data?
ggplot (nestmale, aes(x=PropCann))+geom_histogram(binwidth=0.005)+labs(x="Proportion surface canniblaised")
#let\s have a closer look at the data first> do we have a lot of true zeroes_ what about the very high ouliers like 60% cannibalism ? check
#looking at the plots, I see a few true zeroes in cannibalism, and a lot in overgrowth, so we should probably include zero inflation
#there might be a few extreme upwards outliers to investigate, in cannibalism and overgrowth

#betareg... maybe not ideal because I have zero inflation in my data? I actually don't have a lot of true zeros... there is always a few eggs missing
#problem with betareg: no random effect

#data transformation that is sort of not sneered upon by statisticians
length(nestanalysis[!is.na(nestanalysis$PropCann),]$PropGrown)
nestanalysis$PropCaoff<-(nestanalysis$PropCann*193+0.5)/194
beta2 <- betareg(PropCaoff~scale(DDF)+scale(poplength)+scale(popcon)+Nest_size+scale(withlength)+scale(withcon), data = nestanalysis[nestanalysis$Timeseason=="Late season",])
#for cannibalism, there are very few zeroes and when can question whether they are true zeroes so I am happy with the offset betareg except that it lacks random structure
vif(beta2, type="predictor") #it's fine
Anova(beta2, type=2)
summary(beta2)

#What about glmmTMB? still with offset data but includes the random effect
Tmbeta3<-glmmTMB(PropCaoff~scale(DDF)+scale(poplength)+scale(popcon)+Nest_size+scale(withlength)+scale(withcon)+(1|location), 
                 data = nestanalysis[nestanalysis$Timeseason=="Late season",],
                 family = beta_family(link = "logit") )
summary(Tmbeta3)
Anova(Tmbeta3)


#use a linear model to calculate effect size to report on the path diagram
momod2<-lmer((PropCann)~scale(DDF)+scale(poplength)+scale(popcon)+Nest_size+scale(withlength)+scale(withcon)+(1|area), data=nestanalysis[nestanalysis$Timeseason=="Late season",])
#Boundary fit warning: random effects are very small but the model converged. Use area instead of location!
summary(momod2)
plot(effect("DDF",momod2))






###Model 3: Overgrowth is the response----------------------------------------------
#here there are a lot of true zeroes in the north, we NEED zero inflation model I think + still a problem with no random structure
#In Bert's paper they mention a common transformation to deal with the zeroes in betareg. Let's have a look
hist(nestanalysis$PropGrown, breaks=150)
length(nestanalysis[!is.na(nestanalysis$PropGrown),]$PropGrown)
ggplot(nestanalysis, aes(x=PropGrown, y=(PropGrown*196+0.5)/197))+geom_point()

#it's for sure better than offsetting,or?
#see what gives with betareg
beta3.2 <- betareg(PropGoff~scale(DDF)+scale(poplength)+scale(popcon)+relevel(Nest_size, ref="Small nest")+scale(withlength)+scale(withcon), data = nestanalysis[nestanalysis$Timeseason=="Late season",])
vif(beta3.2, type="predictor") #it's fine
summary(beta3.2)  
Anova(beta3.2, type=2)

#What about glmmTMB? still with offset data but includes the random effect
Tmbeta2<-glmmTMB(PropGoff~scale(DDF)+scale(poplength)+scale(popcon)+relevel(Nest_size, ref="Small nest")+scale(withlength)+scale(withcon)+(1|location), 
        data = nestanalysis[nestanalysis$Timeseason=="Late season",],
        family = beta_family(link = "logit") )
summary(Tmbeta2)
Anova(Tmbeta2)


#what about a hurdle model to really account for the zeroes?
#binomial model to account for zeroes
#1. create the binomial variable
nestanalysis$GrownBin<-NA
nestanalysis[which(nestanalysis$PropGrown>0),]$GrownBin<-"y"
nestanalysis[which(nestanalysis$PropGrown==0),]$GrownBin<-"n"
table(nestanalysis$GrownBin)

ggplot(nestanalysis[!is.na(nestanalysis$GrownBin),],aes(x=area, fill=GrownBin))+geom_bar()

#2. run the binomial model
hurdlbin<-glmer(factor(GrownBin)~scale(DDF)+scale(poplength)+scale(popcon)+relevel(Nest_size, ref="Small nest")+scale(withlength)+scale(withcon)+(1|location), 
                data = nestanalysis[nestanalysis$Timeseason=="Late season",],
                family=binomial(link = "logit"))
summary(hurdlbin)
Anova(hurdlbin)
# relative male lengths affects if there is growth or not
#betareg for strictly positive values
Tmbeta2.2<-glmmTMB(PropGoff~scale(DDF)+scale(poplength)+scale(popcon)+relevel(Nest_size, ref="Small nest")+scale(withlength)+scale(withcon)+(1|location), 
                 data = nestanalysis[nestanalysis$Timeseason=="Late season" & nestanalysis$GrownBin=="y",],
                 family = beta_family(link = "logit") )
summary(Tmbeta2.2)
Anova(Tmbeta2)

# DDF and nest size affect how much growth there is. Not Huge difference from the simple offset model. But cool :)

momod3<-lmer(PropGrown~scale(DDF)+scale(poplength)+scale(popcon)+relevel(Nest_size, ref="Small nest")+scale(withlength)+scale(withcon)+(1|area), data=nestanalysis[nestanalysis$Timeseason=="Late season",])
summary(momod3)
Anova(momod3)
plot(effect("DDF",momod3))
plot(effect("Nest_size",momod3))
plot(effect("scale(poplength)",momod3))



###model 4: Nest size is the response ----
momod4<-glmer(family=binomial(link = "logit"),relevel(Nest_size, ref="Small nest")~scale(poplength)+scale(popcon)+scale(withlength)+scale(withcon)+(1|area), data=nestanalysis[nestanalysis$Timeseason=="Late season",])
#singular fit? estimation of random effect failed. Compare to non mixed model?
vif(momod4, type="predictor") #no problem
summary(momod4)
Anova(momod4)
r.squaredGLMM(momod4)


###model 5: DDF on Length ----
momod5<-lmer(scale(length)~scale(DDF)+(1|location),data=nestmale[nestmale$Timeseason=="Late season",] )
summary(momod5)
Anova(momod5, type=2)
r.squaredGLMM(momod5)

###model 6: DDF on condition ----
momod6<-lmer(scale(con_all)~scale(DDF)+(1|location),data=nestanalysis[nestmale$Timeseason=="Late season",] )
summary(momod6)
Anova(momod6, type=2)
plot(effect("DDF",momod6))
r.squaredGLMM(momod6)



#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# 4. Supporting figures and models for the mansucript, in addition to path analysis ----
#___---___---___---___---___---___---___---___---___---___---___----__
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


## Figure 2 ----
# DDF and latitude? plus temperature curves...
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&#first panel with relationship between temperature and population
minidat<-data.frame(cbind(latitude=unique(nestdat$latitude),DDF=unique(nestdat$DDF), area=as.character(unique(nestdat$area))))
minidat$latitude<-as.numeric(minidat$latitude)
minidat$DDF<-as.numeric(minidat$DDF)
minidat$area<- factor(minidat$area, levels=c('Kristineberg', 'Arendal', 'Austevoll', 'Hitra', 'Helligvaer', 'Ringstad'))
ddfz<-ggplot(data=minidat, aes(x=latitude,y=DDF))+
  geom_point(size=5, aes(fill=factor(area)), shape=21)+
  scale_fill_manual("Population",values=popcols)+
  theme( text = element_text(size=15),plot.margin = margin(0.2,0.2,0.2,1, "cm"))+
    labs(title="b.",x=expression("Latitude " ( degree*N)),y="Thermal index (DDF)")
ddfz
#NOTE: for this next plot, the temperature data must be loaded with teh "final_script_REPRO.R" file from the reproductive strategy manucript
Envz<-ggplot(data=ALLmonthly, aes(x=Monthdate, y=Temperature, color=Area))+
  scale_x_date(date_labels = "%B",breaks="month")+
  geom_line(data=ALLfit[ALLfit$Date>"2022-04-01" & ALLfit$Date<"2022-08-31" & ALLfit$Area!="Northern Limit",],
            aes(y=Temperature,x=Date),linewidth=3, alpha=0.55)+
  theme(axis.text.x = element_text(angle=22.5, hjust = 1), text = element_text(size=15),plot.margin = margin(0.2,0.2,0.2,1, "cm"))+
  scale_colour_manual("Population",values =popcols)+
  labs(title="a.",x="Month", y=expression("Sea surface temeprature " ( degree*C)))
Envz
ggarrange(Envz,ddfz, common.legend = T, legend = "right")
ggsave("figures/Figure2.pdf", width=11, height=4.5)

## Figure 3 ----
#Nest occupancy figure.
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
nestsamp<-as.data.frame(table(nestdat$area, nestdat$Timeseason))
colnames(nestsamp)<-c("area","Timeseason","N")
nestsamp$text<-"N="
nestsamp

ggplot(data=nestdat[!(nestdat$area=="Hitra" & nestdat$Timepoint==2),], aes(y=area))+
  facet_wrap(~Timeseason)+
  #facet_wrap(~Nest_size:Timeseason)+
  #facet_grid(cols=vars(Timeseason))+
  scale_fill_manual("Nest status",
                    labels = c("Empty","Male only", "Eggs only",  "Male & Eggs"),
                    values =c("grey45","darkorchid4","orange1","tomato2"))+
  geom_bar(aes(y=area, fill=Live_eggs_analysis:male),position="fill")+
  theme( axis.text = element_text( size = 14 ),
         axis.text.x = element_text( size = 12,angle=45, hjust = 1 ),
         axis.title = element_text( size = 15 ),
         strip.text = element_text(size = 15),
         legend.text = element_text(size=12),
         legend.title = element_text(size=14),
         panel.background = element_rect(fill = 'white', colour = 'white') )+
  geom_bar(data=nestdat[!(nestdat$location=="Feoya f4-5" & nestdat$Timepoint==2) &
                          nestdat$area=="Hitra" & nestdat$Timepoint==2,],aes(y=area),fill="grey70",alpha=0.9,position="fill")+
    geom_text(data=nestsamp,aes(y=area,x=0.90,label=paste(text,N)), size=5, color="white")+
  labs(title="", x="Proportion of nests", y="Population")
ggsave("figures/Figure3.pdf", width=11, height=4.5)





#model for nest occupancy with DDF
#interactions where not significant so removed from the final model
modoccup1<-glmer(data=nestdat, Live_eggs_analysis~scale(DDF)+Timeseason+Nest_size+(1|location), family="binomial")
summary(modoccup1)
Anova(modoccup1, type=2)
#translate logit into human language
#intercept (peak season)
exp(0.297596)/(1+exp(0.297596))
#SE
exp(0.297596)/(1+exp(0.297596))-exp(0.297596-0.175920)/(1+exp(0.297596-0.175920))
exp(0.297596+0.175920)/(1+exp(0.297596+0.175920))-exp(0.297596)/(1+exp(0.297596))
#late season
exp(0.297596-0.658348)/(1+exp(0.297596-0.658348))
#SE
exp(0.297596-0.658348)/(1+exp(0.297596-0.658348))-exp(0.297596-0.658348-0.132121)/(1+exp(0.297596-0.658348-0.132121))
#slope for DDF
exp(0.0009127     )/(1+exp(0.0009127     ))


plot(effect("DDF",modoccup1 ), multiline=T)
plot(effect("Timeseason",modoccup1 ), multiline=T)
#This next line is to get the model fit into a nice data frame you can use in the plots
fittime<-as.data.frame(effect("DDF",modoccup1, xlevels=20))
fittime[20,1]-fittime[1,1]
fittime[20,2]-fittime[1,2] # increase of probability to occupy nest
(fittime[20,2]-fittime[1,2])/fittime[1,2]
minidat

## Effect cannibalism on RS // Suppl figure X ----
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
nestcols<-c("royalblue1","orange2")

ggplot(nestanalysis, aes(x=withCan, y=Eggs,color=Nest_size))+
  geom_point(aes(color=Nest_size), size=2, alpha=0.6)+
  geom_line(aes(), stat="smooth",method = "lm", size =1.5,alpha = 0.4)+
  facet_grid(cols=vars(Timeseason))+
  theme( text = element_text(size=15))+
  scale_color_manual("Nest size", values=nestcols)+
  labs(x="Cannibalism index, relative to population mean",y="Reproductive success")
ggsave("figures/Suppl_cannib.pdf", width=10, height=4.5)

median(nestanalysis[nestanalysis$Timeseason=="Peak season",]$PropCann, na.rm=T)
max(nestanalysis[nestanalysis$Timeseason=="Peak season",]$PropCann, na.rm=T)
median(nestanalysis[nestanalysis$Timeseason=="Late season",]$PropCann, na.rm=T)
max(nestanalysis[nestanalysis$Timeseason=="Late season",]$PropCann, na.rm=T)

## Fig 5. Effect of male length and nest size on RS // Figure 5 ----
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

peakz<-ggplot(nestanalysis[nestanalysis$Timeseason=="Peak season",], aes(x=withlength, y=Eggs,color=Nest_size, group=Nest_size, fill=Nest_size))+
  geom_point(size=2, alpha=0.6)+
  #facet_grid(rows=vars(Timeseason))+
  scale_color_manual("Nest size", values=nestcols)+
  scale_fill_manual("Nest size", values=nestcols)+
  geom_line(aes(), stat="smooth",method = "lm", size =1.5,alpha = 0.4)+
  labs(x="Relative male length",y="Reproductive success", title="a. Peak season")+
  theme_light()+theme(   legend.position = "none" ,text = element_text(size=15),plot.title = element_text(size=15),,panel.border = element_rect(colour = "black", fill=NA, linewidth=0.8))
peakz
 peaksM<- ggMarginal(peakz,type="boxplot", margins="x",groupColour = T,groupFill= T)
 peaksM
 latz<-ggplot(nestanalysis[nestanalysis$Timeseason=="Late season",], aes(x=withlength, y=Eggs,color=Nest_size, group=Nest_size, fill=Nest_size))+
   geom_point( size=2, alpha=0.6)+
   #facet_grid(rows=vars(Timeseason))+
   scale_color_manual("Nest size", values=nestcols)+
   scale_fill_manual("Nest size", values=nestcols)+
   geom_line(aes(), stat="smooth",method = "lm", size =1.5,alpha = 0.4)+
   labs(x="Relative male length",y="", title="b. Late season")+
   theme_light()+theme( text = element_text(size=15),plot.title = element_text(size=15),,panel.border = element_rect(colour = "black", fill=NA, linewidth=0.8))
 latz
 latzM<- ggMarginal(latz,type="boxplot", margins="x",groupColour = T,groupFill= T)
 latzM
 
 ggarrange(peaksM, latzM, widths = c(0.75,1))
 ggsave("figures/Figure 5.pdf", width=10, height=4.5)
 
 #further analyse the relationship between male relative length and reproductive success?
 #is it okay now to fit a model with interaction with nest size? Post hoc model
 
 postRS<-lmer(data=nestanalysis,Eggs~Nest_size+withlength+Timeseason+
                Nest_size:withlength+
                +(1|location))
 vif(postRS)
 Anova(postRS,type=3)
 plot(effect("Nest_size:scale(withlength)",postRS))
 plot(effect("scale(withlength)",postRS2))
 plot(effect("Nest_size",postRS2))
 #compare the slopes of male length and nest size with zero
 m.lst <- lstrends(postRS,~Nest_size|Timeseason, var="withlength",adjust = "sidak")
 summary(m.lst, infer=c(TRUE,TRUE),null=0)
 
 
## Fig 6 and 7 Effect of DDF on male traits and nest status variables // Figure 6 & 7 ----
 #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 #This next line is to get the model fit into a nice data frame you can use in the plots

 
 #model 5: DDF on Length
 modlenglate<-lmer(length~DDF+(1|location),data=nestmale[nestmale$Timeseason=="Late season",] )
 fitlenglate<-as.data.frame(effect("DDF",modlenglate,x.var="DDF", xlevels=20))
 modlengpeak<-lmer(length~DDF+(1|location),data=nestmale[nestmale$Timeseason=="Peak season",] )
 fitlengpeak<-as.data.frame(effect("DDF",modlengpeak,x.var="DDF", xlevels=20))
 
 modconlate<-lmer(con_all~DDF+(1|location),data=nestmale[nestmale$Timeseason=="Late season",] )
 fitconlate<-as.data.frame(effect("DDF",modconlate,x.var="DDF", xlevels=20))
 modconpeak<-lmer(con_all~DDF+(1|location),data=nestmale[nestmale$Timeseason=="Peak season",] )
 fitconpeak<-as.data.frame(effect("DDF",modconpeak,x.var="DDF", xlevels=20))
 
  beta3.2 <- betareg(PropGoff~DDF+scale(poplength)+scale(popcon)+Nest_size+scale(withlength)+scale(withcon), data = nestanalysis[nestanalysis$Timeseason=="Late season",])
 fitgrolate<-as.data.frame(effect("DDF",beta3.2,x.var="DDF", xlevels=20))
 beta3.2.2 <- betareg(PropGoff~DDF+scale(poplength)+scale(popcon)+Nest_size+scale(withlength)+scale(withcon), data = nestanalysis[nestanalysis$Timeseason=="Peak season",])
 fitgropeak<-as.data.frame(effect("DDF",beta3.2.2,x.var="DDF", xlevels=20))
 
 beta2 <- betareg(PropCaoff~DDF+scale(poplength)+scale(popcon)+Nest_size+scale(withlength)+scale(withcon), data = nestanalysis[nestanalysis$Timeseason=="Late season",])
 fitcanlate<-as.data.frame(effect("DDF",beta2,x.var="DDF", xlevels=20))
 beta2.2 <- betareg(PropCaoff~DDF+scale(poplength)+scale(popcon)+Nest_size+scale(withlength)+scale(withcon), data = nestanalysis[nestanalysis$Timeseason=="Peak season",])
 fitcanpeak<-as.data.frame(effect("DDF",beta2.2,x.var="DDF", xlevels=20))

 
 

lengz<-ggplot(nestanalysis, aes(x=DDF))+
  scale_fill_manual("Populations",values=popcols)+
  scale_color_manual("Time period", values=c("black","grey80"))+
  scale_linewidth_manual("Time period",values=c(1,0.8))+
  geom_boxplot(alpha=0.5, aes( group=area:Timeseason, fill=area,, y=length, color=Timeseason, linewidth = Timeseason))+
  geom_point(alpha=0.2,position=position_jitterdodge(jitter.width =80, dodge.width = 50), aes( group=area:Timeseason, fill=area, y=length))+
  #This next line does the confidence interval around the model fit
  #geom_ribbon(data=fitlengpeak, aes(x = DDF, ymin=fit-se, ymax=fit+se), alpha=0.2) + 
  #This next line does the model fit
  #geom_line(data=fitlengpeak,aes(x=DDF, y=fit), size=1.5, color="black",alpha=0.5)+
  #This next line does the confidence interval around the model fit
  #geom_ribbon(data=fitlenglate, aes(x = DDF, ymin=fit-se, ymax=fit+se), alpha=0.1) + 
  #This next line does the model fit
  #geom_line(data=fitlenglate,aes(x=DDF, y=fit), size=1.5, color="grey50",alpha=0.5)+
  theme_light()+theme( text = element_text(size=15),plot.title = element_text(size=15),,panel.border = element_rect(colour = "black", fill=NA, linewidth=0.8))+
  #facet_wrap(~Timeseason)+
  labs(title="a. ",x="Thermal environment (DDF)", y="Male body length (mm)")
lengz

condz<-ggplot(nestanalysis, aes(x=DDF))+
  scale_fill_manual("Populations",values=popcols)+
  scale_color_manual("Time period", values=c("black","grey80"))+
  scale_linewidth_manual("Time period",values=c(1.,0.8))+
  geom_boxplot(alpha=0.5, outlier.shape=NA,aes( group=area:Timeseason, fill=area,, y=con_all, color=Timeseason, linewidth = Timeseason))+
  geom_point(alpha=0.2,position=position_jitterdodge(jitter.width =80, dodge.width = 50), aes( group=area:Timeseason, fill=area, y=con_all))+
  #This next line does the confidence interval around the model fit
  #geom_ribbon(data=fitconpeak, aes(x = DDF, ymin=fit-se, ymax=fit+se), alpha=0.2) + 
  #This next line does the model fit
  #geom_line(data=fitconpeak,aes(x=DDF, y=fit), size=1.5, color="black",alpha=0.5)+
  #This next line does the confidence interval around the model fit
  #geom_ribbon(data=fitconlate, aes(x = DDF, ymin=fit-se, ymax=fit+se), alpha=0.1) + 
  #This next line does the model fit
  #geom_line(data=fitconlate,aes(x=DDF, y=fit), size=1.5, color="grey50",alpha=0.5)+
  theme_light()+theme( text = element_text(size=15),plot.title = element_text(size=15),panel.border = element_rect(colour = "black", fill=NA, linewidth=0.8))+
  #facet_wrap(~Timeseason)+
  labs(title="b. ",x="", y="Male condition")
condz


Caz<-ggplot(nestanalysis, aes(x=DDF))+
  scale_fill_manual("Populations",values=popcols)+
  scale_color_manual("Time period", values=c("black","grey80"))+
  scale_linewidth_manual("Time period",values=c(1,0.8))+
  geom_boxplot(alpha=0.5, outlier.shape=NA, aes( group=area:Timeseason, fill=area,, y=PropCann*100, color=Timeseason, linewidth = Timeseason))+
  geom_point(alpha=0.2,position=position_jitterdodge(jitter.width =80, dodge.width = 50), aes( group=area:Timeseason, fill=area, y=PropCann*100))+
  #This next line does the confidence interval around the model fit
 #geom_ribbon(data=fitcanpeak, aes(x = DDF, ymin=fit-se, ymax=fit+se), alpha=0.2) + 
  #This next line does the model fit
  #geom_line(data=fitcanpeak,aes(x=DDF, y=fit), size=1.5, color="black",alpha=0.5)+
  #This next line does the confidence interval around the model fit
  #geom_ribbon(data=fitcanlate, aes(x = DDF, ymin=fit-se, ymax=fit+se), alpha=0.1) + 
  #This next line does the model fit
  #geom_line(data=fitcanlate,aes(x=DDF, y=fit), size=1.5, color="grey90",alpha=0.5)+
  theme_light()+theme( text = element_text(size=15),plot.title = element_text(size=15),panel.border = element_rect(colour = "black", fill=NA, linewidth=0.8))+
  #facet_wrap(~Timeseason)+
  labs(title="a. ",x="Thermal environment (DDF)", y="Cannibalism index")
Caz

Grz<-ggplot(nestanalysis, aes(x=DDF))+
  scale_fill_manual("Populations",values=popcols)+
  scale_color_manual("Time period", values=c("black","grey80"))+
  scale_linewidth_manual("Time period",values=c(1,0.8))+
  geom_boxplot(alpha=0.5, outlier.shape=NA, aes( group=area:Timeseason, fill=area,, y=PropGrown*100, color=Timeseason, linewidth = Timeseason))+
  geom_point(alpha=0.2,position=position_jitterdodge(jitter.width =80, dodge.width = 50), aes( group=area:Timeseason, fill=area, y=PropGrown*100))+
  #This next line does the confidence interval around the model fit
  #geom_ribbon(data=fitgropeak, aes(x = DDF, ymin=fit-se, ymax=fit+se), alpha=0.2) + 
  #This next line does the model fit
  #geom_line(data=fitgropeak,aes(x=DDF, y=fit), size=1.5, color="black",alpha=0.5)+
  #This next line does the confidence interval around the model fit
  #geom_ribbon(data=fitgrolate, aes(x = DDF, ymin=fit-se, ymax=fit+se), alpha=0.1) + 
  #This next line does the model fit
  #geom_line(data=fitgrolate,aes(x=DDF, y=fit), size=1.5, color="grey50",alpha=0.5)+
  theme_light()+theme( text = element_text(size=15),plot.title = element_text(size=15),panel.border = element_rect(colour = "black", fill=NA, linewidth=0.8))+
  #facet_wrap(~Timeseason)+
  labs(title="b. ",x="", y="Overgrowth index")
Grz



ggarrange(lengz+coord_flip()+scale_x_reverse(), condz+coord_flip()+scale_x_reverse(),ncol=2,common.legend = T, legend = "right")
ggsave("figures/Figure 6.pdf", width=10, height=6.5)
ggarrange(Caz+coord_flip()+scale_x_reverse(), Grz+coord_flip()+scale_x_reverse(),ncol=2,common.legend = T, legend = "right")
ggsave("figures/Figure 7.pdf", width=10, height=6.5)


## Fig 8 RS and egg density across latitude, egg size // Figure 8 ----
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
# RS alone could be a supplementary figure
Eggz<-ggplot(nestanalysis, aes(x=DDF))+
  scale_fill_manual("Populations",values=popcols)+
  scale_color_manual("Time period", values=c("black","grey80"))+
  scale_linewidth_manual("Time period",values=c(1,0.8))+
  geom_boxplot(alpha=0.5, outlier.shape=NA, aes( group=area:Timeseason, fill=area,, y=Eggs, color=Timeseason, linewidth = Timeseason))+
  geom_point(alpha=0.2,position=position_jitterdodge(jitter.width =80, dodge.width = 50), aes( group=area:Timeseason, fill=area, y=Eggs))+
  #This next line does the confidence interval around the model fit
  #geom_ribbon(data=fitgropeak, aes(x = DDF, ymin=fit-se, ymax=fit+se), alpha=0.2) + 
  #This next line does the model fit
  #geom_line(data=fitgropeak,aes(x=DDF, y=fit), size=1.5, color="black",alpha=0.5)+
  #This next line does the confidence interval around the model fit
  #geom_ribbon(data=fitgrolate, aes(x = DDF, ymin=fit-se, ymax=fit+se), alpha=0.1) + 
  #This next line does the model fit
  #geom_line(data=fitgrolate,aes(x=DDF, y=fit), size=1.5, color="grey50",alpha=0.5)+
  theme_light()+theme( text = element_text(size=15),plot.title = element_text(size=15),panel.border = element_rect(colour = "black", fill=NA, linewidth=0.8))+
  facet_wrap(~Nest_size)+
  labs(title=" ",x="Thermal environment (DDF)", y="Eggs")+coord_flip()+scale_x_reverse()
Eggz
ggsave("figures/SupplS4.pdf", width=10, height=6.5)


Denz<-ggplot(nestanalysis, aes(x=DDF))+
  scale_fill_manual("Populations",values=popcols)+
  scale_color_manual("Time period", values=c("black","grey80"))+
  scale_linewidth_manual("Time period",values=c(1,0.8))+
  geom_hline(yintercept=175, color="grey40")+
  geom_boxplot(alpha=0.5, outlier.shape=NA, aes( group=area:Timeseason, fill=area,, y=Eggs/BroodCm, color=Timeseason, linewidth = Timeseason))+
  geom_point(alpha=0.2,position=position_jitterdodge(jitter.width =80, dodge.width = 50), aes( group=area:Timeseason, fill=area, y=Eggs/BroodCm))+
  #This next line does the confidence interval around the model fit
  #geom_ribbon(data=fitgropeak, aes(x = DDF, ymin=fit-se, ymax=fit+se), alpha=0.2) + 
  #This next line does the model fit
  #geom_line(data=fitgropeak,aes(x=DDF, y=fit), size=1.5, color="black",alpha=0.5)+
  #This next line does the confidence interval around the model fit
  #geom_ribbon(data=fitgrolate, aes(x = DDF, ymin=fit-se, ymax=fit+se), alpha=0.1) + 
  #This next line does the model fit
  #geom_line(data=fitgrolate,aes(x=DDF, y=fit), size=1.5, color="grey50",alpha=0.5)+
  theme_light()+theme( text = element_text(size=15),plot.title = element_text(size=15),panel.border = element_rect(colour = "black", fill=NA, linewidth=0.8))+
  #facet_wrap(~Nest_size)+
  labs(title="a. ",x="Thermal environment (DDF)", y=bquote('Egg density'~(Eggs/Cm^2)))+coord_flip()+scale_x_reverse()
Denz

Egsiz<-ggplot(eggiz, aes(x=DDF))+
  scale_fill_manual("Populations",values=popcols)+
  scale_color_manual("Time period", values=c("black","grey80"))+
  scale_linewidth_manual("Time period",values=c(1,0.8))+
  geom_hline(yintercept=0.68, color="grey40")+
  geom_boxplot(alpha=0.5, outlier.shape=NA, aes( group=Area:Timeseason, fill=Area,, y=meansizeMm  , color=Timeseason, linewidth = Timeseason))+
  geom_point(alpha=0.2,position=position_jitterdodge(jitter.width =80, dodge.width = 50), aes( group=Area:Timeseason, fill=Area, y=meansizeMm ))+
  #This next line does the confidence interval around the model fit
  #geom_ribbon(data=fitgropeak, aes(x = DDF, ymin=fit-se, ymax=fit+se), alpha=0.2) + 
  #This next line does the model fit
  #geom_line(data=fitgropeak,aes(x=DDF, y=fit), size=1.5, color="black",alpha=0.5)+
  #This next line does the confidence interval around the model fit
  #geom_ribbon(data=fitgrolate, aes(x = DDF, ymin=fit-se, ymax=fit+se), alpha=0.1) + 
  #This next line does the model fit
  #geom_line(data=fitgrolate,aes(x=DDF, y=fit), size=1.5, color="grey50",alpha=0.5)+
  theme_light()+theme( text = element_text(size=15),plot.title = element_text(size=15),panel.border = element_rect(colour = "black", fill=NA, linewidth=0.8))+
  #facet_wrap(~Nest_size)+
  labs(title="b. ",x="", y="Egg diamemter (mm)")+coord_flip()+scale_x_reverse()
Egsiz

ggarrange(Denz,Egsiz,ncol=2,common.legend = T, legend = "right")
ggsave("figures/Figure 8.pdf", width=11, height=6.5)

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
##Statistical models for egg density and egg size ----
#egg size
modegg<-lmer(meansizeMm~DDF+Nest_size_pic+Timeseason+
             (1|Location), data=eggiz)
Anova(modegg,type=2)
summary(modegg)
plot(effect("DDF",modegg ), multiline=T)
plot(effect("Timeseason",modegg ), multiline=T)
mean(eggiz$meansizeMm)
emmeans(modegg, list(pairwise ~ Timeseason), adjust = "tukey")
(mean(eggiz[eggiz$Area=="Ringstad",]$meansizeMm)-mean(eggiz[eggiz$Area=="Kristineberg",]$meansizeMm))/mean(eggiz[eggiz$Area=="Kristineberg",]$meansizeMm)*100

#egg density
hist(nestanalysis$Eggs/nestanalysis$BroodCm)
modens<-lmer(Eggs/BroodCm~DDF+Nest_size+Timeseason+
               (1|location), data=nestanalysis)
#no signif interactions
Anova(modens,type=2)
summary(modens)
plot(effect("DDF",modens ), multiline=T)
plot(effect("Nest_size",modens ), multiline=T)
#DDF effect size?
emmeans(modens, list(pairwise ~ Nest_size), adjust = "tukey")
nestanalysis$density<-nestanalysis$Eggs/nestanalysis$BroodCm
(mean(nestanalysis[nestanalysis$area=="Ringstad",]$density, na.rm=T )/mean(nestanalysis[nestanalysis$area=="Kristineberg",]$density, na.rm=T )-1)*100


#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
## Supplementary figure for overgrowth affected by nest size and male length ----
str(nestanalysis)

#nest size and overgrowth
grovio<-ggplot(nestanalysis, aes(y=PropGrown*100, x=Nest_size))+
  geom_violin(alpha=0.3)+
  geom_point(aes(color=area), position = position_jitterdodge())+
  scale_color_manual("Population", values=popcols)+
  labs(title="a.",x="Nest size",y="Overgrowth index %")

#peak season population length
grolen<-ggplot(nestanalysis[nestanalysis$Timeseason=="Peak season",], aes(y=PropGrown*100, x=length))+
  geom_point(aes(color=area), size=2, alpha=0.6)+
  scale_color_manual("Population", values=popcols)+
  facet_wrap(~Timeseason)+
  labs(title="b.",x="Male length (mm)",y="Overgrowth index %")

#Late season, overgrowth is affected by the realtive length of guarding male?
growi<-ggplot(nestanalysis[nestanalysis$Timeseason=="Late season",], aes(y=PropGrown*100, x=withlength))+
  geom_point(aes(color=area), size=2, alpha=0.6)+
  scale_color_manual("Population", values=popcols)+
  facet_wrap(~Timeseason)+
  labs(title="c.",x="Male relative length",y="Overgrowth index %")

ggarrange(grovio,grolen,growi, ncol=3, common.legend = T, legend = "right")

  







# OLD STUFF ----

#-----------------------------------
#LATE SEASON
#-----------
#Model 1: Reproductive success explained by EVRYTHING
#scale variables to get meaningful effect sizes
momod1.3<-lmer(scale(Eggs)~scale(DDF)+scale(poplength)+scale(popcon)+popCan+popGro+relevel(Nest_size, ref="Small nest")+scale(withlength)+scale(withcon)+withCan+withGro+(1|location), data=nestmale[nestmale$Timeseason=="Late season",])
summary(momod1.3)
Anova(momod1.3, type=2)
#visualise significant effects to make sure I interpret effect sizes properly
plot(effect("withCan",momod1.2))
plot(effect("withlength_sc",momod1.2))
plot(effect("Nest_size",momod1.2))
ggplot(nestmale,aes(x=scale(withlength), y=scale(Eggs)))+geom_point()
#Adjusted R square
r.squaredGLMM(momod1.3) #very similar R squared value to the model without random effect here



#Model 2: Cannibalism is the response
#Problem: how to handle these percentage (true proportion ) data?
ggplot (nestmale, aes(x=PropCann))+geom_histogram(binwidth=0.005)+labs(x="Proportion surface canniblaised")
#
momod2.3<-lmer((PropCann+0.001)~scale(DDF)+scale(poplength)+scale(popcon)+Nest_size+scale(withlength)+scale(withcon)+(1|area), data=nestmale[nestmale$Timeseason=="Late season",])
#Boundary fit warning: random effects are very small but the model converged. Use area instead of location!
summary(momod2.3)
Anova(momod2)
#strange that removing the random effects makes the fixed effects more significant.... not what I expected.
momod2.2<-lm((PropCann+0.001)~scale(DDF)+scale(poplength)+scale(popcon)+Nest_size+scale(withlength)+scale(withcon), data=nestmale[nestmale$Timeseason=="Late season",])
summary(momod2.2)
Anova(momod2.2)

#betareg... maybe not ideal because I have zero inflation in my data? I actually don't have a lot of true zeros... there is always a few eggs missing
#problem with betareg: no random effect
beta2.2 <- betareg(I(PropCann+0.001)~scale(DDF)+scale(poplength)+scale(popcon)+Nest_size+scale(withlength)+scale(withcon), data = nestmale[nestmale$Timeseason=="Late season",])
Anova(beta2.2, type=2)
summary(beta2.2)
#no significant effect, only poplength marginally


#Model 3: Overgrowth is the response
momod3.3<-lmer((PropGrown+0.001)~scale(DDF)+scale(poplength)+scale(popcon)+Nest_size+scale(withlength)+scale(withcon)+(1|location), data=nestmale[nestmale$Timeseason=="Late season",])
summary(momod3.3)
Anova(momod3)
r.squaredGLMM(momod3)
plot(effect("Nest_size",momod3))
plot(effect("scale(poplength)",momod3))

momod3.2<-lm(PropGrown~scale(DDF)+scale(poplength)+scale(popcon)+Nest_size+scale(withlength)+scale(withcon), data=nestmale[nestmale$Timeseason=="Late season",])
summary(momod3.2)
Anova(momod3.2)
#Warning the poplength effect disappears if I have the random effect, what does that mean? ask Yimen?


beta3.2 <- betareg(I(PropGrown+0.001)~scale(DDF)+scale(poplength)+scale(popcon)+relevel(Nest_size, ref="Small nest")+scale(withlength)+scale(withcon), data = nestmale[nestmale$Timeseason=="Late season",])
summary(beta3.2)
Anova(beta3.2, type=2)


#model 4: Nest size is the response
momod4.3<-glmer(family=binomial(link = "logit"),relevel(Nest_size, ref="Small nest")~scale(poplength)+scale(popcon)+scale(withlength)+scale(withcon)+(1|area), data=nestmale[nestmale$Timeseason=="Late season",])
#singular fit? estimation of random effect failed. Compare to non mixed model?
summary(momod4.3)
Anova(momod4.3)
r.squaredGLMM(momod4.3)

momod4.4<-glm(family=binomial(link = "logit"),relevel(Nest_size, ref="Small nest")~scale(poplength)+scale(popcon)+scale(withlength)+scale(withcon), data=nestmale[nestmale$Timeseason=="Late season",])
#singular fit? estimation of random effect failed. Compare to non mixed model?
summary(momod4.4)
Anova(momod4.2)
r.squaredGLMM(momod4.2)

#model 5: DDF on Length
momod5.2<-lmer(scale(length)~scale(DDF)+(1|location),data=nestmale[nestmale$Timeseason=="Late season",] )
summary(momod5.2)
Anova(momod5.2, type=2)
r.squaredGLMM(momod5.2)

#model 6: DDF on condition
momod6.2<-lmer(scale(con_all)~scale(DDF)+(1|location),data=nestmale[nestmale$Timeseason=="Late season",] )
summary(momod6.2)
Anova(momod6.2, type=2)
r.squaredGLMM(momod6.2)








########################################################################
# 2. Graphical exploration
#######################################################################
#Color scheme for populations
popcols<-c("red3",  "brown1", "darkorchid4",  "mediumpurple1", "blue3",  "steelblue1", "grey15")

#Nest occupancy-----------------------------------------------------------
#NOTE data from Hitra:2 Feoya should be removed as empty nests weren't collected (storm)
ggplot(data=nestdat[!(nestdat$location=="Feoya f4-5" & nestdat$Timepoint==2),], aes(x=area, fill=male:Live_eggs_analysis))+
  #facet_wrap(~Timeseason)+
  #facet_wrap(~Nest_size:Timeseason)+
  facet_grid(cols=vars(Timeseason), rows=vars(Nest_size))+
  scale_fill_manual("Nest status",
                    labels = c("Empty", "Eggs only", "Male only", "Male & Eggs"),
                    values =c("grey","orange1","darkorchid4","tomato2"))+
  geom_bar()+
  theme( axis.text = element_text( size = 14 ),
         axis.text.x = element_text( size = 12,angle=45, hjust = 1 ),
         axis.title = element_text( size = 15 ),
         strip.text = element_text(size = 15),
         legend.text = element_text(size=12),
         legend.title = element_text(size=14))+
  labs(title="", x="Population", y="Nest count")


#Nest fullness of nests that have live eggs only-------------------------------
ggplot(data=nestdat[nestdat$Live_eggs_field=="y",], aes(x=area,y=TotBroodFrac))+
  geom_point(aes(fill=area), shape=21, color="grey90",alpha=0.5, position = position_jitterdodge(jitter.width = 0.8))+
  geom_boxplot(aes(fill=area), alpha=0.5, outlier.shape=NA)+
  scale_fill_manual("Population", values=popcols)+
  theme( axis.text = element_text( size = 14 ),
         axis.text.x = element_text( size = 12,angle=45, hjust = 1 ),
         axis.title = element_text( size = 15 ),
         strip.text = element_text(size = 15),
         legend.text = element_text(size=12),
         legend.title = element_text(size=14))+
  facet_grid(cols=vars(Timeseason), rows=vars(Nest_size))+
  labs(title="", x="Population", y="Nest Fullness")
#Questions: what are the nests that have low fullness? are they hatched? are they old, new, overgown?
#are they cannibalised?, ovegrown?, hatching
nestdat[nestdat$Live_eggs_field=="y"
        & nestdat$TotBroodFrac<=0.25,]
table(nestdat$Live_eggs_field, nestdat$Live_eggs_analysis)
#it is not the fact that they are considered as dead eggs either in the filed or later looking at pictures.
ggplot(nestdat, aes(y=TotBroodFrac, x=area))+
  geom_boxplot()+ facet_grid(cols=vars(Live_eggs_field), rows=vars(Live_eggs_analysis))
#If anything, with eyed eggs the nests are a bit more full-> hatching does not explain low brood fraction
ggplot(nestdat, aes(y=TotBroodFrac, x=area))+
  geom_boxplot()+ facet_grid(cols=vars(Eyed_analysis))
#Nor high cannibalism explains it
ggplot(nestdat, aes(y=TotBroodFrac, x=CannibFrac))+geom_point()
#Nor high canniblaism
ggplot(nestdat, aes(y=TotBroodFrac, x=OvergrownFrac))+geom_point()
#need to do a visual check on these nests, here they are
na.omit(nestdat[nestdat$Live_eggs_analysis=="y"  & nestdat$TotBroodFrac<=0.25,c("area","location","Timepoint","nest_ID", "TotBroodFrac", "male")])
# Conclusion: they mostly correspond to healthy nests with only one clutch. I think they belong to the data :)


#Cannibalism as a fraction of the total brood area---------------------------------------------
ggplot(data=nestdat[nestdat$Live_eggs_field=="y",], aes(x=area, y=(CannibCm/TotBroodCm)*100))+
  geom_point(aes(fill=area), shape=21, color="grey90",alpha=0.5, position = position_jitterdodge(jitter.width = 0.8))+
  geom_boxplot(aes(fill=area), alpha=0.5, outlier.shape=NA)+
  scale_fill_manual("Population", values=popcols)+
  theme( axis.text = element_text( size = 14 ),
         axis.text.x = element_text( size = 12,angle=45, hjust = 1 ),
         axis.title = element_text( size = 13 ),
         strip.text = element_text(size = 15),
         legend.text = element_text(size=12),
         legend.title = element_text(size=14))+
  facet_grid(cols=vars(Timeseason), rows=vars(Nest_size))+
  scale_y_continuous(trans='sqrt')+
  labs(title="", x="Population", y="Brood cannibalism (% of brood surface)")

hist(nestdat[nestdat$Live_eggs_field=="y",]$CannibCm/nestdat[nestdat$Live_eggs_field=="y",]$TotBroodCm)

#Overgrowth as a fraction of the total brood area----------------------------------------------
ggplot(data=nestdat, aes(x=area, y=(OvergrownCm/TotBroodCm)*100))+
  geom_point(aes(fill=area), shape=21, color="grey90",alpha=0.5, position = position_jitterdodge(jitter.width = 0.8))+
  geom_boxplot(aes(fill=area), alpha=0.5, outlier.shape=NA)+
  scale_fill_manual("Population", values=popcols)+
  theme( axis.text = element_text( size = 14 ),
         axis.text.x = element_text( size = 12,angle=45, hjust = 1 ),
         axis.title = element_text( size = 13 ),
         strip.text = element_text(size = 15),
         legend.text = element_text(size=12),
         legend.title = element_text(size=14))+
  facet_grid(cols=vars(Timeseason), rows=vars(Nest_size))+
  scale_y_continuous(trans='sqrt')+
  labs(title="", x="Population", y="Overgrowth (% of brood surface)")


#Number of eggs/reproductive success ----------------------------------------------
#This considers only the brood area that is not overgrown or cannibalised
ggplot(data=nestdat, aes(x=area, y=Eggs, group=area:Nest_size))+
  geom_point(aes(fill=area), shape=21, color="grey90",alpha=0.5, position = position_jitterdodge(jitter.width = 0.8))+
  geom_boxplot(aes(fill=area, linewidth = Nest_size), alpha=0.5, outlier.shape=NA)+
  scale_fill_manual("Population", values=popcols)+
  scale_linewidth_manual("Nest size",values=c(0.8,0.1))+
  theme( axis.text = element_text( size = 14 ),
         axis.text.x = element_text( size = 12,angle=45, hjust = 1 ),
         axis.title = element_text( size = 15 ),
         strip.text = element_text(size = 15),
         legend.text = element_text(size=12),
         legend.title = element_text(size=14))+
  facet_wrap(~Timeseason)+
  #facet_grid(cols=vars(Timeseason), rows=vars(Nest_size))+
 # scale_y_continuous(trans='log')+
  labs(title="", x="Population", y="Egg number")

#sanity check Hitra has higher RS than anyone else, is that a true pattern?
#Hitra also has large eggs like other northern populations
#Brood surface? higher than southern pops but not crazy high, mostly, no low values..TotBroodCm
#Nest surface in number of squares? very average-> different grid pattern did not cause the effect
#Cannib is quite low
#density is not that high, but egg number is..yes becasue the main driver of the Hitra success is brood area.
#variable TotbroodFrac is very high in Hitra ()
#nestdat$TotBroodFrac<-nestdat$TotBrood_surface/(nestdat$X1sq_surface*nestdat$nestsurf)
#Totbrood surface is mostly responsible for the pattenr

#I did visul checks, and nests from Hitra look very full from edge to edge.

#Egg density------------------------------------------------------------------
ggplot(data=nestdat, aes(x=area, y=Eggs/BroodCm, group=area:Nest_size))+
  geom_point(aes(fill=area), shape=21, color="grey90",alpha=0.5, position = position_jitterdodge(jitter.width = 0.8))+
  geom_boxplot(aes(fill=area, linewidth = Nest_size), alpha=0.5, outlier.shape=NA)+
  scale_fill_manual("Population", values=popcols)+
  scale_linewidth_manual("Nest size",values=c(0.8,0.1))+
  theme( axis.text = element_text( size = 14 ),
         axis.text.x = element_text( size = 12,angle=45, hjust = 1 ),
         axis.title = element_text( size = 15 ),
         strip.text = element_text(size = 15),
         legend.text = element_text(size=12),
         legend.title = element_text(size=14))+
  facet_wrap(~Timeseason)+
  #facet_grid(cols=vars(Timeseason), rows=vars(Nest_size))+
  # scale_y_continuous(trans='log')+
  labs(title="", x="Population", y="Egg density (egg/cm^2)")


#Show the complex relationship between male length, nest size, and reproductive success
#reproductive success as eggs counted 
ggplot(nestdat, aes(x=length, y=Eggs))+
  geom_point(aes(fill=Nest_size), shape=21)+
  geom_line(aes( group=Nest_size, color=Nest_size), stat="smooth",method = "lm", size = 1.5,alpha = 0.4)+
  geom_line( stat="smooth",method = "lm", size = 2.,alpha = 0.3)+
  #ggpubr::stat_regline_equation(label.y = c(12000,5000), label.x=c(40,60), aes(group=Nest_size,label = ..rr.label..))+
  facet_grid(cols=vars(Timeseason))+
  ylim(0,15000)+
  labs(x="Male total length (mm)",y="Eggs in the nest")

ggplot(nestdat, aes(x=length, y=Eggs))+
  geom_point(aes(fill=Nest_size), shape=21)+
  geom_line(aes( group=Nest_size, color=Nest_size), stat="smooth",method = "lm", size = 1.5,alpha = 0.4)+
  geom_line( stat="smooth",method = "lm", size = 2.,alpha = 0.3)+
  #ggpubr::stat_regline_equation(label.y = c(12000,5000), label.x=c(40,60), aes(group=Nest_size,label = ..rr.label..))+
  facet_grid(cols=vars(area))+
  ylim(0,15000)+
  labs(x="Male total length (mm)",y="Eggs in the nest")

#reproductive success as maximum brood surface
ggplot(nestdat, aes(x=length, y=TotBroodCm))+
  geom_point(aes(fill=Nest_size), shape=21)+
  geom_line(aes( group=Nest_size, color=Nest_size), stat="smooth",method = "lm", size = 1.5,alpha = 0.4)+
  geom_line( stat="smooth",method = "lm", size = 2.,alpha = 0.3)+
  #ggpubr::stat_regline_equation(label.y = c(12000,5000), label.x=c(40,60), aes(group=Nest_size,label = ..rr.label..))+
  facet_grid(cols=vars(Timeseason))+
  ylim(0,100)

#overgrown area
ggplot(nestdat, aes(x=100*weight/length^3, y=(OvergrownCm/TotBroodCm)*100))+
  geom_point(aes(fill=Nest_size), shape=21)+
  geom_line(aes( group=Nest_size, color=Nest_size), stat="smooth",method = "lm", size = 1.5,alpha = 0.4)+
  geom_line( stat="smooth",method = "lm", size = 2.,alpha = 0.3)+
  #ggpubr::stat_regline_equation(label.y = c(12000,5000), label.x=c(40,60), aes(group=Nest_size,label = ..rr.label..))+
  facet_grid(cols=vars(area))

#cannibalised area vs condition and length
ggplot(nestdat, aes(x=100*weight/length^3, y=(CannibCm/TotBroodCm)*100))+
  geom_point(aes(fill=Nest_size), shape=21)+
  geom_line(aes( group=Nest_size, color=Nest_size), stat="smooth",method = "lm", size = 1.5,alpha = 0.4)+
  geom_line( stat="smooth",method = "lm", size = 2.,alpha = 0.3)+
  #ggpubr::stat_regline_equation(label.y = c(12000,5000), label.x=c(40,60), aes(group=Nest_size,label = ..rr.label..))+
  facet_grid(cols=vars(area))

ggplot(nestdat, aes(x=length, y=(CannibCm/TotBroodCm)*100))+
  geom_point(aes(fill=Nest_size), shape=21)+
  geom_line(aes( group=Nest_size, color=Nest_size), stat="smooth",method = "lm", size = 1.5,alpha = 0.4)+
  geom_line( stat="smooth",method = "lm", size = 2.,alpha = 0.3)+
  #ggpubr::stat_regline_equation(label.y = c(12000,5000), label.x=c(40,60), aes(group=Nest_size,label = ..rr.label..))+
  facet_grid(cols=vars(area))


#I see that mid season larger males get larger nests which also gets them more eggs, but not larger broods?
#Overgrowth is driven by latitude only, but apparent size effect due to larger males in the North


# Looking for correlations with male traits!------------------------------------
#consider doing this per time point? Per population?

#first make a smaller dataset with only the interesting variables plus a couple of new ones
nestdat$CannibProp<-(nestdat$CannibCm/nestdat$TotBroodCm)*100
nestdat$GrowthProp<-(nestdat$OvergrownCm/nestdat$TotBroodCm)*100
nestdat$fultonC<-(nestdat$weight/nestdat$length^3)*100
nestmales<-nestdat[nestdat$male=="y",c("area", "Timeseason","surfcm","Eggs","TotBroodFrac","length","fultonC","CannibProp","GrowthProp")]


corrplot(cor(nestmales[nestmales$Timeseason=="Mid-season",3:9], use="complete.obs"), 
         type="upper",col = COL2('RdBu', 10),addCoef.col = 'black')

corrplot(cor(nestmales[nestmales$Timeseason=="Late season",3:9], use="complete.obs"), 
         type="upper",col = COL2('RdBu', 10),addCoef.col = 'black')



cor(nestmales, use="complete.obs")

summary(nestdat[nestdat$male=="y",])


nestdat[!is.na(nestdat$OvergrownCm) & nestdat$OvergrownCm/nestdat$TotBroodCm>0.2,]


nestdat[!is.na(nestdat$TotBroodFrac) & nestdat$TotBroodFrac>1,]
summary(nestdat$TotBroodFrac)
str(nestdat)

########################################################
# Male size and condition
########################################################
str(nestdat)

ggplot(nestdat, aes(x=length))+
  facet_grid( rows=vars(area))+
  geom_histogram( color="grey", linewidth=1, alpha=0.2,position="identity",
                  aes(fill=Nest_size,group=Nest_size))



#%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&
#%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&
#    ~~~~~~    Statistical analysis and publication figures         ~~~~~~          #
#####################################################################################
#%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&
#%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&
#%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&


#%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&
#1 .  Nest Occupancy------------------------------------
#########################################################

#First, we are interested in the proportion of nests with live eggs. A nest with live eggs also likely has a male
#even if he is not observed on that sampling occasion

#Live egg model----------------------------------------------
#interactions where not significant so removed from the final model
modoccup1<-glm(data=nestdat, Live_eggs_analysis~latitude+Timeseason+Nest_size, family="binomial")
summary(modoccup1)
Anova(modoccup1, type=2)

plot(effect("latitude",modoccup1 ), multiline=T)
plot(effect("Timeseason",modoccup1 ), multiline=T)

modoccup1.2<-glm(data=nestdat, Live_eggs_analysis~DDF+Timeseason+Nest_size, family="binomial")
summary(modoccup1.2)
Anova(modoccup1.2, type=2)



# Figure 2 Nest occupancy------------------------------------
levels(nestdat$Live_eggs_analysis)


compact<-as.data.frame(table(nestdat[!(nestdat$location=="Feoya f4-5" & nestdat$Timepoint==2),c("area","Timeseason","Live_eggs_analysis","male")]))
compact2<-cbind(compact[1:12,1:2],"Empty"=compact[1:12,"Freq"])
compact2

ggplot(data=nestdat[!(nestdat$location=="Feoya f4-5" & nestdat$Timepoint==2),], aes(y=area, fill=Live_eggs_analysis:male))+
  #facet_wrap(~Timeseason)+
  #facet_wrap(~Nest_size:Timeseason)+
  facet_grid(cols=vars(Timeseason))+
  scale_fill_manual("Nest status",
                    labels = c("Empty","Male only", "Eggs only",  "Male & Eggs"),
                    values =c("grey","darkorchid4","orange1","tomato2"))+
 geom_bar(aes(y=area, fill=Live_eggs_analysis:male),position="fill")+
  theme( axis.text = element_text( size = 14 ),
         axis.text.x = element_text( size = 12,angle=45, hjust = 1 ),
         axis.title = element_text( size = 15 ),
         strip.text = element_text(size = 15),
         legend.text = element_text(size=12),
         legend.title = element_text(size=14))+
  geom_bar(data=nestdat[!(nestdat$location=="Feoya f4-5" & nestdat$Timepoint==2) &
                     nestdat$area=="Hitra" & nestdat$Timepoint==2,],aes(y=area),fill="grey95",alpha=0.3,position="fill")+
  labs(title="", x="Proportion of nests", y="Population")



#%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&
#2. Nest Fullness-----------------------------------------
#########################################################

#In this case there is no linear gradient with latitude so we may be better off with a population effect.
#All interactions were non-sign so removed from the final model.
modfull1<-lm(data=nestdat, TotBroodFrac~area+Timeseason+Nest_size
            # +area:Timeseason+ area:Nest_size+Timeseason:Nest_size
            )
summary(modfull1)
Anova(modfull1, type=2)


plot(effect("area",modfull1 ), multiline=T)
plot(effect("Timeseason",modfull1 ), multiline=T)
posthoc<-emmeans(modfull1, list(pairwise ~ area), adjust = "tukey")
data.frame(posthoc$`pairwise differences of area`)
write.table(format(data.frame(posthoc$`pairwise differences of area`), digits=1,scientific=F), "posthoc_fullness.txt", sep=",", row.names = F)

# Figure 3 nest fullness----------------------
ggplot(data=nestdat[nestdat$Live_eggs_field=="y",], aes(y=area,x=TotBroodFrac*100))+
  geom_point(aes(fill=area), shape=21, size=2,color="grey90",alpha=0.5, position = position_jitterdodge(jitter.width = 2))+
  geom_boxplot(aes(fill=area),color="black", alpha=0.5, outlier.shape=NA)+
  scale_fill_manual("Population", values=popcols)+
  theme( axis.text = element_text( size = 14 ),
         axis.text.x = element_text( size = 11,angle=45, hjust = 1 ),
         axis.title = element_text( size = 15 ),
         strip.text = element_text(size = 15),
         legend.position="none")+
  facet_grid(cols=vars(Timeseason))+
  labs(title="", y="Population", x="Nest Fullness (%)")


#%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&
#3. Brood cannibalism and overgrowth---------------------
#########################################################

#These are percentage of surfaces, close to zero. What is a good distribution for such data?
#Ideally I need a zero inflated beta reg
#To keep things simple, I will use betareg and offsetting the zeros.
#For cannibalism, there are few true zeroes so that is really fine.

#Whether to use population or latitude as a predictor is unclear here.
#BIC preferes latitude as it is a simpler model, but it does not capture the entirety of differences 
#between populations. At the same time, the latitude is what we are interested in.

#Model with latitude instead and see which one has lowest AIC?
betacab2 <- betareg(I(PropCann+0.001)~latitude+Timeseason+Nest_size
                   #+latitude:Timeseason+latitude:Nest_size+Timeseason:Nest_size
                     #+latitude:Timeseason:Nest_size
                     , data = nestdat)
summary(betacab2)
Anova(betacab2, type=2)
plot(effect("latitude",betacab2 ), multiline=T)


#Overgrowth: 1 sign interaction
betagro2 <- betareg(I(PropGrown+0.01)~latitude+Timeseason+Nest_size
                      +latitude:Nest_size
                     # +latitude:Timeseason:Nest_size
                    , data = nestdat)
Anova(betagro2, type=3)
plot(effect("latitude:Nest_size",betagro2 ), multiline=T)
plot(effect("Nest_size",betagro2 ), multiline=T)
plot(effect("latitude",betagro2 ), multiline=T)
AIC(betagro,betagro2)
BIC(betagro,betagro2)

# Figure 4 overgrowth and cannibalism----------------------
gro<-ggplot(data=nestdat[nestdat$Live_eggs_field=="y",], aes(y=area,x=(OvergrownCm/TotBroodCm)*100))+
  geom_point(aes(fill=area), shape=21, size=2,color="grey90",alpha=0.5, position = position_jitterdodge(jitter.width = 2))+
  geom_boxplot(aes(fill=area),color="black", alpha=0.5, outlier.shape=NA)+
  scale_fill_manual("Population", values=popcols)+
  theme( axis.text = element_text( size = 14 ),
         axis.text.x = element_text( size = 11,angle=45, hjust = 1 ),
         axis.title = element_text( size = 15 ),
         strip.text = element_text(size = 15),
         legend.position="none",
         axis.title.y=element_blank(),
         axis.text.y=element_blank())+
  scale_x_sqrt(breaks=c(0,1,2,5,10,20,30))+
  #facet_grid(cols=vars(Nest_size))+
  labs(title="b. Overgrowth", y="Population", x="Overgrown (%)")
gro
cani<-ggplot(data=nestdat[nestdat$Live_eggs_field=="y",], aes(y=area,x=PropCann*100))+
  geom_point(aes(fill=area), shape=21, size=2,color="grey90",alpha=0.5, position = position_jitterdodge(jitter.width = 2))+
  geom_boxplot(aes(fill=area),color="black", alpha=0.5, outlier.shape=NA)+
  scale_fill_manual("Population", values=popcols)+
  theme( axis.text = element_text( size = 14 ),
         axis.text.x = element_text( size = 11,angle=45, hjust = 1 ),
         axis.title = element_text( size = 15 ),
         strip.text = element_text(size = 15),
         legend.position="none",
         axis.title.y=element_blank(),
         axis.text.y=element_blank())+
         scale_y_discrete(position = "left")+
  #scale_x_continuous()+
  scale_x_sqrt(breaks=c(0,2,5,10,20,30,40,60))+
 # facet_grid(cols=vars(Timeseason))+
  labs(title="a. Cannibalism", y="Population", x="Cannibalised (%)")
cani
#how does it look with logit transofrmation?
cani2<-ggplot(data=nestdat[nestdat$Live_eggs_field=="y",], aes(y=area,x=logit(PropCann+0.001)))+
  geom_point(aes(fill=area), shape=21, size=2,color="grey90",alpha=0.5, position = position_jitterdodge(jitter.width = 2))+
  geom_boxplot(aes(fill=area),color="black", alpha=0.5, outlier.shape=NA)+
  scale_fill_manual("Population", values=popcols)+
  theme( axis.text = element_text( size = 14 ),
         axis.text.x = element_text( size = 11,angle=45, hjust = 1 ),
         axis.title = element_text( size = 15 ),
         strip.text = element_text(size = 15),
         legend.position="none",
         axis.title.y=element_blank(),
         axis.text.y=element_blank())+
  scale_y_discrete(position = "left")+
  #scale_x_continuous()+
 # scale_x_sqrt(breaks=c(0,2,5,10,20,30,40,60))+
  # facet_grid(cols=vars(Timeseason))+
  labs(title="a. Cannibalism", y="Population", x="Cannibalised (%)")
cani2
#How does it look in relation to DDF?
ggplot(data=nestdat[nestdat$Live_eggs_field=="y",], aes(x=DDF,y=logit(PropCann+0.001)))+
  geom_point(aes(fill=area), shape=21, size=2,color="grey90",alpha=0.5, position = position_jitter(width = 0.1))+
  #geom_boxplot(aes(fill=area),color="black", alpha=0.5, outlier.shape=NA)+
  #scale_fill_manual("Population", values=popcols)+

 # scale_y_discrete(position = "left")+
  #scale_x_continuous()+
  # scale_x_sqrt(breaks=c(0,2,5,10,20,30,40,60))+
  facet_grid(cols=vars(Timeseason))+
  labs(title="a. Cannibalism", y="logit cannibalism", x="DDF")
#ok I acutally believe that there is no DDF effect on cannibalism
ggplot(data=nestmale, aes(x=popcon,y=logit(PropCann+0.001)))+
  geom_point(aes(fill=area), shape=21, size=2,color="grey90",alpha=0.5)+
  #geom_boxplot(aes(fill=area),color="black", alpha=0.5, outlier.shape=NA)+
  #scale_fill_manual("Population", values=popcols)+
  
  # scale_y_discrete(position = "left")+
  #scale_x_continuous()+
  # scale_x_sqrt(breaks=c(0,2,5,10,20,30,40,60))+
  facet_grid(cols=vars(Timeseason))+
  labs(title="a. Cannibalism", y="logit cannibalism", x="condition (pop mean)")
#same with growth?
ggplot(data=nestmale, aes(x=DDF,y=logit(PropGrown+0.001)))+
  geom_point(aes(fill=area), shape=21, size=2,color="grey90",alpha=0.5, position = position_jitter(width = 30))+
  geom_boxplot(aes(fill=area),color="black", alpha=0.5, outlier.shape=NA)+
  #scale_fill_manual("Population", values=popcols)+
  
  # scale_y_discrete(position = "left")+
  #scale_x_continuous()+
  # scale_x_sqrt(breaks=c(0,2,5,10,20,30,40,60))+
  facet_grid(cols=vars(Timeseason))+
  labs(title="a. Cannibalism", y="logit cannibalism", x="DDF")
#there is a problem in peak season, little data in the north, why?
table(nestmale$PropGrown,nestmale$area)

#multiplot
ggarrange(cani,gro, ncol = 2,widths = c(2.,2.))


#%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&
#4. Size and condition of nest guarding males------------
#########################################################
#we will work only with nest that have a male measured
nestmale<-nestdat[!is.na(nestdat$length),]

#---------------------------------------------------------------------------------
#male total length in the different nest types, populations and time of the season.

#only 1 sign. interaction
modnestmal<-lm(data=nestmale, length~Nest_size+Timeseason+area+Nest_size:area)
summary(modnestmal)
Anova(modnestmal, type=3)
plot(effect("Nest_size:area",modnestmal ), multiline=T)
emmeans(modnestmal, list(pairwise ~ Nest_size|area), adjust = "tukey")
plot(effect("Timeseason",modnestmal ), multiline=T)
plot(effect("Nest_size",modnestmal ), multiline=T)

str(nestmale)

# Figure 5a size of nest guarding males-------------------------------------------------
lengz<-ggplot(nestmale, aes(x=length, y= area))+
    geom_vline(xintercept = 50, linewidth=1.5, alpha=0.3)+
  geom_point(aes(fill=area, group=Nest_size), shape=21, size=2,color="grey90",alpha=0.5, position = position_jitterdodge(jitter.width = 1))+
  geom_boxplot(aes(fill=area, linewidth = Nest_size),color="black", alpha=0.5, outlier.shape=NA)+
  scale_linewidth_manual("Nest size",values=c(0.8,0.1))+
  scale_fill_manual("Population", values=popcols)+
  theme( axis.text = element_text( size = 14 ),
         axis.text.x = element_text( size = 11,angle=45, hjust = 1 ),
         axis.title = element_text( size = 15 ),
         strip.text = element_text(size = 15))+
  facet_grid(cols=vars(Timeseason))+
  labs(title="a.", y="Population", x="Male size (mm)")


#--------------------------------------------------------------------------------
#Male condition
#First, we compute male condition as the residual of a size weight regression.
ggplot(nestmale, aes(x=length^3,y=weight, color=area))+geom_point()+
  geom_line(stat="smooth",method = "lm",  size = 1.5, alpha = 0.6)
#regression over all pops
#lmcondall<-lm(weight~I(length^3), data=nestmale)
lmcondall<-lm(log(weight)~log(length), data=nestmale)
plot(effect("length",lmcondall))

ggplot(nestmale, aes(x=length^3,y=weight))+geom_point()+geom_smooth(method="lm")
ggplot(nestmale, aes(x=length,y=weight))+geom_point()+geom_smooth(method="lm")
ggplot(nestmale, aes(x=log(length),y=log(weight)))+geom_point()+geom_smooth(method="lm")
#add residuals to the dataset
hist(residuals (lmcondall))
nestmale$con_all<-residuals (lmcondall)

# analyze condition data
modcond<-lm(data=nestmale, con_all~Nest_size+Timeseason+latitude
            +latitude:Timeseason)
summary(modcond)
Anova(modcond, type=3)

plot(effect("Timeseason:latitude",modcond ), multiline=T)
m.lst.cond <- lstrends(modcond, ~Timeseason, var="latitude",adjust = "sidak") #get slopes for each sex
summary(m.lst.cond, infer=c(TRUE,TRUE),null=0) #compare each slope to zero 
pairs(m.lst.cond)
emmeans(modcond, list(pairwise ~ Timeseason|area), adjust = "tukey")

plot(effect("latitude",modcond ), multiline=T)
plot(effect("Timeseason",modcond ), multiline=T)
plot(effect("Nest_size",modcond ), multiline=T)

# Figure 5b condition of nest guarding males-------------------------------------------------
condz<-ggplot(nestmale, aes(x=con_all, y= area))+
  geom_vline(xintercept = 0, linewidth=1.5, alpha=0.3)+
  geom_point(aes(fill=area, group=Nest_size), shape=21, size=2,color="grey90",alpha=0.5, position = position_jitterdodge(jitter.width = 1))+
  geom_boxplot(aes(fill=area, linewidth = Nest_size),color="black", alpha=0.5, outlier.shape=NA)+
  scale_linewidth_manual("Nest size",values=c(0.8,0.1))+
  scale_fill_manual("Population", values=popcols)+
  theme( axis.text = element_text( size = 14 ),
         axis.text.x = element_text( size = 11,angle=45, hjust = 1 ),
         axis.title = element_text( size = 15 ),
         strip.text = element_text(size = 15))+
  facet_grid(cols=vars(Timeseason))+
  labs(title="b.", y="Population", x="Male condition ")

ggarrange(lengz,condz, nrow=2, common.legend = T, legend = "right")

#%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&
#4.b Path analysis for latitude / male quality /nest status
#########################################################
#Correlation with non-std variables
# We start by a correlation matrix to understand how male size, male condition, cannibalism etc covary with egg number in the nest
# We do this separately by nest size of course, but combining all populations.
str(nestmale)
malecor<-nestmale[,c("area","Timeseason","Nest_size","length","latitude","TotBroodCm","PropCann","PropGrown","Eggs","con_all")]

#overall corrplot
corrplot(cor(malecor[,c(4,5,7:10)], use="complete.obs"), 
         type="upper",col = COL2('RdBu', 10),addCoef.col = 'black',tl.col="grey40", diag=F, outline=F)
#Per nest size
corrplot(cor(malecor[malecor$Nest_size=="Large nest",c(4,5,7:10)], use="complete.obs"), 
         type="upper",col = COL2('RdBu', 10),addCoef.col = 'black',tl.col="grey40", diag=F, outline=T)
corrplot(cor(malecor[malecor$Nest_size=="Small nest",c(4,5,7:10)], use="complete.obs"), 
         type="upper",col = COL2('RdBu', 10),addCoef.col = 'black',tl.col="grey40", diag=F, outline=T)

#Per season
corrplot(cor(malecor[malecor$Timeseason=="Peak season",c(4,5,7:10)], use="complete.obs"), 
         type="upper",col = COL2('RdBu', 10),addCoef.col = 'black',tl.col="grey40", diag=F, outline=T)
corrplot(cor(malecor[malecor$Timeseason=="Late season",c(4,5,7:10)], use="complete.obs"), 
         type="upper",col = COL2('RdBu', 10),addCoef.col = 'black',tl.col="grey40", diag=F, outline=T)

#per nest size and season
corrplot(cor(malecor[malecor$Nest_size=="Large nest"
                     & malecor$Timeseason=="Peak season"
                     ,4:10], use="complete.obs"), 
         type="upper",col = COL2('RdBu', 10),addCoef.col = 'black',tl.col="grey40", diag=F, outline=T)


corrplot(cor(malecor[malecor$Nest_size=="Small nest"
                     & malecor$Timeseason=="Mid-season",4:11], use="complete.obs"), 
         type="upper",col = COL2('RdBu', 10),addCoef.col = 'grey80')

corrplot(cor(malecor[malecor$Nest_size=="Large nest"
                     & malecor$Timeseason=="Late season",4:11], use="complete.obs"), 
         type="upper",col = COL2('RdBu', 10),addCoef.col = 'grey80')


corrplot(cor(malecor[malecor$Nest_size=="Small nest"
                     & malecor$Timeseason=="Late season",4:11], use="complete.obs"), 
         type="upper",col = COL2('RdBu', 10),addCoef.col = 'grey80')

#Path analysis. The response variables of growth and cannibalism require a betareg model to test for significance.
#We will use a lm for effect size maybe.
#has to be done by sampling time... and by nest size?

#model 1 => cannibalism
modpath1.1<-lm(PropCann*100~length+con_all+latitude, data=nestmale)
summary(modpath1.1)
plot(effect("latitude",modpath1.1))

betagcan1.1 <- betareg(I(PropCann+0.001)~length+con_all+latitude, data = nestmale)
Anova(betagcan1.1, type=2)

#model 2 => overgrowth
modpath2.1<-lm(100*PropGrown~length+con_all+latitude, data=nestmale)
summary(modpath2.1)

betagro2.1 <- betareg(I(PropGrown+0.001)~length+con_all+latitude, data = nestmale)
Anova(betagro2.1, type=2)

#model 3 => body length
modpath3<-lm(length~latitude, data=nestmale)
summary(modpath3)

#model 4 => body condition
modpath4<-lm(con_all~latitude, data=nestmale)
summary(modpath4)

plot (effect("latitude",modpath4))




#%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&
#5. Reproductive success and egg density-----------------
#########################################################

#Number of live eggs
modegg<-lm(data=nestdat, Eggs~area+Nest_size+Timeseason
           +area:Timeseason+Nest_size:area)
Anova(modegg, type=3)
plot(effect("area:Timeseason",modegg ), multiline=T)
plot(effect("area:Nest_size",modegg ), multiline=T)
plot(effect("area",modegg ), multiline=T)
plot(effect("Nest_size",modegg ), multiline=T)
emmeans(modegg, list(pairwise ~ Nest_size), adjust = "tukey")
posthoc<-emmeans(modegg, list(pairwise ~ area), adjust = "tukey")
emmeans(modegg, list(pairwise ~ Timeseason|area), adjust = "tukey")
emmeans(modegg, list(pairwise ~ Nest_size|area), adjust = "tukey")
#export post hoc table
write.table(format(data.frame(posthoc$`pairwise differences of area`), digits=1,scientific=F), "posthoc_eggs.txt", sep=",", row.names = F)



#Number of eggs/reproductive success Figure 6.a----------------------------------------------
#This considers only the brood area that is not overgrown or cannibalised
Eggz<-ggplot(data=nestdat, aes(y=area, x=Eggs, group=area:Nest_size))+
  geom_point(aes(fill=area), shape=21, color="grey90",alpha=0.5, position = position_jitterdodge(jitter.width = 0.8))+
  geom_boxplot(aes(fill=area, linewidth = Nest_size), alpha=0.5, outlier.shape=NA)+
  geom_vline(xintercept = 10000, linewidth=1.5, alpha=0.3)+
  scale_fill_manual("Population", values=popcols)+
  scale_linewidth_manual("Nest size",values=c(0.8,0.1))+
  theme( axis.text = element_text( size = 14 ),
         axis.text.x = element_text( size = 12,angle=45, hjust = 1 ),
         axis.title = element_text( size = 15 ),
         strip.text = element_text(size = 15),
         legend.text = element_text(size=12),
         legend.title = element_text(size=14))+
  facet_wrap(~Timeseason)+
  #facet_grid(cols=vars(Timeseason), rows=vars(Nest_size))+
  # scale_y_continuous(trans='log')+
  labs(title="", y="Population", x="Egg number")
Eggz
#Second version of the egg figure changing panel variable
Eggz<-ggplot(data=nestdat, aes(y=area, x=Eggs, group=area:Timeseason))+
  geom_point(aes(fill=area), shape=21, color="grey90",alpha=0.5, position = position_jitterdodge(jitter.width = 0.8))+
  geom_boxplot(aes(fill=area, linewidth = Timeseason), alpha=0.5, outlier.shape=NA)+
  geom_vline(xintercept = 10000, linewidth=1.5, alpha=0.3)+
  scale_fill_manual("Population", values=popcols)+
  scale_linewidth_manual("Sampling time",values=c(0.8,0.1))+
  theme( axis.text = element_text( size = 14 ),
         axis.text.x = element_text( size = 12,angle=45, hjust = 1 ),
         axis.title = element_text( size = 15 ),
         strip.text = element_text(size = 15),
         legend.text = element_text(size=12),
         legend.title = element_text(size=14))+
  facet_wrap(~Nest_size)+
  #facet_grid(cols=vars(Timeseason), rows=vars(Nest_size))+
  # scale_y_continuous(trans='log')+
  labs(title="a.", y="Population", x="Egg number")
Eggz

#Egg density
modens<-lm(data=nestdat, I(Eggs/BroodCm)~area+Nest_size+Timeseason
           +area:Timeseason+Timeseason:Nest_size)
Anova(modens, type=3)
plot(effect("area:Timeseason",modens ),x.var="Timeseason", multiline=T)
plot(effect("Nest_size:Timeseason",modens ), multiline=T)
plot(effect("area",modens ), multiline=T)
plot(effect("Nest_size",modens ), multiline=T)
emmeans(modens, list(pairwise ~ area), adjust = "tukey")
emmeans(modens, list(pairwise ~ Nest_size|Timeseason), adjust = "tukey")
emmeans(modens, list(pairwise ~ Timeseason|area), adjust = "tukey")

levels(nestdat$location)

#Egg density FIgure 6.b------------------------------------------------------------------
Denz<-ggplot(data=nestdat, aes(y=area, x=Eggs/BroodCm, group=area:Nest_size))+
  geom_point(aes(fill=area), shape=21, color="grey90",alpha=0.5, position = position_jitterdodge(jitter.width = 0.8))+
  geom_boxplot(aes(fill=area, linewidth = Nest_size), alpha=0.5, outlier.shape=NA)+
  scale_fill_manual("Population", values=popcols)+
  geom_vline(xintercept = 175, linewidth=1.5, alpha=0.3)+
  scale_linewidth_manual("Nest size",values=c(0.8,0.1))+
  theme( axis.text = element_text( size = 14 ),
         axis.text.x = element_text( size = 12,angle=45, hjust = 1 ),
         axis.title = element_text( size = 15 ),
         strip.text = element_text(size = 15),
         legend.text = element_text(size=12),
         legend.title = element_text(size=14))+
  facet_wrap(~Timeseason)+
  #facet_grid(cols=vars(Timeseason), rows=vars(Nest_size))+
  # scale_y_continuous(trans='log')+
  labs(title="", y="Population", x=bquote('Egg density'~(Eggs/Cm^2)))
#version with nest size on panels
#Egg density FIgure 6.b------------------------------------------------------------------
Denz<-ggplot(data=nestdat, aes(y=area, x=Eggs/BroodCm, group=area:Timeseason))+
  geom_point(aes(fill=area), shape=21, color="grey90",alpha=0.5, position = position_jitterdodge(jitter.width = 0.8))+
  geom_boxplot(aes(fill=area, linewidth = Timeseason), alpha=0.5, outlier.shape=NA)+
  scale_fill_manual("Population", values=popcols)+
  geom_vline(xintercept = 175, linewidth=1.5, alpha=0.3)+
  scale_linewidth_manual("Sampling time",values=c(0.8,0.1))+
  theme( axis.text = element_text( size = 14 ),
         axis.text.x = element_text( size = 12,angle=45, hjust = 1 ),
         axis.title = element_text( size = 15 ),
         strip.text = element_text(size = 15),
         legend.text = element_text(size=12),
         legend.title = element_text(size=14))+
  facet_wrap(~Nest_size)+
  #facet_grid(cols=vars(Timeseason), rows=vars(Nest_size))+
  # scale_y_continuous(trans='log')+
  labs(title="b.", y="Population", x=bquote('Egg density'~(Eggs/Cm^2)))

ggarrange(Eggz,Denz, nrow=2, common.legend = T, legend = "right")

#%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&
#5.b Egg size                          -----------------
#########################################################

eggiz<-read.delim("eggsize.txt", header=T)
eggiz$Area<-as.factor(eggiz$Area)
eggiz$Area<- factor(eggiz$Area, levels=c('Kristineberg', 'Arendal','Austevoll','Hitra','Helligvaer',  'Ringstad'))
eggiz$Nest_size_pic<-as.factor(eggiz$Nest_size_pic)
levels(eggiz$Nest_size_pic)<-c("Large","Small")
eggiz$eyed<-as.factor(eggiz$eyed)
levels(eggiz$eyed)<-c("No","Yes")
eggiz$Timeseason<-as.factor(eggiz$Timeseason)
eggiz$Timeseason<- factor(eggiz$Timeseason, levels=c('Peak season', 'Late season'))
eggiz$meansize<-(eggiz[,7]+eggiz[,8]+eggiz[,9]+eggiz[,10]+eggiz[,11]+eggiz[,12]+eggiz[,13]+eggiz[,14]+eggiz[,15]+eggiz[,16])/10
eggiz$meansizeMm<-eggiz$meansize/eggiz$scale*10
str(eggiz)

ggplot(eggiz, aes(y=meansize, x=eyed))+geom_boxplot()+geom_point(position = position_jitter())
#okay, beware there is an effect of dev stage... should have expected that.

modegg<-lm(meansizeMm~Area+Nest_size_pic+Timeseason
           +Area:Nest_size_pic, data=eggiz)
Anova(modegg,type=3)
plot(effect("Area",modegg ), multiline=T)
plot(effect("Timeseason",modegg ), multiline=T)
plot(effect("Area:Nest_size_pic",modegg ), multiline=T)
emmeans(modegg, list(pairwise ~ Area), adjust = "tukey")
emmeans(modegg, list(pairwise ~ Timeseason), adjust = "tukey")
emmeans(modegg, list(pairwise~Nest_size_pic|Area), adjust = "tukey")

ggplot(eggiz, aes(x=meansizeMm, y=Area, fill=Area, group=Area:Nest_size_pic))+
  scale_fill_manual('Population',values=popcols)+
# scale_shape_manual('Eyed eggs',values=c(21,24))+
  geom_vline(xintercept = 0.7, color="grey30")+
  geom_boxplot(alpha=.5,outlier.shape=NA, aes(linewidth=Nest_size_pic))+
  geom_point(position = position_jitterdodge(jitter.width = 0.5), size=2, alpha=0.5, shape=21,aes(fill=Area))+
  scale_linewidth_manual("Nest size",values=c(0.8,0.1))+
  facet_wrap(~Timeseason)+labs(y="Population",x="Egg diameter (mm)")

#simpler figure for conference talk
ggplot(eggiz, aes(x=meansizeMm, y=Area, fill=Area, group=Area))+
  scale_fill_manual('Population',values=popcols)+
  # scale_shape_manual('Eyed eggs',values=c(21,24))+
  #geom_vline(xintercept = 0.7, color="grey30")+
  geom_boxplot(alpha=.5,outlier.shape=NA)+
  geom_point(position = position_jitterdodge(jitter.width = 0.5), size=2, alpha=0.5, shape=21,aes(fill=Area))+
  #scale_linewidth_manual("Nest size",values=c(0.8,0.1))+
  theme( axis.text = element_text( size = 14 ),
         axis.text.x = element_text( size = 12,angle=45, hjust = 1 ),
         axis.title = element_text( size = 15 ),
         strip.text = element_text(size = 15),
         legend.text = element_text(size=12),
         legend.title = element_text(size=14))+
  #facet_wrap(~Timeseason)+
  labs(y="Population",x="Egg diameter (mm)")

#%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&%¤%&
#6. Determinants of reproductive success-----------------
#########################################################

#Here we are interested in sexual selection dynamics that may be common across populations
#thus we scale the variables within each population.
#One possible scaling is to divide by the mean, so that the scale becomes percentages

nestmale$RelEggs<-NA
nestmale$RelLength<-NA
nestmale$RelCann<-NA
nestmale$RelGro<-NA
nestmale$RelCond<-NA
nestmale$RelCond2<-NA
for (i in levels(nestmale$area)) {nestmale[nestmale$area==i,]$RelEggs<-nestmale[nestmale$area==i,]$Eggs/mean(nestmale[nestmale$area==i,]$Eggs, na.rm=T)}
for (i in levels(nestmale$area)) {nestmale[nestmale$area==i,]$RelLength<-nestmale[nestmale$area==i,]$length/mean(nestmale[nestmale$area==i,]$length, na.rm=T)}
for (i in levels(nestmale$area)) {nestmale[nestmale$area==i,]$RelCann<-nestmale[nestmale$area==i,]$PropCann/mean(nestmale[nestmale$area==i,]$PropCann, na.rm=T)}
for (i in levels(nestmale$area)) {nestmale[nestmale$area==i,]$RelGro<-nestmale[nestmale$area==i,]$PropGrown/mean(nestmale[nestmale$area==i,]$PropGrown, na.rm=T)}
for (i in levels(nestmale$area)) {nestmale[nestmale$area==i,]$RelCond2<-nestmale[nestmale$area==i,]$con_all/mean(nestmale[nestmale$area==i,]$con_all, na.rm=T)}
str(nestmale)
##For condition, we can estimate a regression per pop and then get the residuals from that
lmcondpop<-lm(log(weight)~log(length)*area, data=nestmale)
#add residuals to the dataset
nestmale$RelCond<-residuals (lmcondpop)


ggplot(nestmale, aes(x=area,y=RelCond))+geom_boxplot()


#we could also do the same thing per study site.
nestmale$LocRelEggs<-NA
nestmale$LocRelLength<-NA
nestmale$LocRelCann<-NA
nestmale$LocRelGro<-NA
nestmale$LocRelCond<-NA
for (i in levels(nestmale$location)) {nestmale[nestmale$location==i,]$LocRelEggs<-nestmale[nestmale$location==i,]$Eggs/mean(nestmale[nestmale$location==i,]$Eggs, na.rm=T)}
for (i in levels(nestmale$location)) {nestmale[nestmale$location==i,]$LocRelLength<-nestmale[nestmale$location==i,]$length/mean(nestmale[nestmale$location==i,]$length, na.rm=T)}
for (i in levels(nestmale$location)) {nestmale[nestmale$location==i,]$LocRelCann<-nestmale[nestmale$location==i,]$PropCann/mean(nestmale[nestmale$location==i,]$PropCann, na.rm=T)}
for (i in levels(nestmale$location)) {nestmale[nestmale$location==i,]$LocRelGro<-nestmale[nestmale$location==i,]$PropGrown/mean(nestmale[nestmale$location==i,]$PropGrown, na.rm=T)}
for (i in levels(nestmale$location)) {nestmale[nestmale$location==i,]$LocRelCond<-nestmale[nestmale$location==i,]$con_all/mean(nestmale[nestmale$location==i,]$con_all, na.rm=T)}


# now we are ready for path analysis.
#Some graphical exploration

ggplot(nestmale, aes(x=RelLength, y=RelEggs))+
  geom_point(aes(fill=Nest_size), shape=21)+
  facet_wrap(~Timeseason)
str(nestmale)

#Model 1 of path analysis => EGGS---------

modp1<-lm(RelEggs~RelLength+RelCond+Nest_size+RelCann+RelGro, data=nestmale[nestmale$Timeseason=="Peak season",])
summary(modp1)
plot(effect("RelLength",modp1))

modp2<-lm(RelEggs~RelLength+RelCond+Nest_size+RelCann+RelGro, data=nestmale[nestmale$Timeseason=="Late season",])
summary(modp2)

#Model 2 of path analysis => nest size---------
str(nestmale)
nestmale$nsizebin<-NA
nestmale[nestmale$Nest_size=="Small nest",]$nsizebin<-0
nestmale[nestmale$Nest_size=="Large nest",]$nsizebin<-1
hist(nestmale$nsizebin)
nestmale[50,c("Nest_size","nsizebin")]
nestmale[300,c("Nest_size","nsizebin")]

modp1.2<-glm(family=binomial(link = "logit"),relevel(Nest_size, ref="Small nest")~RelLength+RelCond, data=nestmale[nestmale$Timeseason=="Peak season",])
summary(modp1.2)
plot(effect("RelLength",modp1.2 ))
plot(effect("RelCond",modp1.2 ))
plot(nestmale$nsizebin~nestmale$RelLength)
plot(modp1.2$fitted.values~nestmale[nestmale$Timeseason=="Peak season",]$RelLength)
plot(modp1.2$fitted.values/(1-modp1.2$fitted.values)~nestmale[nestmale$Timeseason=="Peak season",]$RelLength)


modp1.2.2<-glm(family="binomial",nsizebin~RelLength+RelCond, data=nestmale[nestmale$Timeseason=="Peak season",])
summary(modp1.2.2)
plot(effect("RelLength",modp1.2 ))
plot(effect("RelCond",modp1.2 ))


modp2.2<-glm(family="binomial",relevel(Nest_size, ref="Small nest")~RelLength+RelCond, data=nestmale[nestmale$Timeseason=="Late season",])
summary(modp2.2)

#Model 3 of path analysis => cannibalism--------
modp1.3<-lm(RelCann~RelLength+RelCond+Nest_size, data=nestmale[nestmale$Timeseason=="Peak season",])
summary(modp1.3)

betagro1.3 <- betareg(I(PropCann+0.01)~area+RelLength+RelCond+Nest_size, data = nestmale[nestmale$Timeseason=="Peak season",])
Anova(betagro1.3, type=2)
#++
modp2.3<-lm(RelCann~RelLength+RelCond+Nest_size, data=nestmale[nestmale$Timeseason=="Late season",])
summary(modp1.3)

betagro2.3 <- betareg(I(PropCann+0.01)~area+RelLength+RelCond+Nest_size, data = nestmale[nestmale$Timeseason=="Late season",])
Anova(betagro2.3, type=2)

#Model 4 of path analysis => overgrowth--------
modp1.4<-lm(RelGro~RelLength+RelCond+Nest_size, data=nestmale[nestmale$Timeseason=="Peak season",])
summary(modp1.4)
ggplot(nestmale,aes(x=Nest_size,y=RelGro))+geom_boxplot()

betagro1.4 <- betareg(I(PropGrown+0.01)~area+RelLength+RelCond+Nest_size, data = nestmale[nestmale$Timeseason=="Peak season",])
Anova(betagro1.4, type=2)
plot(effect("Nest_size",betagro1.4 ), multiline=T)
#++
modp2.4<-lm(RelGro~RelLength+RelCond+relevel(Nest_size, ref="Small nest"), data=nestmale[nestmale$Timeseason=="Late season",])
summary(modp2.4)

betagro2.4 <- betareg(I(PropGrown+0.01)~area+RelLength+RelCond+Nest_size, data = nestmale[nestmale$Timeseason=="Late season",])
Anova(betagro2.4, type=2)

#Suppl figure on male size/reproductive success
ggplot(nestmale,aes(x=RelLength,y=RelEggs))+geom_point(aes(color=Nest_size))+
  geom_line(aes(group=Nest_size, color=Nest_size), stat="smooth",method = "lm",  size = 1.5, alpha = 0.6)+
  #geom_line( stat="smooth",color="grey50",method = "lm",  size = 1.5, alpha = 0.6, linetype="dashed")+
  #facet_grid(cols=vars(Timeseason))+
  labs(x="Relative male length",y="Relative egg number")

modrelsize<-lm(RelEggs~RelLength+Nest_size+Timeseason
               +RelLength:Nest_size
               , data=nestmale)
Anova(modrelsize,type=3)
plot(effect("RelLength:Nest_size",modrelsize ), multiline=T)

m.lst.size<- lstrends(modrelsize, ~Nest_size, var="RelLength",adjust = "sidak") #get slopes for each sex
summary(m.lst.size, infer=c(TRUE,TRUE),null=0) #compare each slope to zero 







################################
#  SIDE QUESTS
###############################


str(nestmale)
ggplot(nestmale, aes(x=Eyed_analysis, y=RelCond))+geom_boxplot()+facet_wrap(~area)
ggplot(nestmale, aes(x=Eyed_field, y=RelCond))+geom_boxplot()+facet_wrap(~area)
ggplot(nestmale, aes(x=Eyed_analysis, y=con_all))+geom_boxplot()+facet_wrap(~area)

#############################
############################
##############################
##############################

hist(nestmale$RelCann)
str(nestmale)
#I do not have the power to test for these correlation separately for each pop:nest size:time point combination
#however, it makes sense to get rid of the latitude effect that can hide within pop patterns.
#-------------------------------
#With standardized variables

malecorstd<-nestmale[,c("area","Timeseason","Nest_size","RelLength","RelCann","RelGro","con_pop","RelEggs")]

#overall corrplot
colnames(malecorstd)<-c("area" ,"Timeseason", "Nest_size","Male length","Cannibalism index","Overgrowth index","Male condition", "Egg number")
corrplot(cor(malecorstd[,4:8], use="complete.obs"), 
         type="upper",col = COL2('RdBu', 10),addCoef.col = 'black',tl.col="grey40", diag=F, outline=T)
#-> male size and egg number positively corr
#-> Cannibalism index and egg number negatively corr

#by nest size and season
corrplot(cor(malecorstd[malecorstd$Nest_size=="Large nest"
                     & malecorstd$Timeseason=="Mid-season",4:8], use="complete.obs"), 
         type="upper",col = COL2('RdBu', 10),addCoef.col = 'black',tl.col="grey40", diag=F, outline=T)


corrplot(cor(malecorstd[malecorstd$Nest_size=="Small nest"
                     & malecorstd$Timeseason=="Mid-season",4:8], use="complete.obs"), 
         type="upper",col = COL2('RdBu', 10),addCoef.col = 'black',tl.col="grey40", diag=F, outline=T)

corrplot(cor(malecorstd[malecorstd$Nest_size=="Large nest"
                     & malecorstd$Timeseason=="Late season",4:8], use="complete.obs"),
         type="upper",col = COL2('RdBu', 10),addCoef.col = 'black',tl.col="grey40", diag=F, outline=T)


corrplot(cor(malecorstd[malecorstd$Nest_size=="Small nest"
                     & malecorstd$Timeseason=="Late season",4:8], use="complete.obs"), 
         type="upper",col = COL2('RdBu', 10),addCoef.col = 'black',tl.col="grey40", diag=F, outline=T)

#Relationship between male size and relative egg number
ggplot(nestmale,aes(x=RelLength,y=RelEggs))+geom_point(aes(color=Nest_size))+
  geom_line(aes(group=Nest_size, color=Nest_size), stat="smooth",method = "lm",  size = 1.5, alpha = 0.6)+
  geom_line( stat="smooth",color="grey50",method = "lm",  size = 1.5, alpha = 0.6, linetype="dashed")+
  facet_grid(cols=vars(Timeseason))+
  labs(x="Relative male length",y="Relative egg number")


ggplot(nestmale,aes(x=con_pop,y=RelCann, color=Nest_size, group=Nest_size))+geom_point()+
  geom_line(stat="smooth",method = "lm",  size = 1.5, alpha = 0.6)+
  facet_grid(cols=vars(Timeseason))+  #scale_y_sqrt()+
  labs(x="Relative male length",y="Relative cannibalism")

ggplot(nestmale,aes(x=con_pop,y=RelGro, color=Nest_size, group=Nest_size))+geom_point()+
  geom_line(stat="smooth",method = "lm",  size = 1.5, alpha = 0.6)+
  facet_grid(cols=vars(Timeseason))+  #scale_y_sqrt()+
  labs(x="Relative male length",y="Relative growth")


### ATTEMPTING PATH ANALYSIS


library(lavaan)
library(semPlot)


#Latitude effect on male size and condition, effect of latitude, male size and condition on cannibalism and growth
#problem: need a betareg for cannibalism and growth...

nestmale$slat<-scale(nestmale$latitude)
nestmale$slength<-scale(nestmale$length)
nestmale$scond<-scale(nestmale$con_all)
nestmale$sCann<-scale(nestmale$PropCann)
nestmale$sGro<-scale(nestmale$PropGrown)


path2<-'
sCann ~ slat+slength+scond
sGro ~ slat+slength+scond
slength~slat
scond~slat
 '
path2.mod<-sem(path2, data=nestmale)

summary(path2.mod)

semPaths(path2.mod, what='est',whatLabels = 'est', rotation=1,
         edge.label.cex=1, curvePivot = T, style="mx",layout='spring',residuals=T,
         nCharNodes=0,sizeMan=10,sizeMan2=8,shapeMan='circle',#nodeLabels=c("Eggs","Length","Lat."),
         fade=FALSE)

#second question: reproductive success!
#Instead of using relative data, I can make models with area or locality as a random effect!
#This allows each area/locality to have it's own intercept.
#I can use the betareg models to test for significance when the Cannibalism and overgrowth index are the response variable

#a. Repro success~ Nest status + Male quality + 1|location
ggplot(nestmale,aes(y=Eggs, x=area, group=location))+geom_boxplot() #location seems a good level for a random effect

modRS<-lmer(Eggs~PropGrown+PropCann+length+)

str(nestmale)

#now we will use relative data per population, and run one model for each timepoint and nest size

nestmale$sRelEggs<-scale(nestmale$RelEggs)
nestmale$sRelLength<-scale(nestmale$RelLength)
nestmale$sRelCann<-scale(nestmale$RelCann)
nestmale$sRelGro<-scale(nestmale$RelGro)
nestmale$scon_pop<-scale(nestmale$con_pop)

path3<-'
sRelEggs~sRelLength+scon_pop+sRelCann+sRelGro
sRelCann~scon_pop+sRelLength
sRelGro~scon_pop+sRelLength
 '
path3.mod<-sem(path3, data=nestmale[nestmale$Timeseason=="Mid-season" & nestmale$Nest_size=="Large nest",])

summary(path3.mod)

semPaths(path3.mod, what='est',whatLabels = 'est', rotation=1,
         edge.label.cex=1, curvePivot = T, style="mx",layout='spring',residuals=T,
         nCharNodes=0,sizeMan=10,sizeMan2=8,shapeMan='circle',#nodeLabels=c("Eggs","Length","Lat."),
         fade=FALSE)
#¤%¤%¤
path4.mod<-sem(path3, data=nestmale[nestmale$Timeseason=="Late season" & nestmale$Nest_size=="Large nest",])

summary(path4.mod)

semPaths(path4.mod, what='est',whatLabels = 'est', rotation=1,
         edge.label.cex=1, curvePivot = T, style="mx",layout='spring',residuals=T,
         nCharNodes=0,sizeMan=10,sizeMan2=8,shapeMan='circle',#nodeLabels=c("Eggs","Length","Lat."),
         fade=FALSE)

#¤%¤%¤5
path5.mod<-sem(path3, data=nestmale[nestmale$Timeseason=="Mid-season" & nestmale$Nest_size=="Small nest",])

summary(path5.mod)

semPaths(path5.mod, what='est',whatLabels = 'est', rotation=1,
         edge.label.cex=1, curvePivot = T, style="mx",layout='spring',residuals=T,
         nCharNodes=0,sizeMan=10,sizeMan2=8,shapeMan='circle',#nodeLabels=c("Eggs","Length","Lat."),
         fade=FALSE)

#¤%¤%¤5
path6.mod<-sem(path3, data=nestmale[nestmale$Timeseason=="Late season" & nestmale$Nest_size=="Small nest",])

summary(path6.mod)

semPaths(path6.mod, what='est',whatLabels = 'est', rotation=1,
         edge.label.cex=1, curvePivot = T, style="mx",layout='spring',residuals=T,
         nCharNodes=0,sizeMan=10,sizeMan2=8,shapeMan='circle',#nodeLabels=c("Eggs","Length","Lat."),
         fade=FALSE)


#`·.¸¸.·´¯`·.¸¸.·´¯`·.¸ ><((((º> ¸.·´¯`·.¸¸.·´¯`·.¸¸.·´¯`·.¸ ><((((º>. ¸.·´¯`·.¸¸.·´¯`·.¸¸.·´¯`·.¸ ><((((º>.¸.·´¯`·.¸¸.·´¯`·.¸¸.·´¯`·.¸¸.·
#,¸.·´¯`·.¸¸.·´¯`·.¸ ¸.·´¯`· ><((((º>.¸¸.·´¯`·.¸¸.·´¯`·.¸ ><((((º>. ¸.·´¯`·.¸¸.·´¯`·.¸¸.·´¯`·.¸ ><((((º>.¸.·´¯`·.¸¸.·´¯`·.¸¸.·´¯`·.¸¸.·

########################################################################
# Appendix. test panel
#######################################################################
# In order to test the nest picture analysis protocol, 20 nest images from
#different locations and time periods were selected and analysed twice
#The aim of this section is to test for repetability
#To do this, I propose to run a linear model on the variables of interest
#the model will have Nest Size as a fixed effect (otherwise this will create
#a strong effect on total egg number) and replicate as a random effect
#We then compare the variance attributed to replicate to the residual variance

#-------------------
#import the data 
testnest<-read.delim(header=T,"testpanel.txt")
testnest$Area<-as.factor(testnest$Area)
testnest$NestSize<-as.factor(testnest$NestSize)
testnest$replicate<-as.factor(testnest$replicate)
str(testnest)

#--------------------
#we need to make some calculations to obtain the variables we want to analyse
# The total brood surface in square unit (comparable across nests)
testnest$TotBrood<-testnest$TotBrood_surface/testnest$X1sq_surface
# The cannibalised area and overgrown area in square unit
testnest$Cannib<-testnest$Cann_surface/testnest$X1sq_surface
testnest$Overgrown<-testnest$Grown_surface/testnest$X1sq_surface
# The area still covered with eggs within the brood area
testnest$Brood<-testnest$TotBrood-(testnest$Cannib+testnest$Overgrown)
# Egg counting area in square unit (5 areas counted)
testnest$Count_Area<-(5*testnest$Egg_surface)/testnest$X1sq_surface
# Egg number in total brood area
testnest$Totegg<-(testnest$TotBrood/testnest$Count_Area)*testnest$Egg_Nr
# Egg number left after cannibalism and overgrowth
testnest$Eggs<-(testnest$Brood/testnest$Count_Area)*testnest$Egg_Nr
# Cannibalised proportion of brood
testnest$Prop_cann<-testnest$Cann_surface/testnest$TotBrood_surface
# Overgrown proportion of brood
testnest$Prop_grown<-testnest$Grown_surface/testnest$TotBrood_surface

#Proportion of nest covered: requires nr of squares avaiable for large and small nests

#----------------------
#Test REPETABILITY
#graphical exploration
#Create a data frame that allows to plot the two replicates against each other.
plottestpanel<-rbind(
  data.frame(replicate1=testnest[testnest$replicate==1,]$TotBrood,replicate2=testnest[testnest$replicate==2,]$TotBrood, Var="TotBrood"),
  data.frame(replicate1=testnest[testnest$replicate==1,]$Brood,replicate2=testnest[testnest$replicate==2,]$Brood, Var="Brood"),
  data.frame(replicate1=testnest[testnest$replicate==1,]$Eggs,replicate2=testnest[testnest$replicate==2,]$Eggs, Var="Eggs"),
  data.frame(replicate1=testnest[testnest$replicate==1,]$Prop_cann,replicate2=testnest[testnest$replicate==2,]$Prop_cann, Var="Prop_cann"),
  data.frame(replicate1=testnest[testnest$replicate==1,]$Prop_grown,replicate2=testnest[testnest$replicate==2,]$Prop_grown, Var="Prop_grown")
  )
str(plottestpanel)

ggplot(data=plottestpanel,aes(x=replicate1,y=replicate2))+geom_point()+geom_abline(slope=1,intercept=0)+
  facet_wrap(~Var, scales = "free", nrow=1)

#Total brood surface (in unit of square)
#Including the fixed effect of nest size will be necessary so as to not inflate the residual variance
#simple model with fixed effects and ANOVA
mtotbrood<-lm(TotBrood~NestSize+replicate, data=testnest)
Anova(mtotbrood)
# Using the formula for approximative repeatbility based on F value from Lessells and Boag 1987
Fv<-0.0128
(1/Fv-1)/(1/Fv-1+(39/38))

#Surface with eggs (brood left)
mbrood<-lm(Brood~NestSize+replicate, data=testnest)
Anova(mbrood)
Fv<-0.0035
(1/Fv-1)/(1/Fv-1+(39/38))

#Total number of eggs
meggs<-lm(Eggs~NestSize+replicate, data=testnest)
Anova(meggs)
Fv<-0.003
(1/Fv-1)/(1/Fv-1+(39/38))

#Proportion cannibalised
#Ideally it would be a zero inflated betareg, but I cannot calculate repeatability easily from these
mpcan<-gamlss(Prop_cann ~ NestSize+replicate,  family = BEZI, data = testnest, trace = F)
#Instead I suggest to not use the proportion of the nest that is cannibalised, but the absolute surface in squares
#In further analysis, I will use the proportion, but to test repeatability it isn't necessary
modcan<-lm(Cannib ~NestSize+replicate, data=testnest)
Anova(modcan)
Fv<-0.1927
(1/Fv-1)/(1/Fv-1+(39/38))

#Surface  overgrown
modgrown<-lm(Overgrown ~NestSize+replicate, data=testnest)
Anova(modgrown)
Fv<-0.1278
(1/Fv-1)/(1/Fv-1+(39/38))



# JUNK DNA

#Singular fit with lme4 suggests the variance due to replicate is too small to be estimated
modtotbrood<-lmer(TotBrood~NestSize+(1|replicate), data=testnest)
#Turn to bayesian
baytotbrood<-brm(TotBrood~NestSize+(1|replicate), data=testnest,
                  init=0,iter=10000, cores=8,control = list(adapt_delta = 0.95)) #technical issue with Stan
summary(baytotbrood)
plot(baytotbrood)



#Surface with eggs (brood left)
#Singular fit with lme4 suggests the variance due to replicate is too small to be estimated
modbrood<-lmer(Brood~NestSize+(1|replicate), data=testnest)
#Turn to bayesian
baybrood<-brm(Brood~NestSize+(1|replicate), data=testnest,
                 init=0,iter=10000, cores=8,control = list(adapt_delta = 0.95)) #technical issue with Stan
summary(baybrood)
plot(baybrood)


#Total number of eggs
#Singular fit with lme4 suggests the variance due to replicate is too small to be estimated
modeggs<-lmer(Eggs~NestSize+(1|replicate), data=testnest)
#Turn to bayesian
bayeggs<-brm(Eggs~NestSize+(1|replicate), data=testnest,
              init=0,iter=10000, cores=8,control = list(adapt_delta = 0.95)) #technical issue with Stan
summary(bayeggs)
plot(bayeggs)


#Proportion of brood cannibalised
#The response being a proportion not resulting from bernoulli trials, I need a beta regression.
#This is not possible as far as I know with lme4, and betareg won't do mixed models
#Turn to bayesian
baypcan<-brm(Prop_cann~NestSize+(1|replicate), data=testnest,family=zero_inflated_beta(link = "logit", link_phi = "log", link_zi = "logit"),
             init=0,iter=10000, cores=8,control = list(adapt_delta = 0.95)) #technical issue with Stan
summary(baypcan)
plot(baypcan)

#NOTE: if doing more zero inflated beta reg, go to this page: https://www.andrewheiss.com/blog/2021/11/08/beta-regression-guide/

#Proportion of brood overgrown
baypbrood<-brm(Prop_grown~NestSize+(1|replicate), data=testnest,family=zero_inflated_beta(link = "logit", link_phi = "log", link_zi = "logit"),
             iter=10000, cores=8,control = list(adapt_delta = 0.95)) #technical issue with Stan
summary(baypbrood)
plot(baypbrood)




############### junk dna
#Part of the model takes care of the zeroes in a binomial way (the nest is overgrown or not)
#part of the model asks, how much is it overgrown when it is.
#zoib package
nestdatclean<-nestdat[ !is.na(nestdat$CannibCm) ,]
str(nestdatclean)
cani.bay <- zoib(PropCann ~ Timeseason +Nest_size+ area| Timeseason +Nest_size+ area | Timeseason +Nest_size+ area, data = nestdatclean, joint = T,
                 random = 0, EUID=1:nrow(nestdatclean) , zero.inflation = T, 
                 one.inflation = F, n.iter = 1050, n.thin = 5, n.burn = 50)
#About the zoib formula: on the right hand side, there can be 5 parts. separated by |
# First estimated the mean of the beta dist for specificed effects
# Precision parameter of the beta distribution for specified effects
# Zero inflation for specificed parameters
# same for 1 inflation
# random effects
#So the formula specified abofe had the mean of the beta distribution and the zero infl parameters estimated for the 3 effects
#Putting "1" in the middle section would give only an intercept to the precision parameter

sample1 <- cani.bay$coeff
# check  convergence of the MCMC chains
traceplot(sample1)
autocorr.plot(sample2)
gelman.diag(sample1)
summary(sample1)$stat[, 1]
zoib::paraplot(cbind(summary(sample1)$stat[, 1],summary(sample1)$quant[, c(1, 5)]))




#What about glmmTMB? still with offset data but includes the random effect
Tmbeta1<-glmmTMB(PropGoff~scale(DDF)+scale(poplength)+relevel(Nest_size, ref="Small nest")+scale(withlength)+scale(withcon)+(1|location), 
                 data = nestanalysis[nestanalysis$Timeseason=="Peak season",],
                 family = beta_family(link = "logit") )
summary(Tmbeta1)
Anova(Tmbeta1)
