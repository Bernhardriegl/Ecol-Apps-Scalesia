#VERSION FOR SUBMISSION
#Scalesia_file; Bernhard Riegl, started on 5/27/2019, first submitted 3/21/21; fixed, updated and cleaned by 10/19/21
#resubmitted 11/4/2021, modified after review 1/10/2022 and  check 3/30/2022
# next change for resubmission 8/9/2022: added DBH analyses, nLMEr and filepath simplification

#Map of the Galapagos for MAIN TEXT FIGURE 1
library(maps)
library(marmap)
library(mapproj)
# This makes the world map 
windows(5,5)
map("world", regions=".",exact=F, boundary=T,interior=T,proj="ortho",orientation=c(10,-80,0),fill=T,wrap=c(0,360))
#Galapagos map
Gal<-marmap::getNOAA.bathy(lon1=-92.5,lon2=-88.5,lat1=-2,lat2=2,resolution=1)
greys<-c("grey80", "grey85","grey90", "grey95")
windows(10,10)
plot(Gal, image = TRUE, land = TRUE, lwd = 0.005,col="grey",axes=FALSE,xlab="",ylab="",
     bpal = list(c(0, max(Gal), "grey60"),c(min(Gal),0,greys)))
plot(Gal, deep = 0, shallow=0, step=0,lwd=2,col="grey50",add=TRUE)

#now for the Scalesia data and their evaluation
#1=Treatment, 2=Control
library(tidyverse)
library(ggpubr)
library(nlme)
library(lme4)

#PART 1: GROWTH MODEL FOR SCALESIA PEDUNCULATA
setwd("G:/University/PAPERS/E-Pacific/Scalesia/Data/")
DBH<-read.table("Scalesia_DBF-TH.txt",header=T) # DBH and Height at random locations in forest prior to beginning of experiment in 2014
DBH_Treat<-read.table("Scalesia-DBH-Height_Treatment.txt",header=T)# BH and Height of haphazardly chosen thees in Treatment plots in 2020 (Rubus removed)
DBH_Control<-read.table("Scalesia-DBH-Height_Control.txt",header=T)# BH and Height of haphazardly chosen thees in Control plots in 2020 (Rubus remains)

#calculated for the first DBH dataset, taken from random spots within forest IN 2014
# calculate the parameters for the non-linear model, which is y=a-b*exp(-cx)
# Note: this is standard expression of the asymptotic regression model a-(a-b)exp(-cX) 
#under constraint that b=0, it turns into the negative exponential equation Y=a*(1-exp(-cX))
# the negative exponential equation passes through the origin!
#a=approximate asymptotic max height=1000
a=1000
b=920
y=200 #where curve is rising most steeply, height is greatest...next line
x=2 #where x=2
c=(log((a-y)/b))/x
#===>APPENDIX 1, TABLE S-1
model1<-nls(Height_cm ~ a-b*exp(-c*DBH_cm),data=DBH, start=list(a=1000,b=920,c=0.12))#3 parameter asymptotic exponential
summary(model1)
#===>APPENDIX 1, TABLE S-2
model2<-nls(Height_cm ~ a*(1-exp(-c*DBH_cm)),data=DBH,start=list(a=1000,c=0.12))#2 parameter asymptotic exponential
summary(model2)
AIC(model1)
anova(model1,model2)
#===>APPENDIX 1, TABLE S-3
#asymptotic exponential uses more parameters, while the 2-parameter asymptotic correctly puts the curve through the origin
av<-seq(0,30,0.29)
bv1<-predict(model1,list(DBH_cm=av))
bv2<-predict(model2,list(DBH_cm=av))
#now do the same for the treatment area and the control area DATA IN 2020

model3<-nls(Height_cm ~ a-b*exp(-c*DBH_cm),data=DBH_Treat, start=list(a=1000,b=920,c=0.12))#3 parameter asymptotic exponential
model4<-nls(Height_cm ~ a*(1-exp(-c*DBH_cm)),data=DBH_Treat,start=list(a=1000,c=0.12))#2 parameter asymptotic exponential
anova(model4,model3)#specify the reduced model first to avoid negative values

model5<-nls(Height_cm ~ a-b*exp(-c*DBH_cm),data=DBH_Control, start=list(a=1000,b=920,c=0.12))#3 parameter asymptotic exponential
model6<-nls(Height_cm ~ a*(1-exp(-c*DBH_cm)),data=DBH_Control,start=list(a=1000,c=0.12))#2 parameter asymptotic exponential
anova(model6,model5)
#now calculate a t-test between Treatment and Control in height and DBH, in main text
# NOTE: t.test in R is Welch's, for unequal variances!
t.test(DBH_Treat$DBH_cm,DBH_Control$DBH_cm)
t.test(DBH_Treat$Height_cm,DBH_Control$Height_cm)

#=====MAIN TEXT, FIGURE 2=====
av<-seq(0,30,0.29)
bv3<-predict(model3,list(DBH_cm=av))
bv4<-predict(model4,list(DBH_cm=av))
bv5<-predict(model5,list(DBH_cm=av))
bv6<-predict(model6,list(DBH_cm=av))
windows(10,5)
plot(DBH$DBH_cm,DBH$Height_cm,pch=21,col="blue",bg="blue",xlab="DBH in cm",ylab="Height in cm",xlim=c(0,30),ylim=c(0,1000))
#lines(av,bv1,col="blue",lwd=2,lty=2)
lines(av,bv2,col="blue",lwd=2)
text(27.5,1000, "pre-experiment - 2014", col="blue")
par(new=T)
plot(DBH_Treat$DBH_cm,DBH_Treat$Height_cm,pch=21,col=rgb(0,0.6,0.5),bg=rgb(0,0.6,0.5),xlab="DBH in cm",ylab="Height in cm",xlim=c(0,30),ylim=c(0,1000))
lines(av,bv3,col=rgb(0,0.6,0.5),lwd=2)
lines(av,bv4,col=rgb(0,0.6,0.5))
text(27,850,expression(paste(italic('R. niveus')," removed - 2020")),col=rgb(0,0.6,0.5))
par(new=T)
plot(DBH_Control$DBH_cm,DBH_Control$Height_cm,pch=21,col=rgb(0.9,0.6,0),bg=rgb(0.9,0.6,0),xlab="DBH in cm",ylab="Height in cm",xlim=c(0,30),ylim=c(0,1000))
lines(av,bv5,col=rgb(0.9,0.6,0),lwd=2)
lines(av,bv6,col=rgb(0.9,0.6,0))
text(27,700,expression(paste(italic('R. niveus')," present - 2020")),col=rgb(0.9,0.6,0))

#Step 1: pool data from all plots and see if there is a difference pre/post experiment and Treatment/Control  w.r.t. height and DBH
DB.H=cbind(rep(1,dim(DBH)[1]),rep(1,dim(DBH)[1]),DBH);colnames(DB.H)<-c("Phase","Group","DBH","Height") #before experiment
DB.H_Treat=cbind(rep(2,dim(DBH_Treat)[1]),rep(2,dim(DBH_Treat)[1]),DBH_Treat[,c(3,4)]); colnames(DB.H_Treat)<-c("Phase","Group","DBH","Height") #Rubus removed
DB.H_Control=cbind(rep(2,dim(DBH_Control)[1]),rep(3,dim(DBH_Control)[1]),DBH_Control[,c(3,4)]); colnames(DB.H_Control)<-c("Phase","Group","DBH","Height") #Rubus present

HTS<-rbind(DB.H,DB.H_Treat,DB.H_Control)#all the heights and DBHs with an identifier
HTS1<-lapply(HTS,as.numeric)#make all colums numeric, last one was integer
#exploration plot just to see how the groups will fall
plot(HTS$Height~HTS$DBH,bty="n",type="n")
text(HTS$DBH,HTS$Height,HTS$Group)
#now, since the model is clearly nonlinear (but we knew that already...it saturates), use nlme to evaluate interaction of
#experimental phase (before/after begin of removal) and treatments (and)Control versus Treatment), nonlinear regression model, least squares
mymodel<-nlsList(Height~a*(1-exp(-c*DBH))|Group,data=HTS,start=c(a=1000,c=0.12) )
summary(mymodel)

HTS2<-groupedData(Height~DBH | Phase,data=HTS)
#HTS2<-HTS[-which(HTS2$Height==0),]#if you want to remove the 0-heights
#set up control list for LME
control.list <- lmeControl(maxIter = 500, msMaxIter = 500, msMaxEval=500,tolerance = 0.1, msTol = 0.1, sing.tol=1e-20)
#Description: factor a (asymp size) depends on Group, c depends on nothing. a and c very randomly across Phases and Groups
mm2a<-nlme(Height~a*(1-exp(-c*DBH)),fixed=list(a~Group,c~1),random=a+c~1|Phase/Group,data=HTS2,start=c(a=1000,c=0.12,Group=1), method="ML",control=control.list)
#Description: 
#the random variable should be Phase only, Group has a fixed effect on a, c is a constant
mm2a<-nlme(Height~a*(1-exp(-c*DBH))-1,fixed=list(a~Group,c~1),random=a+c~1|Phase,data=HTS2,start=c(a=1000,c=0.12,Group=1), method="ML",control=control.list)
#here, Group has a fixed effect on a, but a and c show random variability from Phase to Phase and Group to Group
mm2b<-nlme(Height~a*(1-exp(-c*DBH))-1,fixed=list(a~Group,c~1),random=a+c~1|Phase/Group,data=HTS2,start=c(a=1000,c=0.12,Group=1), method="ML",control=control.list)
#here, no fixed effect on a or c, but a and c show random variability from Phase to Phase and Group to Group
mm2c<-nlme(Height~a*(1-exp(-c*DBH))-1,fixed=a+c~1,random=a~1|Phase,data=HTS2,start=c(a=1000,c=0.12), method="ML",control=control.list)
AIC(mm2a,mm2b,mm2c)
anova(mm2a,mm2b,mm2c)

summary(mm2a)$tTable
VarCorr(mm2a)

#==>APPENDIX 1: FIGURE S1====================================
#below plot shows that the three groups pr-, post-control, post-treatment differ in slope and asymptote
windows(10,5)
plot(augPred(mm2b))
#===========================================================
summary(mm2d)
#in the summary the Fixed effect values are the means of the parameter values (Crawley p.725), so get each:
coef(mm2a)
#this shows: trees were larger prior to experiment, larger when Rubus removes, smallest with Rubus present
#but it is also clear that basically no recruitment was observed post-2014 (while 2014 had many)
#Solution cut the dataset at 5cm DBH

#Step 2, not in paper: since there are no recruits after 2015, see if allometry of big trees changes
HTS_c<-HTS[-which(HTS$DBH<8),]
plot(HTS_c$Height~HTS_c$DBH,bty="n",type="n")
text(HTS_c$DBH,HTS_c$Height,HTS$Group)
HTS_c2<-groupedData(Height~DBH | Phase,data=HTS)
mymodel_c2<-nlme(Height~a*(1-exp(-c*DBH)),fixed=a+c~1,random=a+c~1|Phase/Group,data=HTS_c2,start=c(a=1000,c=0.12) )
#below plot shows that the three groups pr-, post-control, post-treatment differ in slope and asymptote
windows(10,5)
plot(augPred(mymodel_c2))
summary(mymodel_c2)

#=>RESULT: above analyses show that across the experiment, tree allometry has changed, and that the effect is greatest if Rubus is present


#======> PART 2: GROWTH RATES (DBH-RATIOS) AND THE INFLUENCE OF RUBUS COMPETITION
SCAL.T<-read.table("Scalesia_DBH_Treatment2.txt", header=T)
SCAL.RC<-read.table("Scalesia_DBH_RubusControl2.txt", header=T)
# remove Sep-14 and Aug-15
SCAL.T<-SCAL.T[,-c(4,6)]#cut out September and August measurements, so that only annual February is left
SCAL.T[is.na(SCAL.T)]<-0 #the dead trees are coded NA, make 0
SCAL.RC<-SCAL.RC[,-c(4,6)]
SCAL.RC[is.na(SCAL.RC)]<-0
# add a counter 1=TREATMENT, 2=RUBUS CONTROL
group=c(rep(1,dim(SCAL.T)[1]),rep(2,dim(SCAL.RC)[1]))
SCAL<-rbind(SCAL.T,SCAL.RC)
SCAL<-cbind(group,SCAL)
# calculate the growth rates
SC1<-SCAL[,c(4:10)]
SC2<-SCAL[,c(5:11)]
GRW<-SC2-SC1
GRW[GRW<0]<-0
GRW<-cbind(SCAL$group,SCAL$Plot,SCAL$tree,GRW)
colnames(GRW)<-c("group","Plot","tree","Feb15","Feb16","Feb17","Feb18","Mar19","Mar20","Mar21")

Grw<-GRW%>%
  gather('Feb15','Feb16','Feb17','Feb18','Mar19','Mar20','Mar21',key="year",value="Growth.rate")
Scal<-SCAL%>%
  gather('Feb14','Feb15','Feb16','Feb17','Feb18','Mar19','Mar20','Mar21',key="year",value="DBH")
#=======THIS IS GROWTH AS DBHt2-DBHt1===================
ScalGrw<-left_join(Scal,Grw)
#=======================================================
Grw.pos<-filter(ScalGrw,Growth.rate>0)
#for tree-growth model we look at the growth differently. We say "a tree at DBH x, will grow x cm this year"
#that way, we don't have a growth rate for the last measured year, which we prefer over not having it for the first measured year
SC1<-SCAL[,c(4:10)]
SC2<-SCAL[,c(5:11)]
GRW<-SC2/SC1 #the ratio is year2/year1; this is now relative growth
GRW[GRW<0]<-0
GRW<-cbind(SCAL$group,SCAL$Plot,SCAL$tree,GRW)
colnames(GRW)<-c("group","Plot","tree","Feb14","Feb15","Feb16","Feb17","Feb18","Mar19","Mar20")
Grw2<-GRW%>%
  gather('Feb14','Feb15','Feb16','Feb17','Feb18','Mar19','Mar20',key="year",value="Growth.rate")
Scal2<-SCAL%>%
  gather('Feb14','Feb15','Feb16','Feb17','Feb18','Mar19','Mar20','Mar21',key="year",value="DBH")
#========THIS IS THE DBH-ratio=======================================
ScalGrw2<-left_join(Scal2,Grw2)

#=========data for FIGURE 4 MAIN TEXT====================================
ScalGrwT<-filter(ScalGrw2,group==1 & year!="Mar21")#TREATMENT
ScalGrwC<-filter(ScalGrw2,group==2 & year!="Mar21")#WITH RUBUS CONTROL
#there are many NAs and zeros in there, they make no sense as growth rate
ScalGrwT[is.na(ScalGrwT)]<-0
ScalGrwC[is.na(ScalGrwC)]<-0
ScalGrwT[ScalGrwT==Inf]<-0
ScalGrwC[ScalGrwC==Inf]<-0
SG_T<-ScalGrwT[-which(ScalGrwT$Growth.rate==0),]
SG_C<-ScalGrwC[-which(ScalGrwC$Growth.rate==0),]

#Relationship of DHB-ratio with diameter as second-order polynomial in main text
# APPENDIX 1, TABLE S-6
Gr_mod1<-lm(Growth.rate~poly(DBH,2),data=SG_C)
summary.lm(Gr_mod1)
summary.aov(Gr_mod1)
#or fit exponential decay
# Select an approximate $\theta$, since theta must be lower than min(y), and greater than zero
theta.0 <- min(SG_C$Growth.rate) * 0.5  
# Estimate the rest parameters using a linear model. Note: this concerns only the pooled DBH values overall; not year, plot, etc...
model.0 <- lm(log(Growth.rate - theta.0) ~ DBH, data=SG_C)  
alpha.0 <- exp(coef(model.0)[1])
beta.0 <- coef(model.0)[2]
# Starting parameters
start <- list(alpha = alpha.0, beta = beta.0, theta = theta.0)
#if error is additive, then it is constant with x-axis and we should not use log-scale
Gr_mod2 <- nls(Growth.rate ~ alpha * exp(beta * DBH) + theta , data = SG_C, start = start)
summary(Gr_mod2)
#lines(SG_C$DBH,predict(Gr_mod5,list(x=SG_C$DBH)),col='skyblue',lwd=3)
AIC(Gr_mod2,Gr_mod1)

Gr_mod3<-lm(Growth.rate~poly(DBH,2),data=SG_T)
summary.lm(Gr_mod3)
#Gr_mod4<-lm(log(Growth.rate)~DBH,data=SG_T)#if error multiplicative, grows with x-axis, then it is constant on log-scale
theta.0 <- min(SG_C$Growth.rate) * 0.5  
model.0 <- lm(log(Growth.rate - theta.0) ~ DBH, data=SG_C)  
alpha.0 <- exp(coef(model.0)[1])
beta.0 <- coef(model.0)[2]
start1 <- list(alpha = alpha.0, beta = beta.0, theta = theta.0)
start2<-list(alpha = alpha.0, beta = beta.0)
Gr_mod4 <- nls(Growth.rate ~ alpha * exp(beta * DBH) + theta , data = SG_C, start = start1)
Gr_mod5 <- nls(Growth.rate ~ alpha * exp(beta * DBH) , data = SG_C, start = start2)
anova(Gr_mod4,Gr_mod5)#model simplification not justified, keep theta or model is signifcantly worse
AIC(Gr_mod4,Gr_mod3)

#==> APPENDIX 1, FIGURE S2 <====
# function for computing mean, DS, max and min values
#note, the correct equation for the 95%SI would be mean+/-SI*t(0.975,n-1)
min.mean.sd.max <- function(x) {
  r <- c(min(x), mean(x) - sd(x), mean(x), mean(x) + sd(x), max(x))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r}
#excluding the extremes, a zero-ratio makes no sense. 
p1<-ggplot(data=filter(ScalGrwC, Growth.rate>0.5 & Growth.rate<1.5),aes(x=year,y=Growth.rate))+
  stat_summary(fun.data = min.mean.sd.max, geom = "boxplot")+
  geom_jitter(position=position_jitter(width=.2), size=1)+
  ggtitle(c(expression(italic("R. niveus")*" present, Mean, SD, min and max")))+
  facet_wrap(~Plot,nrow=3)+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90))+
  labs(y="DBH-ratios")
p2<-ggplot(data=filter(ScalGrwT, Growth.rate>0.5 & Growth.rate<1.5),aes(x=year,y=Growth.rate))+
  stat_summary(fun.data = min.mean.sd.max, geom = "boxplot")+
  geom_jitter(position=position_jitter(width=.2), size=1)+
  ggtitle(c(expression(italic("R. niveus")*" removed, Mean, SD, min and max")))+
  facet_wrap(~Plot,nrow=3)+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90))+
  labs(y="DBH-ratio")
#now do mean growth rate against DBH
p3<-ggplot(data=filter(ScalGrwC, Growth.rate>0.5 & Growth.rate<1.5),aes(x=DBH,y=Growth.rate))+
  geom_jitter(position=position_jitter(width=.2), size=1)+
  stat_smooth(mapping=aes(x=DBH,y=Growth.rate),method="lm",formula=y~poly(x,2))+
  ggtitle(c(expression(italic("R. niveus")*" present")))+
  facet_wrap(~Plot,nrow=3)+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90))+
  labs(y="DBH-ratio")
p4<-ggplot(data=filter(ScalGrwT, Growth.rate>0.5 & Growth.rate<1.5),aes(x=DBH,y=Growth.rate))+
  geom_jitter(position=position_jitter(width=.2), size=1)+
  stat_smooth(mapping=aes(x=DBH,y=Growth.rate),method="lm",formula=y~poly(x,2))+
  ggtitle(c(expression(italic("R. niveus")*" removed")))+
  facet_wrap(~Plot,nrow=3)+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90))+
  labs(y="DBH-ratio")
windows(25,20)
ggarrange(p1,p2,p3,p4,labels=c("A","B","C","D"),ncol=2,nrow=2)

# ====> MAIN TEXT, Figure 3  <========
windows(11,5)
par(mfrow=c(1,2))
plot(SG_T$DBH,SG_T$Growth.rate,xlab="DBH",ylab="DBH-ratio",ylim=c(0.9,1.55),col="green4")
text(22.8,1.5,expression(paste(italic('R. niveus')," removed")),col="green4",cex=1.5)
x<-0:35
y<-predict(Gr_mod4,list(DBH=x))
lines(x,y,col="green4",lwd=3)
mtext('A',side=3,line=1,at=0.6,cex=2)   

plot(SG_C$DBH,SG_C$Growth.rate,xlab="DBH",ylab="DBH-ratio",ylim=c(0.9,1.55),col="orange3")
text(23.6,1.5,expression(paste(italic('R. niveus')," present")),col="orange3",cex=1.5)
x<-0:35
y<-predict(Gr_mod2,list(DBH=x))
lines(x,y,col="orange3",lwd=3)
mtext('B',side=3,line=1,at=0.6,cex=2)


#==> MAIN TEXT, FIGURE 4<====
# function for computing mean, DS, max and min values
#note, the correct equation for the 95%SI would be mean+/-SI*t(0.975,n-1)
min.mean.sd.max <- function(x) {
  r <- c(min(x), mean(x) - sd(x), mean(x), mean(x) + sd(x), max(x))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r}
#excluding the extremes, a zero-ratio makes no sense. 
#leave the facet_wrap active to see every plot seperately
p1<-ggplot(data=filter(ScalGrwC, Growth.rate>0.5 & Growth.rate<1.5),aes(x=year,y=Growth.rate))+
  stat_summary(fun.data = min.mean.sd.max, geom = "boxplot")+
  geom_jitter(color=rgb(0.9,0.6,0),position=position_jitter(width=.2), size=1)+
  ggtitle(c(expression(italic("R. niveus")*" present")))+
  #facet_wrap(~Plot,nrow=3)+
  theme_bw()+ylim(c(0.7,1.5))+
  theme(axis.text.x=element_text(angle=90))+
  labs(y="DBH-ratio")
p2<-ggplot(data=filter(ScalGrwT, Growth.rate>0.5 & Growth.rate<1.5),aes(x=year,y=Growth.rate))+
  stat_summary(fun.data = min.mean.sd.max, geom = "boxplot")+
  geom_jitter(color=rgb(0,0.6,0.5), position=position_jitter(width=.2), size=1)+
  ggtitle(c(expression(italic("R. niveus")*" removed")))+
  #facet_wrap(~Plot,nrow=3)+
  theme_bw()+ylim(c(0.7,1.5))+
  theme(axis.text.x=element_text(angle=90))+
  labs(y="DBH-ratio")
windows(25,15)
ggarrange(p2,p1,labels=c("A","B"),ncol=2,nrow=1)    

# Q1: Does relationship DBH-ratio/DBH, over all datapoints irrespective of plot and year, change with treatment?
# REMOVE # APPENDIX 1, TABLE S-5<================

#Use all data in Linear Mixed Effects Model with bin as Random effect, Treatment as fixed effect. 
#Q2: does Treatment have an effect on DBH-ratio overall? At this point not interested in "year".
#log(Growth rate) bc. we have established exponential relationship!
m2<-lmer(log(Growth.rate)~DBH*group+(1|group/year/Plot),REML=FALSE,data=Grw.pos)#no REML to allow Anova comparison
m3<-lmer(log(Growth.rate)~DBH+group+year+(1|group/year/Plot),REML=FALSE,data=Grw.pos)#no REML to allow Anova comparison
#try with lme, same thing but different tables

#check best model structure
#random intercept; also: # DBH*Group=>DBH+group+DBH:Group ->two-way interaction
m2_1b<-lme(log(Growth.rate)~1, random=~1|year/Plot, data=Grw.pos, method="ML")
m2_1a<-lme(log(Growth.rate)~DBH+group, random=~1|year/Plot,data=Grw.pos,method="ML")# random intercepts based on group, no DBH-Group interaction
m2_1a<-lme(log(Growth.rate)~DBH*group, random=~1|year/Plot,data=Grw.pos,method="ML")# random intercepts based on group, with DBH-Group interaction

#random Intercept random=~1|bla
m2_1<-lme(log(Growth.rate)~1, random=~1|group/year/Plot, data=Grw.pos, method="ML")
m2_1a<-lme(log(Growth.rate)~DBH+group, random=~1|year/Plot,data=Grw.pos,method="ML")#ANCOVA:DBH and Group have fixed effect; Year and plot have random effect on intercept
m2_1aa<-lme(log(Growth.rate)~DBH+group, random=~1|year,data=Grw.pos,method="ML")#ANCOVA:DBH and Group have fixed effect; Year and plot have random effect on intercept
m2_1b<-lme(log(Growth.rate)~DBH*group, random=~1|year/Plot,data=Grw.pos,method="ML")#ANCOVA:DBH and Group have fixed effect; Year and plot have random effect on intercept
m2_1bb<-lme(log(Growth.rate)~DBH*group, random=~1|year,data=Grw.pos,method="ML")#ANCOVA:DBH and Group have fixed effect; Year and plot have random effect on intercept
AIC(m2_1,m2_1a,m2_1aa,m2_1b,m2_1bb)#m2_1b is better
anova(m2_1,m2_1a,m2_1aa,m2_1b,m2_1bb)

#random intercept and slope model (ANCOVA)
control.list <- lmeControl(maxIter = 500, msMaxIter = 500, msMaxEval=500,tolerance = 0.1, msTol = 0.1, sing.tol=1e-20)
m2_2a<-lme(log(Growth.rate)~DBH+group, random=~1+group|year/Plot, data=Grw.pos, method="ML")#a random intercept AND slope for each group, no DBH-Group interaction
m2_2aa<-lme(log(Growth.rate)~DBH+group, random=~1+group|year, data=Grw.pos, method="ML")#a random intercept AND slope for each group, no DBH-Group interaction
m2_2aaa<-lme(log(Growth.rate)~DBH+group, random=~1+group|Plot, data=Grw.pos, method="ML",control=control.list)#a random intercept AND slope for each group, no DBH-Group interaction
m2_2b<-lme(log(Growth.rate)~DBH*group, random=~1+group|year/Plot, data=Grw.pos, method="ML")#a random intercept AND slope for each group, with DBH-Group interaction
m2_2bb<-lme(log(Growth.rate)~DBH*group, random=~1+group|year, data=Grw.pos, method="ML")#a random intercept AND slope for each group, with DBH-Group interaction
m2_2bbb<-lme(log(Growth.rate)~DBH*group, random=~1+group|Plot, data=Grw.pos, method="ML")#a random intercept AND slope for each group, with DBH-Group interaction

AIC(m2_2a,m2_2aa,m2_2aaa,m2_2b,m2_2bb,m2_2bbb)#m2_2aa has lowest combination AIC df
anova(m2_2a,m2_2aa,m2_2aaa)
anova(m2_2b,m2_2bb,m2_2bbb)


#the chosen model; present with REML
m2_2<-lme(log(Growth.rate)~DBH*group, random=~1+group|year/Plot, data=Grw.pos, method="ML")#a random intercept AND slope for each group
100*vars/sum(vars)
summary(m2_2)
#also do a lmer to get output and calculate where the variance in the random components is (Crowley, p. 704)
mod2_2aa<-lmer(log(Growth.rate)~DBH*group+(1+group|year),REML=T, data=Grw.pos)
#mod2_2aa<-lmer(log(Growth.rate)~DBH+group+(1|year),REML=T, data=Grw.pos)#random intercept only
summary(mod2_2aa)
vars<-c(0.32,0.06,0.47)
100*vars/sum(vars)
summary(m2_2)#APPENDIX 1, TABLE S-7<======================

#multiple variance structure doesn't help
#m2_2b<-lme(log(Growth.rate)~DBH*group, random=~1+group|year/Plot, weights=varIdent(form=~1|year),data=Grw.pos, method="REML")#a random intecept AND slope for each group
#summary(m2_2)

#Q3: Does Treatment have a different effect on DBH-ratio in bins?
# to see patterns more easily make bins reducing the data to six bins only
B1<-hist(Grw.pos$DBH,10,plot=F) #could also use "cut(xxx
#I did it the pedestrian way
#TREATMENT Group, no colnames for ease of later merging
Grw.pos1<-filter(ScalGrw2,Growth.rate>1 & group==1)
sb1_1<-Grw.pos1[which(Grw.pos1$DBH>0 & Grw.pos1$DBH<=5),]
sb1_1<-cbind(rep(1,dim(sb1_1)[1]),sb1_1);colnames(sb1_1)<-NULL#add a counter for the bin
sb1_2<-Grw.pos1[which(Grw.pos1$DBH>5.1 & Grw.pos1$DBH<=10),]
sb1_2<-cbind(rep(2,dim(sb1_2)[1]),sb1_2);colnames(sb1_2)<-NULL
sb1_3<-Grw.pos1[which(Grw.pos1$DBH>10.1 & Grw.pos1$DBH<=15),]
sb1_3<-cbind(rep(3,dim(sb1_3)[1]),sb1_3);colnames(sb1_3)<-NULL
sb1_4<-Grw.pos1[which(Grw.pos1$DBH>15.1 & Grw.pos1$DBH<=20),]
sb1_4<-cbind(rep(4,dim(sb1_4)[1]),sb1_4);colnames(sb1_4)<-NULL
sb1_5<-Grw.pos1[which(Grw.pos1$DBH>20.1 & Grw.pos1$DBH<=25),]
sb1_5<-cbind(rep(5,dim(sb1_5)[1]),sb1_5);colnames(sb1_5)<-NULL
sb1_6<-Grw.pos1[which(Grw.pos1$DBH>25.1 & Grw.pos1$DBH<=30),]
sb1_6<-cbind(rep(6,dim(sb1_6)[1]),sb1_6);colnames(sb1_6)<-NULL
means1<-c(mean(sb1_1[,7]),mean(sb1_2[,7]),mean(sb1_3[,7]),mean(sb1_4[,7]),mean(sb1_5[,7]),mean(sb1_6[,7]))
medians1<-c(median(sb1_1[,7]),median(sb1_2[,7]),median(sb1_3[,7]),median(sb1_4[,7]),median(sb1_5[,7]),median(sb1_6[,7]))
stdevs1<-c(sd(sb1_1[,7]),sd(sb1_2[,7]),sd(sb1_3[,7]),sd(sb1_4[,7]),sd(sb1_5[,7]),sd(sb1_6[,7]))
#CONTROL Group
Grw.pos2<-filter(ScalGrw2,Growth.rate>1 & group==2)
sb2_1<-Grw.pos2[which(Grw.pos2$DBH>0 & Grw.pos2$DBH<=5),]
sb2_1<-cbind(rep(1,dim(sb2_1)[1]),sb2_1);colnames(sb2_1)<-NULL#add a counter for the bin
sb2_2<-Grw.pos2[which(Grw.pos2$DBH>5 & Grw.pos2$DBH<=10),]
sb2_2<-cbind(rep(2,dim(sb2_2)[1]),sb2_2);colnames(sb2_2)<-NULL
sb2_3<-Grw.pos2[which(Grw.pos2$DBH>10 & Grw.pos2$DBH<=15),]
sb2_3<-cbind(rep(3,dim(sb2_3)[1]),sb2_3);colnames(sb2_3)<-NULL
sb2_4<-Grw.pos2[which(Grw.pos2$DBH>15 & Grw.pos2$DBH<=20),]
sb2_4<-cbind(rep(4,dim(sb2_4)[1]),sb2_4);colnames(sb2_4)<-NULL
sb2_5<-Grw.pos2[which(Grw.pos2$DBH>20 & Grw.pos2$DBH<=25),]
sb2_5<-cbind(rep(5,dim(sb2_5)[1]),sb2_5);colnames(sb2_5)<-NULL
sb2_6<-Grw.pos2[which(Grw.pos2$DBH>25 & Grw.pos2$DBH<=30),]
sb2_6<-cbind(rep(6,dim(sb2_6)[1]),sb2_6);colnames(sb2_6)<-NULL
means2<-c(mean(sb2_1[,7]),mean(sb2_2[,7]),mean(sb2_3[,7]),mean(sb2_4[,7]),mean(sb2_5[,7]),mean(sb2_6[,7]))
medians2<-c(median(sb2_1[,7]),median(sb2_2[,7]),median(sb2_3[,7]),median(sb2_4[,7]),median(sb2_5[,7]),median(sb2_6[,7]))
stdevs2<-c(sd(sb2_1[,7]),sd(sb2_2[,7]),sd(sb2_3[,7]),sd(sb2_4[,7]),sd(sb2_5[,7]),sd(sb2_6[,7]))

Ancv<-cbind(c(means1,means2),rep(c(1:6),2),c(rep(1,6),rep(2,6)))
colnames(Ancv)<-c("gr.means","bin","Treatment")
Ancv<-as.data.frame(Ancv)

Binned_T<-rbind(as.matrix(sb1_1),as.matrix(sb1_2),as.matrix(sb1_3),as.matrix(sb1_4),as.matrix(sb1_5),as.matrix(sb1_6))
Binned_C<-rbind(as.matrix(sb2_1),as.matrix(sb2_2),as.matrix(sb2_3),as.matrix(sb2_4),as.matrix(sb2_5),as.matrix(sb2_6))
Binned_all<-as.data.frame(rbind(Binned_C,Binned_T))
colnames(Binned_all)<-c("bin","group","plot","tree","year","DBH","growth.rate")
sapply(Binned_all, class)
#for some reason, all are factors, so reassign DHB and growth.rate as numeric
Binned_all$DBH<-as.numeric(as.character(Binned_all$DBH))
Binned_all$growth.rate<-as.numeric(as.character(Binned_all$growth.rate))

#explore graphically what is going on here
ScalGrwT<-filter(ScalGrw2,group==1)
ScalGrwC<-filter(ScalGrw2,group==2)
# how much zero, fast,negative growth?
ZeroT<-length(which(ScalGrwT$Growth.rate==1))#a ratio of 1=Zero growth
ZeroC<-length(which(ScalGrwC$Growth.rate==1))
NegT<-length(which(ScalGrwT$Growth.rate<1))# a ratio <1 = shrinkage
NegC<-length(which(ScalGrwC$Growth.rate<1))
FastT<-length(which(ScalGrwT$Growth.rate<1.2))# ratio >1.2 is fast grwth
FastC<-length(which(ScalGrwC$Growth.rate<1.2))

dat1<-c(1,1,1,2,2,2)
dat2<-c("Zero","Negative","Fast","Zero","Negative","Fast")
dat3<-c(ZeroT,NegT,FastT,ZeroC,NegC,FastC)
dat<-cbind(dat1,dat2,dat3)
colnames(dat)<-c("Group","Growth","Data")
dat<-as.data.frame(dat)
dat$Data<-as.numeric(as.character(dat$Data)) #Note::as above, line 284

#============> MAIN TEXT, Fig.5A&B<==============================================================
#don't ggplot it, else it's a pain with Fig.A
windows(15,6)
par(mfrow=c(1,2))
par(mar=c(6,5,3,1))
plot(means1,xlab="5cm DBH size bins", ylab="Mean DBH-ratio of growing trees", pch=21, col="black",bg=rgb(0,0.6,0.5),cex=1,ylim=c(1,1.3),lwd=2,xaxt="n")
axis(1,at=seq(1,6,1), labels=c("SC1\n0-5cm","SC2\n5.1-10cm","SC3\n10.1-15cm","SC4\n15.1-20cm","SC5\n20.1-25cm","SC6\n25.1-30cm"),tick=F)
arrows(c(1:6),means1-stdevs1,c(1:6),means1+stdevs1,code=3,angle=90,length=0.1,col="green4",lwd=2)
par(new=T)
plot(means2,pch=21, xlab="", ylab="", col="black",bg=rgb(0.9,0.6,0),cex=1,ylim=c(1,1.3), lwd=1,xaxt="n")
arrows(c(1:6),means2-stdevs2,c(1:6),means2+stdevs2,code=3,angle=90,length=0.1,col=rgb(0.9,0.6,0))
legend(4.55,1.30,legend=c(expression(italic("R. niveus")*" removed"), expression(italic("R. niveus")*" present")),pch=c(16,16),col=c(rgb(0,0.6,0.5),rgb(0.9,0.6,0)))
mtext('A',side=3,line=1,at=0.2,cex=2)
#more fast growers in Treatment....but were there in general more trees? Change to percentage of total
plodat<-dat3
pls1<-sum(plodat[c(1,3,5)])
pls2<-sum(plodat[c(2,4,6)])
plotdat<-c(1303/pls1,193/pls1,93/pls1,857/pls2,166/pls2,45/pls2)#the values are rearranged from plodat 
par(mar=c(6,5,3,1))
barplot(plotdat,col=c(rgb(0,0.6,0.5),rgb(0.9,0.6,0)),ylim=c(0,1),space=c(0.2,0,0.2,0,0.2,0),ylab="Proportion of occurrence")
legend(4.7,0.9, legend=c(expression(italic("R. niveus")*" removed"), expression(italic("R. niveus")*" present")),fill=c(rgb(0,0.6,0.5),rgb(0.9,0.6,0)))
box(lty="solid")
mtext('B',side=3,line=1,at=-0.9,cex=2)
text(c(1.3,3.5,5.5),par("usr")[3]-0.02,srt=45,adj=1,labels=c("Fast Growth","Negative Growth","Zero Growth"),xpd=T)

# these graphics suggest that
# - growth is higher only in C in smallest size-class
# - there is more fast growth and less negative growth in T (but non-significant, see below)

#now test this second statement
chidat<-matrix(c(93,640,1750,45,166,875),nrow=2,byrow=T)
#same data as proportions*100 -> there are more points in Treatment than Control!!
#chidat<-matrix(c(0.05,0.102,0.84,0.044,0.15,0.8),nrow=2,byrow=T)
chidat<-matrix(c(50,102,840,44,150,800),nrow=2,byrow=T)
chisq.test(chidat)

#what is going on with smallest size bin? 
#Full Binned dataset including smallest bin
m4<-lmer(growth.rate~group/bin+(1|year/plot),REML=FALSE,data=Binned_all)#no REML to allow Anova comparison
summary(m4) # APPENDIX 1, TABLE S-7<=========
#Reduced binned dataset; effect still strongest in the smallest bin, now bin 2 
Binned2<-filter(Binned_all,bin!=1)
m4<-lmer(growth.rate~group/bin+(1|year/plot),REML=FALSE,data=Binned2)#no REML to allow Anova comparison

# ===============PART 3 MORTALITIES==================================================
# we want  absolute and relative mortalities
# for this we look at those who survive and plot up their size distributions per year
#propertional mortality every year= number that remain/over number at previous year
Surv1<-filter(ScalGrw,ScalGrw$DBH>0) #this removes all NA in Grw
Surv<-group_by(Surv1,group,year)
S<-summarize(Surv,survivors=length(year)) #this is the table of all survivors per year and treatment

All<-group_by(ScalGrw,group,year)
A<-summarize(All,totals=length(year))
#mean annual mortality
Mort_T<-1-(S$survivors[2:8]/S$survivors[1:7])
Mort_C<-1-(S$survivors[10:16]/S$survivors[9:15])

#calculate how many died every year and do the stats
# these must be GLMs with Gamma error structure
S.T<-S[1:8,]#T=Treatment
S.RC<-S[9:16,]#RC=RubusControl
Abs.deaths.p.y.T<-S.T$survivors[-length(S.T$survivors)]-S.T$survivors[-1]
Abs.deaths.p.y.RC<-S.RC$survivors[-length(S.RC$survivors)]-S.RC$survivors[-1]
#test for differences between the Treatment/Control
#comparing absolute depths makes only marginal sense - better to compare relative because of varuable N
Ab_Surv<-glm(Abs.deaths.p.y.T~Abs.deaths.p.y.RC, family=poisson)#counts have poisson error structure
summary(Ab_Surv)
#test for differences
Rel.survivors.p.y.T<-S.T$survivors[-1]/S.T$survivors[-length(S.T$survivors)]
Rel.deaths.p.y.T<-1-Rel.survivors.p.y.T
Rel.survivors.p.y.RC<-S.RC$survivors[-1]/S.RC$survivors[-length(S.RC$survivors)]
Rel.deaths.p.y.RC<-1-Rel.survivors.p.y.RC
# these are proportions, so the glm error structure is binomial
Rel_Dead1<-glm(Rel.deaths.p.y.T~Rel.deaths.p.y.RC, family=binomial)
summary(Rel_Dead1) #APPENDIX 1, TABLE S-9
Rel_Dead2<-glm(Rel.deaths.p.y.T[c(1,3,4,5)]~Rel.deaths.p.y.RC[c(1,3,4,5)],family=binomial)
summary(Rel_Dead2) #APPENDIX 1, TABLE S-10

#while glm is preferrable, a Welch's t-test should do the same trick
t.test(Rel.survivors.p.y.T,Rel.survivors.p.y.RC)#<- not significant


#======> MAIN TEXT, FIGURE 6 
#proportion of deaths each year
Surv_T<-rep(1,7)-S$survivors[2:8]/S$survivors[1:7]
Surv_C<-rep(1,7)-S$survivors[10:16]/S$survivors[9:15]
Year<-c("2014-15","2015-16","2016-17","2017-18","2018-19","2019-20","2020-21")
Yr<-rep(Year,2)
Grp<-c(rep(1,7),rep(2,7))
Srv<-rbind(Grp,Yr,cbind(t(Surv_T),t(Surv_C)) )
Srv<-as.data.frame(t(Srv))
colnames(Srv)<-c("Group","Year","PropSurvivors")
Srv$PropSurvivors<-as.numeric(as.character(Srv$PropSurvivors))

windows(10,4)
p4<- ggplot(data=Srv,aes(x=Year,y=PropSurvivors,fill=Group))+geom_col(position="dodge",color="black")+
  theme_bw()+ labs(y="Annual proportion of deceased trees (all plots)",x="Interval during which mortality occurred")+
  scale_fill_manual(values=c(rgb(0,0.6,0.5),rgb(0.9,0.6,0)),labels=c(c(expression(italic("R. niveus")*" removed")),c(expression(italic("R. niveus")*" present"))))
p4

S<-as.data.frame(cbind(as.factor(Srv$Group),Srv$PropSurvivors))
t.test(S$V2~S$V1)

#Now look at what size most trees are lost, use the wide data for this SCAL
#first count how many have died, that is seen by turning column Mar19 into 0/1
LD<-SCAL$Mar19
LD[LD>0]<-1
table(complete.cases(subset(LD,group==1)),subset(LD,group==1))#TREATMENT
table(complete.cases(subset(LD,group==2)),subset(LD,group==2))#Rubus control
#the above suggests that more trees died in the control group. Which was the worst year?
#get this from the histograms

#now look at mortality pattern. Which size-classes have died per year
#there is probably a smarter way around massaging the data
SCAL.T<-subset(SCAL,group==1)
SCAL.RC<-subset(SCAL,group==2)
#D.Year calulates how many died
D15.T<-SCAL.T[which(SCAL.T$Feb15==0),]
D16.T<-SCAL.T[which(SCAL.T$Feb16==0 & SCAL.T$Feb15!=0),]
D17.T<-SCAL.T[which(SCAL.T$Feb17==0 & SCAL.T$Feb16!=0),]
D18.T<-SCAL.T[which(SCAL.T$Feb18==0 & SCAL.T$Feb17!=0),]
D19.T<-SCAL.T[which(SCAL.T$Mar19==0 & SCAL.T$Feb18!=0),]
D20.T<-SCAL.T[which(SCAL.T$Mar20==0 & SCAL.T$Mar19!=0),]
D21.T<-SCAL.T[which(SCAL.T$Mar21==0 & SCAL.T$Mar20!=0),]
D15.RC<-SCAL.RC[which(SCAL.RC$Feb15==0),]
D16.RC<-SCAL.RC[which(SCAL.RC$Feb16==0 & SCAL.RC$Feb15!=0),]
D17.RC<-SCAL.RC[which(SCAL.RC$Feb17==0 & SCAL.RC$Feb16!=0),]
D18.RC<-SCAL.RC[which(SCAL.RC$Feb18==0 & SCAL.RC$Feb17!=0),]
D19.RC<-SCAL.RC[which(SCAL.RC$Mar19==0 & SCAL.RC$Feb18!=0),]
D20.RC<-SCAL.RC[which(SCAL.RC$Mar20==0 & SCAL.RC$Mar19!=0),]
D21.RC<-SCAL.RC[which(SCAL.RC$Mar21==0 & SCAL.RC$Mar20!=0),]

#======APPENDIX 1, FIGURE S2
windows(12,12)
par(mfrow=c(4,4))
h1.T<-hist(D15.T$Feb14,breaks=c(0,5,10,15,20,25,30,35),col=rgb(0,0.6,0.5),ylim=c(0,25),xlab="size classes",main="DBH of trees that died 2014-2015")
h2.T<-hist(D16.T$Feb15,breaks=c(0,5,10,15,20,25,30,35),col=rgb(0,0.6,0.5),ylim=c(0,25),main="DBH of trees that died 2015-2016",xlab="size classes")
h3.T<-hist(D17.T$Feb16,breaks=c(0,5,10,15,20,25,30,35),col=rgb(0,0.6,0.5),ylim=c(0,25),main="DBH of trees that died 2016-2017",xlab="size classes")
h4.T<-hist(D18.T$Feb17,breaks=c(0,5,10,15,20,25,30,35),col=rgb(0,0.6,0.5),ylim=c(0,25),main="DBH of trees that died 2017-2018",xlab="size classes")
h5.T<-hist(D19.T$Feb18,breaks=c(0,5,10,15,20,25,30,35),col=rgb(0,0.6,0.5),ylim=c(0,25),main="DBH of trees that died 2018-2019",xlab="size classes")
h6.T<-hist(D20.T$Mar19,breaks=c(0,5,10,15,20,25,30,35),col=rgb(0,0.6,0.5),ylim=c(0,25),main="DBH of trees that died 2019-2020",xlab="size classes")
h7.T<-hist(D21.T$Mar20,breaks=c(0,5,10,15,20,25,30,35),col=rgb(0,0.6,0.5),ylim=c(0,25),main="DBH of trees that died 2020-2021",xlab="size classes")
h1.RC<-hist(D15.RC$Feb14,breaks=c(0,5,10,15,20,25,30,35),col=rgb(0.9,0.6,0),ylim=c(0,25),main="DBH of trees that died 2014-2015",xlab="size classes")
h2.RC<-hist(D16.RC$Feb15,breaks=c(0,5,10,15,20,25,30,35),col=rgb(0.9,0.6,0),ylim=c(0,25),main="DBH of trees that died 2015-2016",xlab="size classes")
h3.RC<-hist(D17.RC$Feb16,breaks=c(0,5,10,15,20,25,30,35),col=rgb(0.9,0.6,0),ylim=c(0,25),main="DBH of trees that died 2016-2017",xlab="size classes")
h4.RC<-hist(D18.RC$Feb17,breaks=c(0,5,10,15,20,25,30,35),col=rgb(0.9,0.6,0),ylim=c(0,25),main="DBH of trees that died 2017-2018",xlab="size classes")
h5.RC<-hist(D19.RC$Feb18,breaks=c(0,5,10,15,20,25,30,35),col=rgb(0.9,0.6,0),ylim=c(0,25),main="DBH of trees that died 2018-2019",xlab="size classes")
h6.RC<-hist(D20.RC$Mar19,breaks=c(0,5,10,15,20,25,30,35),col=rgb(0.9,0.6,0),ylim=c(0,25),main="DBH of trees that died 2019-2020",xlab="size classes")
h7.RC<-hist(D21.RC$Mar20,breaks=c(0,5,10,15,20,25,30,35),col=rgb(0.9,0.6,0),ylim=c(0,25),main="DBH of trees that died 2020-2021",xlab="size classes")


#===> MAIN TEXT, FIGURE 7 using a smarter way than what I did above:
DBH.pos<-filter(ScalGrw,DBH>0)
p1<-ggplot(data=filter(DBH.pos, group==1),aes(x=DBH))+
  geom_histogram(binwidth=5,fill=rgb(0,0.6,0.5),color="black",boundary=0)+ylim(0,120)+xlim(0,35)+
  facet_wrap(~year,nrow=2)+ theme(axis.text.x=element_text(angle=90))+
  labs(y="Count",x="size-classes")+theme_bw()
p2<-ggplot(data=filter(DBH.pos, group==2),aes(x=DBH))+
  geom_histogram(binwidth=5,fill=rgb(0.9,0.6,0),color="black",boundary=0)+ylim(0,120)+
  facet_wrap(~year,nrow=2)+ theme(axis.text.x=element_text(angle=90))+
  labs(y="Count",x="size-classes")+theme_bw()
windows(10,10)
ggarrange(p1,p2,labels=c("A","B"),ncol=1,nrow=2)

#now extract the data from the ggplot
dat1<-ggplot_build(p1)
dat2<-ggplot_build(p2)
# these are now the extracted frequency distributions, read down the column
freqs1<-matrix(dat1$data[[1]][,2],nrow=7)
freqs2<-matrix(dat2$data[[1]][,2],nrow=7)

chisq.test(freqs1[1,],freqs2[1,])

1-pchisq(sum((freqs1[1,]-freqs2[1,])^2/freqs1[1,]),6)
  
#now test this
chidat<-matrix(c(107,120,1030,118,135,706),nrow=2,byrow=T)
chisq.test(chidat)

#obsolete, was turned into a table
#============PART 4: RECRUITMENT<==========================
rec2<-read.table("Scalesia_T_allrecruits.txt",header=T)
rec3<-read.table("Scalesia_C_allrecruits.txt",header=T)

#Now gather the two sizes into a single column 
Rec2<-rec2%>%
  gather('Seedling','Sapling','Tree1','Tree2','Tree3','Tree4',key="Rec.Size",value="Number")
Rec3<-rec3%>%
  gather('Seedling','Sapling','Tree1','Tree2','Tree3','Tree4',key="Rec.Size",value="Number")
#========> OBSOLETE FIGURE, but nice...is now a table <==============================================
p7a<- ggplot(data=Rec2, aes(x=Rec.Size,y=Number,fill=Rec.Size))+geom_col(color="black")+
  facet_wrap(~Year,nrow=1)+theme_bw()+ labs(y="Total count of recruits (all TREATMENT plots)",x="Recruit size-classes")+
  scale_x_discrete(limits=c("Seedling","Sapling","Tree1","Tree2","Tree3","Tree4"))+
  theme(axis.text.x=element_text(angle=30,hjust=0.9,vjust=0.85))+
  #scale_fill_grey(start=.1,end=.7,limits=c("Seedling","Sapling","Tree1","Tree2","Tree3","Tree4"))+
  scale_fill_brewer(palette="Greens",limits=c("Seedling","Sapling","Tree1","Tree2","Tree3","Tree4"))+
  guides(fill=F)
p7b<-ggplot(data=Rec3, aes(x=Rec.Size,y=Number,fill=Rec.Size))+geom_col(color="black")+
  facet_wrap(~Year,nrow=1)+theme_bw()+ labs(y="Total count of recruits (all CONTROL plots)",x="Recruit size-classes")+
  scale_x_discrete(limits=c("Seedling","Sapling","Tree1","Tree2","Tree3","Tree4"))+
  theme(axis.text.x=element_text(angle=30,hjust=0.9,vjust=0.85))+ylim(c(0,200))+
  #scale_fill_grey(start=.1,end=.7,limits=c("Seedling","Sapling","Tree1","Tree2","Tree3","Tree4"))+
  scale_fill_brewer(palette="YlOrRd",limits=c("Seedling","Sapling","Tree1","Tree2","Tree3","Tree4"))+
  guides(fill=F)
windows(25,15)
ggarrange(p7a,p7b,labels=c("A","B"),ncol=1,nrow=2)

#these survival data for the exponential survival decay were calculated in an excel file
#"Natural Scalesia recruitment Controlled Area for Bernhard.xlsx"
surv_dat=as.data.frame(cbind(seq(1:6),c(1,0.3,0.1,0.2,0.2,0.02)))
colnames(surv_dat)<-c("size.class","frequency")

# ==========> MAIN TEXT, FIGURE 8<=====================================
#Note: the survivor function is exp(-t/nu), where a fraction 1/nu dies in the first interval amd subsequently (Crawley 2013, p 874)
# for the data here it calculates to 1.14
p1<-ggplot(surv_dat,aes(x=size.class,y=frequency))+geom_point(shape=21,fill="green4",size=5)+
  theme_bw()+
  geom_smooth(fullrange=T,method="lm", formula=(y~exp(-x/1.14)),se=T,level=0.95,color=1)+
  scale_x_continuous(breaks=c(1:6),labels=c("1"="Seedling\n0-200cm","2"="Sapling\n201-400cm","3"="Tree1\n401-500cm","4"="Tree2\n501-600cm","5"="Tree3\n601-700cm","6"="Tree4\n701-800cm"))+
  labs(y="Prop. frequency of seedling stage = Proxy for rel. survival", x="Size class")+
  #scale_x_discrete(labels=labls)
windows(10,5)
p1

# look at the time-lag towards extinction
#use Crowley's survivor function (exp(-t/mu)) to explore how long a population without recruits will last, mu=fraction of those who die 1/deaths, so 
# this is the wellknown N(t)-N(0)*exp(rt), only expressed as exp(t*-proportion of deaths)
windows(10,5)
#plot(sapply(seq(0,30,by=1),function(i) exp(i*-0.25)),col="blue",bg="blue",xlab="DBH in cm",ylab="Proportional survival",xlim=c(0,30),ylim=c(0,1))
#text(3,0.4, "Hamann (2001)", col="blue")
#lines(1:20,function(i) exp(i*-0.25))
#par(new=T)

#the instantaneous death rate is log(finite survival rate); we should not use the finite death rate here
dr1<-log(0.875)#R. niveus removed 1-0.125, which is discrete death rate
dr2<-log(0.838)#R. niveus present 1-0.162
#this figure not used in the paper
windows(10,5)
plot(sapply(seq(0,30,by=1),function(i) exp(i*dr1)),type="l",col=rgb(0,0.6,0.5),bg=rgb(0,0.6,0.5),xlab="",ylab="Survivors",xlim=c(0,30),ylim=c(0,1))
text(7,0.8,expression(paste(italic('R. niveus')," removed - 2020")),col=rgb(0,0.6,0.5))
par(new=T)
plot(sapply(seq(0,40,by=1),function(i) exp(i*dr2)),type="l",col=rgb(0.9,0.6,0),bg=rgb(0.9,0.6,0),xlab="Years",ylab="",xlim=c(0,30),ylim=c(0,1))
text(6,0.2,expression(paste(italic('R. niveus')," present - 2020")),col=rgb(0.9,0.6,0))
