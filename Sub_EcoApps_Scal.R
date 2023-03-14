# VERSION FOR SUBMISSION
# Scalesia_file; Bernhard Riegl, started on 5/27/2019, first submitted 3/21/21; fixed, updated and cleaned by 10/19/21
# resubmitted 11/4/2021, modified after review 1/10/2022 and checked 3/30/2022
# next change for resubmission 8/9/2022: added DBH analyses, nLMEr and filepath simplification as of 11/23/2022
# accepted paper version 3/10/2023, errors cleaned 3/13/2023;

#Map of the Galapagos for MAIN TEXT FIGURE 1, the final publication figure built by Anna Walentowitz
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
#1=Rubus removed (Treatment), 2=Rubus present (Control)
library(tidyverse)
library(ggpubr)
library(nlme)
library(lme4)
library(MuMIn)

#PART 1: GROWTH MODEL FOR SCALESIA PEDUNCULATA
setwd("D:/University/PAPERS/E-Pacific/Scalesia/Data/")
DBH<-read.table("Scalesia_DBF-TH.txt",header=T) # DBH and Height at random locations in forest prior to beginning of experiment in 2014
DBH_Treat<-read.table("Scalesia-DBH-Height_Treatment.txt",header=T)# BH and Height of haphazardly chosen thees in Treatment plots in 2020 (Rubus removed)
DBH_Control<-read.table("Scalesia-DBH-Height_Control.txt",header=T)# BH and Height of haphazardly chosen thees in Control plots in 2020 (Rubus present)

# calculated for the first DBH dataset, taken from random spots within forest IN 2014
# calculate the parameters for the non-linear model, which is y=a-b*exp(-cx)
# Note: this is standard expression of the asymptotic regression model a-(a-b)exp(-cX) 
# under constraint that b=0, it turns into the negative exponential equation Y=a*(1-exp(-cX))
# the negative exponential equation passes through the origin!
# a=approximate asymptotic max height=1000
a=950
b=820
y=200 #where curve is rising most steeply, height is greatest...next line
x=2 #where x=2
c=(log((a-y)/b))/x
#===>APPENDIX 1, TABLE S-1
model1<-nls(Height_cm ~ a-b*exp(-c*DBH_cm),data=DBH, start=list(a=950,b=820,c=0.12))#3 parameter asymptotic exponential
summary(model1)
#===>APPENDIX 1, TABLE S-2
model2<-nls(Height_cm ~ a*(1-exp(-c*DBH_cm)),data=DBH,start=list(a=950,c=0.12))#2 parameter asymptotic exponential
summary(model2)
anova(model1, model2)

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
#setwd("to wherever you want the png to go")
png("Fig_2.png",width=8.5,height=4.25,units='cm',res=450)
par(mar=c(2,1.8,1.5,0.5)+0.1)
par(cex.lab=0.25,cex=0.7,cex.sub=0.7,cex.axis=0.5)
plot(DBH$DBH_cm,DBH$Height_cm,pch=1,col="blue",bg="blue",cex=0.7,xlim=c(0,30),ylim=c(0,1100),axes=F,frame=T)
#lines(av,bv1,col="blue",lwd=2,lty=2)
lines(av,bv2,col="blue",lwd=2)
text(27.5,1000, "pre-experiment - 2014", col="blue",cex=0.5)
par(new=T)
plot(DBH_Treat$DBH_cm,DBH_Treat$Height_cm,pch=1,col=rgb(0,0.6,0.5),bg=rgb(0,0.6,0.5),cex=0.7,xlim=c(0,30),ylim=c(0,1100),axes=F,frame=T)
lines(av,bv3,col=rgb(0,0.6,0.5),lwd=2)
lines(av,bv4,col=rgb(0,0.6,0.5))
text(27,850,expression(paste(italic('R. niveus')," removed - 2020")),col=rgb(0,0.6,0.5),cex=0.5)
par(new=T)
plot(DBH_Control$DBH_cm,DBH_Control$Height_cm,pch=1,col=rgb(0.9,0.6,0),bg=rgb(0.9,0.6,0),cex=0.7,xlim=c(0,30),ylim=c(0,1100),axes=F,frame=T)
lines(av,bv5,col=rgb(0.9,0.6,0),lwd=2)
lines(av,bv6,col=rgb(0.9,0.6,0))
text(27,700,expression(paste(italic('R. niveus')," present - 2020")),col=rgb(0.9,0.6,0),cex=0.5)
axis(1,seq(0,30,5),labels=c(0,5,10,15,20,25,30),tck=-0.05,mgp=c(1,0.1,0))
axis(2,seq(0,1000,200),labels=c(0,200,400,600,800,1000),tck=-0.03,mgp=c(1,0.4,0),las=1)
title(ylab="Height (cm)",xlab="DBH (cm)",mgp=c(1.1,0.5,0),cex.lab=0.7)
dev.off()

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

# Options (Pinheiro & Bates 2000, p.355): 1) explicitly fit model or 2) use self-starting function (Crawley 2013,p.728; or Pin. & Bates Appendix C)
# first: separate, nonlinear fits by group via nlsList (Pinheiro and Bates 2000, p. 347), we explicitly fit the model
mymodel<-nlsList(Height~a*(1-exp(-c*DBH))|Group,data=HTS,start=c(a=1000,c=0.1))
summary(mymodel)

HTS2<-groupedData(Height~DBH | Phase,data=HTS)
#HTS2<-HTS[-which(HTS2$Height==0),]#if you want to remove the 0-heights
#set up control list for LME
control.list <- lmeControl(maxIter = 500, msMaxIter = 500, msMaxEval=500,tolerance = 0.1, msTol = 0.1, sing.tol=1e-20)
#Description: factor a (asymp size) depends on Group, c depends on nothing. a and c very randomly across Phases (and Groups)
#following advice (https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#testing-significance-of-random-effects) random effect was chosen a priori and kept constant.
#Group has a fixed effect on a, c is a constant, a and c show random variation w. Group
mm2a<-nlme(Height~a*(1-exp(-c*DBH))-1,fixed=list(a~Phase,c~1),random= a+c~1|Group,data=HTS2,start=c(a=950,c=0.12,Group=1), method="ML",control=control.list)
#here, no fixed effect on a or c, but a  random variability from Group to Group on a and c
mm2c<-nlme(Height~a*(1-exp(-c*DBH))-1,fixed=a+c~1,random= a+c~1|Group,data=HTS2,start=c(a=950,c=0.12), method="ML",control=control.list)#chosen model
anova(mm2a,mm2c)#selects mm2c

summary(mm2c)$tTable#in the summary the Fixed effect values are the means of the parameter values (Crawley p.725), so get each with coef(bla), see further below
#below plot shows that the three groups pr-, post-control, post-treatment differ in slope and asymptote
windows(10,5); plot(augPred(mm2c))#==>APPENDIX 1: FIGURE S1<=====

# now play the same with the self-starting function, as a test; it gives the same results as above, except the reported c needs to be exponentiated
mod.lis<-nlsList(Height~SSasympOrig(DBH,a,c)|Group,data=HTS2)
mod.nlme<-nlme(mod.lis)
coef(mod.nlme)# take exp(c) and the results are basically the same as in line 109 and 121
plot(mod.nlme)#this indicates heteroscedasticty, so build a heteroscadestic model
mod2.nlme<-update(mod.nlme,weights=varIdent(form=~1|Phase)) #different variances for each of 2 Phases
mod3.nlme<-update(mod.nlme,weights=varIdent(form=~1|Group)) # different variances for each of 3 experimental groups (2 phases, 2 treatments in Phase 2)
anova(mod.nlme,mod2.nlme,mod3.nlme) #using Phase doesn't help, but with Group, the model improves significantly
anova(mod.nlme,mod3.nlme)#<===this reported in main text, section "Growth Allometry" para 2
coef(mod3.nlme) #===> reported in Table 1<====
exp(coef(mod3.nlme)$c)#===> reported in Table 1<====
#this shows: trees were larger prior to experiment, larger when Rubus removes, smallest with Rubus present
#but it is also clear that basically no recruitment was observed post-2014 (while 2014 had many)
#=>RESULT: above analyses show that across the experiment, tree allometry has changed, and that the effect is greatest if Rubus is present

#================================================================================
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
#=========data for FIGURE 3 MAIN TEXT====================================
#relationship of Growth ratio against DBH is calculated over all data
Scal_Grw<-filter(ScalGrw2, year!="Mar21")
Scal_Grw[is.na(Scal_Grw)]<-0
Scal_Grw[Scal_Grw==Inf]<-0
SG<-Scal_Grw[-which(Scal_Grw$Growth.rate==0),]#eliminate outliers and zeros

#DHB-ratio with diameter as second-order polynomial in main text
# APPENDIX 1, TABLE S-5 and main text Table 2
Gr_mod1<-lm(Growth.rate~poly(DBH,2),data=SG)
summary.lm(Gr_mod1)
summary.aov(Gr_mod1)
#or fit exponential decay
# Select an approximate $\theta$, since theta must be lower than min(y), and greater than zero
theta.0 <- min(SG$Growth.rate) * 0.5  
# Estimate the rest parameters using a linear model. Note: this concerns only the pooled DBH values overall; not year, plot, etc...
model.0 <- lm(log(Growth.rate - theta.0) ~ DBH, data=SG)  
alpha.0 <- exp(coef(model.0)[1])
beta.0 <- coef(model.0)[2]
# Starting parameters
start <- list(alpha = alpha.0, beta = beta.0, theta = theta.0)
start2<-list(alpha = alpha.0, beta = beta.0)
#if error is additive, then it is constant with x-axis and we should not use log-scale
Gr_mod2 <- nls(Growth.rate ~ alpha * exp(beta * DBH) + theta , data = SG, start = start,control = nls.control(maxiter = 1000))
Gr_mod2b<- nls(Growth.rate ~ alpha * exp(beta * DBH), data = SG, start = start2)
anova(Gr_mod2,Gr_mod2b)
summary(Gr_mod2)
#lines(SG_C$DBH,predict(Gr_mod5,list(x=SG_C$DBH)),col='skyblue',lwd=3)
AIC(Gr_mod2,Gr_mod1)

#==> APPENDIX 1, FIGURE S2 <===
#separate evaluation of DBH-ratios in plots with/without R. niveus
ScalGrwT<-filter(ScalGrw2,group==1 & year!="Mar21")#TREATMENT
ScalGrwC<-filter(ScalGrw2,group==2 & year!="Mar21")#WITH RUBUS = CONTROL
#there are many NAs and Infs in there, they make no sense as growth rate
ScalGrwT[is.na(ScalGrwT)]<-0
ScalGrwC[is.na(ScalGrwC)]<-0
ScalGrwT[ScalGrwT==Inf]<-0
ScalGrwC[ScalGrwC==Inf]<-0
SG_T<-ScalGrwT[-which(ScalGrwT$Growth.rate==0),]
SG_C<-ScalGrwC[-which(ScalGrwC$Growth.rate==0),]
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
setwd("E:/University/PAPERS/E-Pacific/Scalesia/Submission_Final/")
png("Fig_3.png",width=8.5,height=6,units='cm',res=450)
par(mar=c(3.6,4,1.5,0.5)+0.1)
par(cex.lab=0.8,cex=0.7,cex.sub=0.7,cex.axis=0.6)
plot(SG[which(SG$group==1),]$DBH,SG[which(SG$group==1),]$Growth.rate,xlab="",ylab="",ylim=c(0.8,1.65),xlim=c(0,30),col=rgb(0,0.6,0.5),cex=0.5,mgp=c(0.1,0.3,0),axes=F,frame=T)
par(new=T)
plot(SG[which(SG$group==2),]$DBH,SG[which(SG$group==2),]$Growth.rate,xlab="",ylab="",ylim=c(0.8,1.65),xlim=c(0,30),col="orange3",cex=0.5,mgp=c(0.1,0.3,0),axes=F,frame=T)
title(ylab="DBH-ratio",mgp=c(1,0.3,0),cex.lab=0.6,xlab="DBH (cm)")
axis(1,seq(0,30,5),labels=c(0,5,10,15,20,25,30),tck=-0.02,mgp=c(1,0.1,0))
axis(2,seq(0.8,1.6,0.1),labels=c(0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6),tck=-0.02,mgp=c(1,0.4,0),las=1)
x<-0:30
y<-predict(Gr_mod2,list(DBH=x))
lines(x,y,col="black",lwd=3)
legend(20,1.6, legend=c(expression(italic("R. niveus")*" removed"), expression(italic("R. niveus")*" present")),pch=21,col=c(rgb(0,0.6,0.5),rgb(0.9,0.6,0)),cex=0.6,bty="n")
box(lty="solid")
dev.off()

#==> MAIN TEXT, FIGURE 4<====
# function for computing mean, DS, max and min values
#note, the correct equation for the 95%SI would be mean+/-SI*t(0.975,n-1)
#setwd("where you want your file")

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
  scale_x_discrete(labels=c("2014","2015","2016","2017","2018","2019","2020"))+
  #facet_wrap(~Plot,nrow=3)+
  theme_bw(base_size=14)+ylim(c(0.7,1.5))+
  theme(axis.text.x=element_text(angle=90,size=14,colour="black"),axis.text.y=element_text(size=14,colour="black"))+
  labs(y="DBH-ratio",x="Sampling year")
p2<-ggplot(data=filter(ScalGrwT, Growth.rate>0.5 & Growth.rate<1.5),aes(x=year,y=Growth.rate))+
  stat_summary(fun.data = min.mean.sd.max, geom = "boxplot")+
  geom_jitter(color=rgb(0,0.6,0.5), position=position_jitter(width=.2), size=1)+
  ggtitle(c(expression(italic("R. niveus")*" removed")))+
  scale_x_discrete(labels=c("2014","2015","2016","2017","2018","2019","2020"))+
  #facet_wrap(~Plot,nrow=3)+
  theme_bw(base_size=14)+ylim(c(0.7,1.5))+
  theme(axis.text.x=element_text(angle=90,size=14,colour="black"),axis.text.y=element_text(size=14,colour="black"))+
  labs(y="DBH-ratio",x="Sampling year")
windows(25,12)
p3<-ggarrange(p2,p1,labels=c("A","B"),font.label=list(size=28,color="black"),ncol=2,nrow=1) 
p3
ggsave("Fig_4.png",p3,width=8.5,height=5,dpi=450)


#Q2: does Treatment have an effect on DBH-ratio overall? At this point not interested in "year".
#log(Growth rate) bc. we have established exponential relationship!
#check best model structure
#random intercept; also: # DBH*Group=>DBH+group+DBH:Group ->two-way interaction
#as dataset use Grw.pos (not SG as in some test versions)
control.list <- lmeControl(maxIter = 500, msMaxIter = 500, msMaxEval=500,tolerance = 0.1, msTol = 0.1, sing.tol=1e-20)
m2_1e<-lme(log(Growth.rate)~1, random=~1|Plot, data=Grw.pos, method="ML")# baseline model, only intercept
m2_1d<-lme(log(Growth.rate)~DBH+group, random=~1|Plot,data=Grw.pos,method="ML")# random intercepts only, based on Plot, no DBH-Group interaction
m2_1c<-lme(log(Growth.rate)~DBH*group, random=~1|Plot,data=Grw.pos,method="ML")# random intercepts only, based on Plot, with DBH-Group interaction
m2_1b<-lme(log(Growth.rate)~DBH+group, random=~1+group|Plot, data=Grw.pos, method="ML",control=control.list)#a random intercept AND slope for each group, no DBH-Group interaction
m2_1a<-lme(log(Growth.rate)~DBH*group, random=~1+group|Plot, data=Grw.pos, method="ML",control=control.list)#full model: a random intercept AND slope for each group, DBH-Group interaction
anova(m2_1a,m2_1b,m2_1c,m2_1d,m2_1e)# this works, because models are ML, with REML anova is not permitted

#Pinheiro and Bates p. 19 "the overall effect of the factor should be assessed with anova, not by t-values of p assoc. with fixed-effect params"
anova(m2_1a) #Appendix: Table S7======
plot(m2_1a,form=resid(.,type="p")~fitted(.)|Plot,abline=0) #[Pinheiro and Bates (2000)p.21]; APPENDIX FIGURE S3=====

#the chosen model is m2_1a; present with REML, for different output tables once with "lme" and once with "lmer"
#Note: using "tree" as the random effect instead of "plot", or "Plot/tree" doesn't make the model fit any better
m2_2<-lme(log(Growth.rate)~DBH*group, random=~1+group|Plot, data=growth, method="REML")#a random intercept AND slope for each group
summary(m2_2)
vars<-c(0.09,0.001,0.01,0.001)
100*vars/sum(vars)
r.squaredGLMM(m2_2)
#also do a lmer to get output and calculate where the variance in the random components is (Crawley, p. 704)
mod2_2aa<-lmer(log(Growth.rate)~DBH+group+(1+group|Plot),REML=T, data=growth)
#mod2_2aa<-lmer(log(Growth.rate)~DBH+group+(1|year),REML=T, data=Grw.pos)#random intercept only
summary(mod2_2aa)
vars<-c(0.11,0.004,0.0054) #these are the fixed effects from REML table
100*vars/sum(vars)
r.squaredGLMM(mod2_2aa)#R2m marginal R2=fixed effects; R2C=conditional, combined fixed and random effects explain X% of variance

#multiple variance structure doesn't help, so abandon
m2_2b<-lme(log(Growth.rate)~DBH*group, random=~1+group|Plot, weights=varIdent(form=~1|Plot),data=Grw.pos, method="ML")
#summary(m2_2b)

#==========Q3: Is there overall more "fast growth", or "Mortality" with R. niveus present/absent?===========================
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

#============> MAIN TEXT, FIGURE 5<==============================================================
#more fast growers in Treatment....but were there in general more trees? Change to percentage of total
plodat<-dat3
pls1<-sum(plodat[c(1,3,5)])
pls2<-sum(plodat[c(2,4,6)])
plotdat<-c(1303/pls1,193/pls1,93/pls1,857/pls2,166/pls2,45/pls2)#the values are rearranged from plodat 

png("Fig_5.png",width=8.5,height=5,units='cm',res=450)
par(mar=c(3.6,4,1.5,0.5)+0.1)
par(cex.lab=0.8,cex=0.7,cex.sub=0.7,cex.axis=0.6)
barplot(plotdat,col=c(rgb(0,0.6,0.5),rgb(0.9,0.6,0)),ylim=c(0,1),space=c(0.2,0,0.2,0,0.2,0),axes=F)
legend(4.7,0.9, legend=c(expression(italic("R. niveus")*" removed"), expression(italic("R. niveus")*" present")),fill=c(rgb(0,0.6,0.5),rgb(0.9,0.6,0)),cex=0.6,bty="n")
box(lty="solid")
text(c(1.3,3.5,5.5),par("usr")[3]-0.02,srt=45,adj=1,labels=c("Fast growth","Negative growth","Zero growth"),xpd=T,cex=0.8)
title(ylab="Proportion of occurrance",mgp=c(1.7,1.2,0))
axis(2,seq(0,1,0.2),labels=c(0,0.2,0.4,0.6,0.8,1),tck=-0.02,mgp=c(1.5,0.4,0),las=1)
dev.off()

#now test what is seen in Fig. 5
chidat<-matrix(c(93,640,1750,45,166,875),nrow=2,byrow=T)
#same data as proportions*100
#chidat<-matrix(c(0.05,0.102,0.84,0.044,0.15,0.8),nrow=2,byrow=T)
chidat<-matrix(c(50,102,840,44,150,800),nrow=2,byrow=T)
chisq.test(chidat)

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
S.T<-S[1:8,]#T=Rubus removed, group 1
S.RC<-S[9:16,]#RC=Rubus present, group 2
Abs.deaths.p.y.T<-S.T$survivors[-length(S.T$survivors)]-S.T$survivors[-1]
Abs.deaths.p.y.RC<-S.RC$survivors[-length(S.RC$survivors)]-S.RC$survivors[-1]
#test for differences between the Rubus removed/Rubus present
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
summary(Rel_Dead1) #APPENDIX 1, TABLE S-8
Rel_Dead2<-glm(Rel.deaths.p.y.T[c(1,3,4,5)]~Rel.deaths.p.y.RC[c(1,3,4,5)],family=binomial)
summary(Rel_Dead2) 

Counter<-c(rep(1,7),rep(2,7))
Surv_dat<-c(Rel.deaths.p.y.T,Rel.deaths.p.y.RC)
Mod<-glm(Surv_dat~Counter,family=quasibinomial)

#while glm is preferrable, a Welch's t-test should do the same trick-> also given in main text
t.test(Rel.survivors.p.y.T,Rel.survivors.p.y.RC)#<- not significant

#======> MAIN TEXT, FIGURE 6 
#proportion of deaths each year,1=Treatment, 2=Rubus Control
Surv_T<-rep(1,7)-S$survivors[2:8]/S$survivors[1:7]
Surv_C<-rep(1,7)-S$survivors[10:16]/S$survivors[9:15]
Year<-c("2014-15","2015-16","2016-17","2017-18","2018-19","2019-20","2020-21")
Yr<-rep(Year,2)
Grp<-c(rep(1,7),rep(2,7))
Srv<-rbind(Grp,Yr,cbind(t(Surv_T),t(Surv_C)) )
Srv<-as.data.frame(t(Srv))
colnames(Srv)<-c("Group","Year","PropSurvivors")
Srv$PropSurvivors<-as.numeric(as.character(Srv$PropSurvivors))

p4<- ggplot(data=Srv,aes(x=Year,y=PropSurvivors,fill=Group))+geom_col(position="dodge",color="black")+
  theme_bw(base_size=14)+ 
  scale_y_continuous(limits=c(0.,0.40),expand=c(0,0))+
  labs(y="Annual proportion of deceased \n trees (all plots)",x="Interval during which mortality occurred")+
  theme(axis.text.x=element_text(size=14,colour="black"),axis.text.y=element_text(size=14,colour="black"))+
  scale_fill_manual(values=c(rgb(0,0.6,0.5),rgb(0.9,0.6,0)),labels=c(c(expression(italic("R. niveus")*" removed")),c(expression(italic("R. niveus")*" present"))))+
  theme(legend.position=c(0.8,0.8),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.title=element_blank())
p4
ggsave("Fig_6.png",p4,width=8.5,dpi=450)

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
p1<-  ggplot(data=filter(DBH.pos, group==1),aes(x=DBH))+
  geom_histogram(binwidth=5,fill=rgb(0,0.6,0.5),color="black",boundary=0)+ylim(0,120)+xlim(0,35)+
  facet_wrap(~year,nrow=2,labeller=labeller(year=c("Feb14"="2014","Feb15"="2015","Feb16"="2016","Feb17"="2017","Feb18"="2018","Mar19"="2019","Mar20"="2020","Mar21"="2021")))+
  theme_bw()+ theme(axis.text.x=element_text(angle=90,size=11,colour="black"),axis.text.y=element_text(size=11,colour="black"))+
  labs(y="Count",x="DBH-size-classes")+theme(axis.title=element_text(size=12))+
  theme(strip.text.x=element_text(size=13,color="black",face="bold"),panel.grid.major=element_blank(),panel.grid.minor=element_blank())
p2<-ggplot(data=filter(DBH.pos, group==2),aes(x=DBH))+
  geom_histogram(binwidth=5,fill=rgb(0.9,0.6,0),color="black",boundary=0)+ylim(0,120)+
  facet_wrap(~year,nrow=2,labeller=labeller(year=c("Feb14"="2014","Feb15"="2015","Feb16"="2016","Feb17"="2017","Feb18"="2018","Mar19"="2019","Mar20"="2020","Mar21"="2021")))+
  theme_bw()+theme(axis.text.x=element_text(angle=90,size=11,colour="black"),axis.text.y=element_text(size=11,colour="black"))+
  labs(y="Count",x="DBH-size-classes")+theme(axis.title=element_text(size=12))+
  theme(strip.text.x=element_text(size=13,color="black",face="bold"),panel.grid.major=element_blank(),panel.grid.minor=element_blank())
windows(10,10)
p5<-ggarrange(p1,p2,labels=c("A","B"),font.label=list(size=20,color="black"),ncol=1,nrow=2)
ggsave("Fig_7.png",p5,width=8.5,dpi=450)
p5

#now extract the data from the ggplot
dat1<-ggplot_build(p1)
dat2<-ggplot_build(p2)
# these are now the extracted frequency distributions, read down the column
freqs1<-matrix(dat1$data[[1]][,2],nrow=7)
freqs2<-matrix(dat2$data[[1]][,2],nrow=7)

chisq.test(freqs1[1,],freqs2[1,])

#=========================RECRUITMENT==============================================
#these survival data for the exponential survival decay were calculated in an outside excel file
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
  labs(y="Prop. frequency of previous height class", x="Height class")+
  theme(axis.text.x=element_text(size=13,colour="black"),axis.text.y=element_text(size=13,colour="black"))+
  theme(axis.title=element_text(size=14))
  #scale_x_discrete(labels=labls)
ggsave("Fig_8.png",p1,width=8.5,dpi=450)
windows(10,5)
p1
# ============= time-lag towards extinction=============
#use Crawley's (2013) survivor function (exp(-t/mu)) to explore how long a population without recruits will last, mu=fraction of those who die 1/deaths, so 
# this is the wellknown N(t)-N(0)*exp(rt), only expressed as exp(t*-proportion of deaths)
