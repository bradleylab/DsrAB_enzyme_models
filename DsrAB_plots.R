plot(f,alphaT.34.a, ylim=c(0.988,1.008), xlim=c(0,1))
plot(f,alpha.ox.34.a,ylim=c(1.000,1.016), xlim=c(0,1))
plot(f,alpha.red.34.a,ylim=c(0.982,0.992), xlim=c(0,1))
d34redvf=lm(alpha.red.34~f)
d34oxavf = lm(alpha.ox.34.a~f)
d34Tavf = lm(alphaT.34.a[2:33]~f[2:33])


plot(f,alphaT.34,ylim=c(0.988,1.008), xlim=c(0,1))
plot(f,alpha.ox.34,ylim=c(1.000,1.016), xlim=c(0,1))
plot(f,alpha.red.34,ylim=c(0.982,0.992), xlim=c(0,1))
d34redavf=lm(alpha.red.34.a~f)
d34oxvf = lm(alpha.ox.34~f)
d34Tvf = lm(alphaT.34~f)


realj = subset(j,j>0.0001)
JvF <-cbind(f,j)
rJvF = subset(JvF, j>0.0001)
plot(rJvF[,1],rJvF[,2], xlab="f",ylab="j", xlim=c(0,1),ylim=c(0,1))
Mjvf = lm(rJvF[,2]~rJvF[,1])

cs=Mjvf$coefficients
abline(cs[1],cs[2])


plot(f,Dsr$d34red)
lnf = log(f)


plot(lnf,Dsr$ThiosulfS,ylim=c(0,10), xlab="ln(f)", ylab="Thiosulfate sulfur")
ThioM = lm(Dsr$ThiosulfS~lnf)
cs=ThioM$coefficients
abline(cs[1],cs[2])


plot(lnf,Dsr$TrithioS, ylim=c(0,10), xlab="ln(f)",ylab="Trithionate sulfur")
TrithM = lm(Dsr$TrithioS~lnf)
cs=TrithM$coefficients
abline(cs[1],cs[2])


#compare slopes
library(lmtest)
coxtest(d34redvf,d34redavf)
jtest(d34redvf,d34redavf)

library(nlme)
lme

#####    ANCOVA of Reduced moiety    #####

RedMod<-data.frame(RResults$f)
RedMod<-cbind(RedMod, RResults$alpha.red.34)
names(RedMod)[1] <- "f"
names(RedMod)[2] <- "alpha.red.34"
RedMod<-cbind(RedMod, 1)  #1 = model uses product
names(RedMod)[3] <- "model"

RedMod.a<-data.frame(RResults.a$f)
RedMod.a<-cbind(RedMod.a, RResults.a$alpha.red.34)
names(RedMod.a)[1] <- "f"
names(RedMod.a)[2] <- "alpha.red.34"
RedMod.a<-cbind(RedMod.a, 0)  #0 = model uses reactant
names(RedMod.a)[3] <- "model"

RedAll<-rbind(RedMod,RedMod.a)
plot(RedAll$f,RedAll$alpha.red.34,ylim=c(0.982,0.992), xlim=c(0,1))

by(RedAll,RedAll$model,summary)

RedLM = lm(RedAll$alpha.red.34~RedAll$f+RedAll$model+RedAll$model:RedAll$f)
anova(RedLM)


#####    ANCOVA of Oxidized moiety    #####
OxMod<-data.frame(RResults$f)
OxMod<-cbind(OxMod, RResults$alpha.ox.34)
names(OxMod)[1] <- "f"
names(OxMod)[2] <- "alpha.ox.34"
OxMod<-cbind(OxMod, 1)  #1 = model uses product
names(OxMod)[3] <- "model"

OxMod.a<-data.frame(RResults.a$f)
OxMod.a<-cbind(OxMod.a, RResults.a$alpha.ox.34)
names(OxMod.a)[1] <- "f"
names(OxMod.a)[2] <- "alpha.ox.34"
OxMod.a<-cbind(OxMod.a, 0)  #0 = model uses reactant
names(OxMod.a)[3] <- "model"

OxAll<-rbind(OxMod,OxMod.a)
plot(OxAll$f,OxAll$alpha.ox.34,ylim=c(1.000,1.016), xlim=c(0,1))

by(OxAll,OxAll$model,summary)

OxLM = lm(OxAll$alpha.ox.34~OxAll$f+OxAll$model+OxAll$model:OxAll$f)
anova(OxLM)


#####    ANCOVA of Total fractionation    #####
TotMod<-data.frame(RResults$f)
TotMod<-cbind(TotMod, RResults$alphaT.34)
names(TotMod)[1] <- "f"
names(TotMod)[2] <- "alphaT.34"
TotMod<-cbind(TotMod, 1)  #1 = model uses product
names(TotMod)[3] <- "model"

TotMod.a<-data.frame(RResults.a$f)
TotMod.a<-cbind(TotMod.a, RResults.a$alphaT.34)
names(TotMod.a)[1] <- "f"
names(TotMod.a)[2] <- "alphaT.34"
TotMod.a<-cbind(TotMod.a, 0)  #0 = model uses reactant
names(TotMod.a)[3] <- "model"

TotAll<-rbind(TotMod,TotMod.a)
plot(TotAll$f,TotAll$alphaT.34,ylim=c(0.988,1.008), xlim=c(0,1))

by(TotAll,TotAll$model,summary)
la
TotLM = lm(TotAll[-34,]$alphaT.34~TotAll[-34,]$f+TotAll[-34,]$model+TotAll[-34,]$model:TotAll[-34,]$f)
anova(TotLM)





###############################
###    Mixing of isotope effects ####

MixMod<-data.frame(Dsr.34Rs$SO3)
MixMod<-cbind(MixMod, Dsr.34Rs$ox)
MixMod<-cbind(MixMod, Dsr.34Rs$red)
names(MixMod)[1] <- "Rso3"
names(MixMod)[2] <- "Rox"
names(MixMod)[3] <- "Rred"

Mix.alpha = 0.98467919

MixMod.complete = MixMod[complete.cases(MixMod),]
Xs <- matrix(seq(0.10,0.99,by=0.01),nrow=1)
XsM <- Xs[rep(1:1,11),]
MixMod.complete<-cbind(MixMod.complete,XsM)

n=dim(MixMod.complete)[2]
MixMod.complete[,4:n] = (MixMod.complete$Rred-MixMod.complete[,4:n]*MixMod.complete$Rso3*Mix.alpha)/(MixMod.complete$Rox*(1-MixMod.complete[,4:n]))

plot(Xs[1,],MixMod.complete[1,4:n], type="n", ylim=c(0.976,1.0),ylab="alpha-unk",xlab="X")
lines(Xs[1,],MixMod.complete[1,4:n])
lines(Xs[1,],MixMod.complete[2,4:n])
lines(Xs[1,],MixMod.complete[3,4:n])
lines(Xs[1,],MixMod.complete[4,4:n])
lines(Xs[1,],MixMod.complete[5,4:n])
lines(Xs[1,],MixMod.complete[6,4:n])
lines(Xs[1,],MixMod.complete[7,4:n])
lines(Xs[1,],MixMod.complete[8,4:n])
lines(Xs[1,],MixMod.complete[9,4:n])
lines(Xs[1,],MixMod.complete[10,4:n])
lines(Xs[1,],MixMod.complete[11,4:n])

#2Dec Plots for Wil
Dropbox/DataShare_Alex_Wil/Figures.2Dec

plot(Reliable(f),Reliable(Results$alpha.ox.34),xlab="f",ylab="alpha.SO3-sulfonate", xlim=rev(range(f)))
plot(Reliable(f),Reliable(Results$alpha.red.34),xlab="f",ylab="alpha.SO3-reduced.S", xlim=rev(range(f)))
plot(Reliable(f),Reliable(Dsr$d34ox),xlab="f",ylab="d34S.sulfonate", xlim=rev(range(f)))
plot(Reliable(f),Reliable(Dsr$d34red),xlab="f",ylab="d34S.reduced.S", xlim=rev(range(f)))




plot(Xs[1,],MixMod.complete[1,4:n], type="n", ylim=c(0.976,1.1),ylab="alpha.secondary",xlab="X")
lines(Xs[1,],MixMod.complete[1,4:n])  #col=cm.colors(9)[1]
lines(Xs[1,],MixMod.complete[2,4:n])
lines(Xs[1,],MixMod.complete[3,4:n])
lines(Xs[1,],MixMod.complete[4,4:n])
lines(Xs[1,],MixMod.complete[5,4:n])
lines(Xs[1,],MixMod.complete[6,4:n])
#lines 7 and 8 excluded, because these are included in calc of alpha for f near 1
#lines(Xs[1,],MixMod.complete[7,4:n])
#lines(Xs[1,],MixMod.complete[8,4:n])
lines(Xs[1,],MixMod.complete[9,4:n])
lines(Xs[1,],MixMod.complete[10,4:n])
lines(Xs[1,],MixMod.complete[11,4:n])

