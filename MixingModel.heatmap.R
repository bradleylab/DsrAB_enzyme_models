###############################
###    Mixing of isotope effects ####

##first to set up the dataframe
MixMod<-data.frame(Dsr.34Rs$SO3)
MixMod<-cbind(MixMod, Dsr.34Rs$ox)
MixMod<-cbind(MixMod, Dsr.34Rs$red)
names(MixMod)[1] <- "Rso3"
names(MixMod)[2] <- "Rox"
names(MixMod)[3] <- "Rred"

## based on product Rayleigh
##Mix.alpha = 0.98467919

##based on reactant Rayleigh
#Mix.alpha = 0.98478247

#set Mix.alpha depending on which Rayleigh model you choose
Mix.alpha = 0.98467919

MixMod.complete = MixMod[complete.cases(MixMod),]
Xs <- matrix(seq(0.10,0.99,by=0.01),nrow=1)
XsM <- Xs[rep(1:1,11),]
MixMod.complete<-cbind(MixMod.complete,XsM)

n=dim(MixMod.complete)[2]
MixMod.complete[,4:n] = (MixMod.complete$Rred-MixMod.complete[,4:n]*MixMod.complete$Rso3*Mix.alpha)/(MixMod.complete$Rox*(1-MixMod.complete[,4:n]))

#plot the results
MixMod.plot=MixMod.complete[-c(7,8),] #lines 7 and 8 excluded, because these are included in calc of alpha for f near 1
MixModf=cbind(MixMod,f)
MixModf.complete = MixModf[complete.cases(MixModf),]
rank(MixModf.complete$f)

plot(Xs[1,],MixMod.plot[1,4:n], type="n", ylim=c(0.976,1.0),ylab="alpha-unk",xlab="X")

for (jj in 1:9) {
lines(Xs[1,],MixMod.plot[jj,4:n], lwd = 2, col=heat.colors(9)[rank(MixModf.complete$f)[jj]])
}


MixModf=cbind(MixMod,f)
MixModf.complete = MixModf[complete.cases(MixModf),]
rank(MixModf.complete$f)