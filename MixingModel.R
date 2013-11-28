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
plot(Xs[1,],MixMod.complete[1,4:n], type="n", ylim=c(0.976,1.0),ylab="alpha-unk",xlab="X")
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

