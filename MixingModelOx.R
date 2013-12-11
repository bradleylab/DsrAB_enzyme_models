###############################
###    Mixing of isotope effects ####

##first to set up the dataframe
MixMod<-data.frame(Dsr.34Rs$SO3)
MixMod<-cbind(MixMod, Dsr.34Rs$ox)
MixMod<-cbind(MixMod, Dsr.34Rs$red)
MixMod<-cbind(MixMod,j)
names(MixMod)[1] <- "Rso3"
names(MixMod)[2] <- "Rox"
names(MixMod)[3] <- "Rred"
names(MixMod)[4]<-"MM.j"

## based on product Rayleigh
##Mix.alpha = 1.00298 #1point

##based on reactant Rayleigh
#Mix.alpha = 1.002748

#set Mix.alpha depending on which Rayleigh model you choose
Mix.alpha = 1.00298

MixMod.complete = MixMod[complete.cases(MixMod),]
Xs <- matrix(seq(0.10,0.99,by=0.01),nrow=1)
XsM <- Xs[rep(1:1,11),]
MixMod.complete<-cbind(MixMod.complete,XsM)

n=dim(MixMod.complete)[2]
MixMod.complete[,5:n] = (-5*MixMod.complete$Rox - 3*Mix.alpha*MixMod.complete$Rox*MixMod.complete$Rso3 + 10*(MixMod.complete$Rox^2)* MixMod.complete[,5:n] - 2*Mix.alpha*MixMod.complete$Rox*MixMod.complete$Rso3*MixMod.complete[,5:n] -
5*(MixMod.complete$Rox^2) + 5 *Mix.alpha*MixMod.complete$Rox*MixMod.complete$Rso3*(MixMod.complete[,5:n]^2) - sqrt(3) * sqrt (3*(Mix.alpha^2)*(MixMod.complete$Rox^2)*(MixMod.complete$Rso3^2)+ 
4* (Mix.alpha^2)*(MixMod.complete$Rox^2)*(MixMod.complete$Rso3^2)*MixMod.complete[,5:n] +18*(Mix.alpha^2)*(MixMod.complete$Rox^2)*(MixMod.complete$Rso3^2)*(MixMod.complete[,5:n]^2) - 
60 * (Mix.alpha^2)*(MixMod.complete$Rox^2)*(MixMod.complete$Rso3^2)*(MixMod.complete[,5:n]^3) + 35*(Mix.alpha^2)*(MixMod.complete$Rox^2)*(MixMod.complete$Rso3^2)*(MixMod.complete[,5:n]^4))) / 
(5*((MixMod.complete$Rox^2) - 2*(MixMod.complete$Rox^2)*MixMod.complete[,5:n] + (MixMod.complete$Rox^2)*(MixMod.complete[,5:n]^2)))

#plot the results
plot(Xs[1,],MixMod.complete[1,5:n], type="n",ylab="alpha-unk",xlab="X", col="red")
lines(Xs[1,],MixMod.complete[1,5:n],col="green")  #col=cm.colors(9)[1]
lines(Xs[1,],MixMod.complete[2,5:n],col="green")
lines(Xs[1,],MixMod.complete[3,5:n],col="green")
lines(Xs[1,],MixMod.complete[4,5:n],col="green")
lines(Xs[1,],MixMod.complete[5,5:n],col="green")
lines(Xs[1,],MixMod.complete[6,5:n],col="green")
#lines 7 and 8 excluded, because these are included in calc of alpha for f near 1
#lines(Xs[1,],MixMod.complete[7,5:n])
#lines(Xs[1,],MixMod.complete[8,5:n])
lines(Xs[1,],MixMod.complete[9,5:n],col="green")
lines(Xs[1,],MixMod.complete[10,5:n],col="green")
lines(Xs[1,],MixMod.complete[11,5:n],col="green")

