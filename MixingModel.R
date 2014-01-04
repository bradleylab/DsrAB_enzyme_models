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

MixMod<-Reliable(MixMod)
MixMod.complete = MixMod[complete.cases(MixMod),]
Xs <- matrix(seq(0.01,0.99,by=0.01),nrow=1)
XsM <- Xs[rep(1:1,length(MixMod.complete[,1])),]
MixMod.complete<-cbind(MixMod.complete,XsM)

f.reliable <-Reliable(f)
f.MM.complete <-f.reliable[complete.cases(MixMod)]

n=dim(MixMod.complete)[2]
MixMod.complete[,4:n] = (MixMod.complete$Rred-MixMod.complete[,4:n]*MixMod.complete$Rso3*Mix.alpha)/(MixMod.complete$Rox*(1-MixMod.complete[,4:n]))

#plot the results
# See Code_for_figures.R