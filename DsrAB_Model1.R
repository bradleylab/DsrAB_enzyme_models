#####
#DsrAB Model 1. ASB for Leavitt et al., ###### 2013.
#####
###################################################################
#FUNCTIONS
##
#function to convert an isotope ratio 34S/32S to a delta34S value
R2Delta <-function(R,CDT)
{
  Delta=((R/CDT)-1)*1000  
  return(Delta)
}

#function to convert a delta34S value to an isotope ratio 34S/32S
Delta2R <-function(Delta, CDT)
{
  Ratio=((Delta/1000)+1)*CDT	
  return(Ratio)
}

#function to convert an alpha value to an epsilon
alpha2epsilon <-function(alpha.in)
{
  epsilon=((alpha.in)-1)*1000  
  return(epsilon)
}

#function to calculate thiosulfate fraction J
CalcJ <-function(ThiosulfS, TrithioS)
{
  J=(ThiosulfS/(ThiosulfS+TrithioS)) 
  return(J)
}

#function to calculate isotope ratio of the total product (Rp)
Calc.Rp <-function(J, Rred, Rox)
{
  Rp=((1/2)*J+(1/3)*(1-J))*Rred + ((1/2)*J+(2/3)*(1-J))*Rox
  return(Rp)
}

#function to calculate isotope mass balance sum for the reduced moieties
Calc.redbal <-function(J,Rediso)
{
  redbal=((1/2)*J+(1/3)*(1-J))*Total.thionate.S*Rediso
  return(redbal)
}

#function to calculate isotope mass balance sum for the oxidized moieties
Calc.oxbal <-function(J,Oxiso)
{
  oxbal=((1/2)*J+(2/3)*(1-J))*Total.thionate.S*Oxiso
  return(oxbal)
}

#function to calculate alpha-total (alphaT)
Calc.alphaT <-function(Rp, Rao, f)
{
  alphaT = log(((Rp/Rao)*(f-1))+1)/(log(f))
  return(alphaT)
}

#function to calculate alpha-total (alphaT)
Calc.alphaT.a <-function(Ra, Rao, f)
{
  alphaT = log(Ra/Rao)/log(f)+1
}

#function to calculate alpha-x (alpha-x where x is reduced or oxidized)
Calc.alphax <-function(Rx, Rao, alphaT, f)
{
  alphax = (Rx/Rao)*(alphaT/((f^alphaT)-1))*(f-1)
  return(alphax)
}

#function to show only reliable results from a calculation. Pass a vector
Reliable <-function(cullresults)
{
  UnReliable <- (is.na(cullresults) | is.nan(cullresults) | Dsr$suspect==TRUE) 
  cullresults[UnReliable]<-NA
  return(cullresults)
}


###### Start of program  ###### 
Dsr<-read.table("/Users/abradley/Documents/Rdata/Wil23Aug.txt",header=TRUE)
attach(Dsr)
#names(Dsr)
numsamples = length(Dsr$ExNo)

cdt = .045005  		#CDT ratio for 34/32. 
cdt3x = 0.007379		#CDT ratio for 3x/32. 

#calculate j,f
j = CalcJ(Dsr$ThiosulfS,Dsr$TrithioS)
f = Dsr$SO3n/Dsr$SO30
Total.thionate.S = (Dsr$ThiosulfS + Dsr$TrithioS)
SO3bydif = Dsr$SO30-Total.thionate.S

#Pull out delta values & convert to Rs
Dsr.34deltas <- Dsr[, c(10:16)]
Dsr.33deltas <- Dsr[,c(17:20)]
Dsr.36deltas <- Dsr[,c(21:24)]
Dsr.34Rs <- Delta2R(Dsr.34deltas,cdt)
colnames(Dsr.34Rs) <- c("SO30", "SO3D","SO3","oxD","ox","redD","red")
Dsr.33Rs <- Delta2R(Dsr.33deltas,cdt3x)
colnames(Dsr.33Rs) <- c("SO30", "SO3", "ox", "red")
Dsr.36Rs <- Delta2R(Dsr.36deltas,cdt3x)
colnames(Dsr.36Rs) <- c("SO30", "SO3", "ox", "red")

############ 34 S ############ 
#Calculate total R of product - using measurements on 253 (default) and on Delta (labeled Delta)
Rp.34 = Calc.Rp(j,Dsr.34Rs$red,Dsr.34Rs$ox)
Rp.34.Delta = Calc.Rp(j,Dsr.34Rs$redD,Dsr.34Rs$oxD)

#calculate alpha total in each case
alphaT.34 = Calc.alphaT(Rp.34, Dsr.34Rs$SO30, f)
alphaT.34.a = Calc.alphaT.a(Rp.34, Dsr.34Rs$SO30, f)
alphaT.34.Delta = Calc.alphaT(Rp.34.Delta,Dsr.34Rs$SO30,f)
alphaT.34.Delta.a = Calc.alphaT.a(Rp.34.Delta,Dsr.34Rs$SO30,f)

#calculate alphas for oxidized and reduced moieties
alpha.red.34 = Calc.alphax(Dsr.34Rs$red, Dsr.34Rs$SO30, alphaT.34, f)
alpha.red.34.Delta = Calc.alphax(Dsr.34Rs$redD, Dsr.34Rs$SO30, alphaT.34.Delta, f)

alpha.red.34.a = Calc.alphax(Dsr.34Rs$red, Dsr.34Rs$SO30, alphaT.34.a, f)
alpha.red.34.Delta.a = Calc.alphax(Dsr.34Rs$redD, Dsr.34Rs$SO30, alphaT.34.Delta.a, f)

alpha.ox.34 = Calc.alphax(Dsr.34Rs$ox, Dsr.34Rs$SO30, alphaT.34, f)
alpha.ox.34.Delta = Calc.alphax(Dsr.34Rs$oxD, Dsr.34Rs$SO30, alphaT.34.Delta, f)

alpha.ox.34.a = Calc.alphax(Dsr.34Rs$ox, Dsr.34Rs$SO30, alphaT.34.a, f)
alpha.ox.34.Delta.a = Calc.alphax(Dsr.34Rs$oxD, Dsr.34Rs$SO30, alphaT.34.Delta.a, f)

############ 33 S ############ 
#Calculate total R of product - using measurements on 253 (default) and on Delta (labeled Delta)
Rp.33 = Calc.Rp(j,Dsr.33Rs$red,Dsr.33Rs$ox)
Rp.33.Delta = Calc.Rp(j,Dsr.33Rs$redD,Dsr.33Rs$oxD)

#calculate alpha total in each case
alphaT.33 <- Calc.alphaT(Rp.33,Dsr.33Rs$SO30,f)
alphaT.33.Delta = Calc.alphaT(Rp.33.Delta,Dsr.33Rs$SO30,f)

alphaT.33.a <- Calc.alphaT.a(Rp.33,Dsr.33Rs$SO30,f)
alphaT.33.Delta.a = Calc.alphaT.a(Rp.33.Delta,Dsr.33Rs$SO30,f)

#calculate alphas for oxidized and reduced moieties
alpha.red.33 = Calc.alphax(Dsr.33Rs$red, Dsr.33Rs$SO30, alphaT.33, f)
alpha.red.33.Delta = Calc.alphax(Dsr.33Rs$redD, Dsr.33Rs$SO30, alphaT.33.Delta, f)

alpha.red.33.a = Calc.alphax(Dsr.33Rs$red, Dsr.33Rs$SO30, alphaT.33.a, f)
alpha.red.33.Delta.a = Calc.alphax(Dsr.33Rs$redD, Dsr.33Rs$SO30, alphaT.33.Delta.a, f)

alpha.ox.33 = Calc.alphax(Dsr.33Rs$ox, Dsr.33Rs$SO30, alphaT.33, f)
alpha.ox.33.Delta = Calc.alphax(Dsr.33Rs$oxD, Dsr.33Rs$SO30, alphaT.33.Delta, f)

alpha.ox.33.a = Calc.alphax(Dsr.33Rs$ox, Dsr.33Rs$SO30, alphaT.33.a, f)
alpha.ox.33.Delta.a = Calc.alphax(Dsr.33Rs$oxD, Dsr.33Rs$SO30, alphaT.33.Delta.a, f)

############ 36 S ############ 
#Calculate total R of product - using measurements on 253 (default) and on Delta (labeled Delta)
Rp.36 = Calc.Rp(j,Dsr.36Rs$red,Dsr.36Rs$ox)
Rp.36.Delta = Calc.Rp(j,Dsr.36Rs$redD,Dsr.36Rs$oxD)

#calculate alpha total in each case
alphaT.36 = Calc.alphaT(Rp.36,Dsr.36Rs$SO30,f)
alphaT.36.Delta = Calc.alphaT(Rp.36.Delta,Dsr.36Rs$SO30,f)

alphaT.36.a = Calc.alphaT.a(Rp.36,Dsr.36Rs$SO30,f)
alphaT.36.Delta.a = Calc.alphaT.a(Rp.36.Delta,Dsr.36Rs$SO30,f)

#calculate alphas for oxidized and reduced moieties
alpha.red.36 = Calc.alphax(Dsr.36Rs$red, Dsr.36Rs$SO30, alphaT.36, f)
alpha.red.36.Delta = Calc.alphax(Dsr.36Rs$redD, Dsr.36Rs$SO30, alphaT.36.Delta, f)

alpha.red.36.a = Calc.alphax(Dsr.36Rs$red, Dsr.36Rs$SO30, alphaT.36.a, f)
alpha.red.36.Delta.a = Calc.alphax(Dsr.36Rs$redD, Dsr.36Rs$SO30, alphaT.36.Delta.a, f)

alpha.ox.36 = Calc.alphax(Dsr.36Rs$ox, Dsr.36Rs$SO30, alphaT.36, f)
alpha.ox.36.Delta = Calc.alphax(Dsr.36Rs$oxD, Dsr.36Rs$SO30, alphaT.36.Delta, f)

alpha.ox.36.a = Calc.alphax(Dsr.36Rs$ox, Dsr.36Rs$SO30, alphaT.36.a, f)
alpha.ox.36.Delta.a = Calc.alphax(Dsr.36Rs$oxD, Dsr.36Rs$SO30, alphaT.36.Delta.a, f)


#######   LAMBDA 33   #######   
##reduced
lambda.red.33 = log(alpha.red.33)/log(alpha.red.34)
lambda.red.33.Delta =  log(alpha.red.33.Delta)/log(alpha.red.34.Delta)

lambda.red.33.a = log(alpha.red.33.a)/log(alpha.red.34.a)
lambda.red.33.Delta.a =  log(alpha.red.33.Delta.a)/log(alpha.red.34.Delta.a)
##oxidized
lambda.ox.33 = log(alpha.ox.33)/log(alpha.ox.34)
lambda.ox.33.Delta =  log(alpha.ox.33.Delta)/log(alpha.ox.34.Delta)

lambda.ox.33.a = log(alpha.ox.33.a)/log(alpha.ox.34.a)
lambda.ox.33.Delta.a =  log(alpha.ox.33.Delta.a)/log(alpha.ox.34.Delta.a)

#######   LAMBDA 36   #######  
##reduced
lambda.red.36 = log(alpha.red.36)/log(alpha.red.34)
lambda.red.36.Delta =  log(alpha.red.36.Delta)/log(alpha.red.34.Delta)

lambda.red.36.a = log(alpha.red.36.a)/log(alpha.red.34.a)
lambda.red.36.Delta.a =  log(alpha.red.36.Delta.a)/log(alpha.red.34.Delta.a)

##oxidized
lambda.ox.36 = log(alpha.ox.36)/log(alpha.ox.34)
lambda.ox.36.Delta =  log(alpha.ox.36.Delta)/log(alpha.ox.34.Delta)

lambda.ox.36.a = log(alpha.ox.36.a)/log(alpha.ox.34.a)
lambda.ox.36.Delta.a =  log(alpha.ox.36.Delta.a)/log(alpha.ox.34.Delta.a)


#### Check mass balance   #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
SBalance = Dsr$SO30 - Dsr$SO3n - Total.thionate.S 
MassBalance.red = Calc.redbal(j,Dsr.34Rs$red)
MassBalance.ox = Calc.oxbal(j,Dsr.34Rs$ox)
MassBalance.reactant = Dsr$SO30*Dsr.34Rs$SO30 
MassBalance.product = Dsr.34Rs$SO3*Dsr.34Rs$SO3 + MassBalance.red + MassBalance.ox
MassBalance.quot = MassBalance.product/MassBalance.reactant  #calc mass balance by quotient
MassBalance.quot.missing.34 = MassBalance.reactant*(1-MassBalance.quot)  #equal to m*R of the missing pool

#also check mass balance via SO3 by difference method
SBalance.dif = Dsr$SO30 - SO3bydif - Total.thionate.S #must be 0 if SO3bydif
MassBalance.product.dif = SO3bydif*Dsr.34Rs$SO3 + MassBalance.red + MassBalance.ox
MassBalance.quot.dif = MassBalance.product.dif/MassBalance.reactant  #calc mass balance by quotient
MassBalance.quot.missing.34.dif = MassBalance.reactant*(1-MassBalance.quot.dif)  #equal to m*R of the missing pool


### determine how much mass is 'missing' at a range of Rs
trial.R.low = 0.95*cdt
trial.R.high = 1.05*cdt
trial.Rs = seq(trial.R.low,trial.R.high, length.out=25) #input the range of R's over which to iterate
trial.alphas = trial.Rs/cdt
trial.deltas = (trial.Rs/cdt - 1)*1000

num.Rs=length(trial.Rs)
Missing.mass = matrix(data=NA, numsamples,num.Rs)
counter = 0
for (tRval in trial.Rs){
  counter = counter + 1
  Missing.mass[,counter]=tRval*MassBalance.quot.missing.34
}

MM.complete = complete.cases(Missing.mass)
Missing.mass.good = Missing.mass[MM.complete,]
plot.new()        #set up plot for mass balance
plot.window(xlim=range(trial.deltas),ylim=range(Missing.mass.good)) #axes of plot
axis(1); axis(2); box()     #draw the axes & box
title(xlab='delta value', ylab='missing mass')    #label axes
                                                #now draw the lines
for (i in seq(1:length(Missing.mass.good[,1]))){    
  lines(trial.deltas,Missing.mass.good[i,])
}
    
###EXCLUDED: Mass balance by difference or as "delta"
#MassBalance.diff = MassBalance.reactant - MassBalance.product #calc mass balance by difference
#MassBalance.diff.R = MassBalance.diff/Dsr$SO30  #normalize the mass*R error to R by dividing by total mass
#MassBalance.diff.delta = R2Delta(abs(MassBalance.R),cdt) + 1000 #I don't know if this makes sense
#MassBalance.diff.percentage = MassBalance.diff/(Dsr$SO30*Dsr$d34SO30) #something like a percentage
#MassBalance.quot.delta = alpha2epsilon(MassBalance.frac) #convert quot to a delta-like number
####  #### #### END of MASS BALANCE #### #### #### #### #### #### #### #### 

####   mean(Reliable(lambda.red.33), na.rm=TRUE)

LineNo = seq(1:numsamples)
Results <- data.frame(LineNo, Dsr$ExNo, Dsr$Temp, Dsr$Hours, f, j, alphaT.34,
                        alphaT.33, alphaT.36, alpha.red.34, alpha.red.33, alpha.red.36, 
                        alpha.ox.34, alpha.ox.33, alpha.ox.36, lambda.red.33, lambda.red.36, 
                        lambda.ox.33, lambda.ox.36)
RResults <- Reliable(Results)
write.table(RResults,file="ReliableResults.csv",sep=",",row.names=F)

poslnf = -log(f)