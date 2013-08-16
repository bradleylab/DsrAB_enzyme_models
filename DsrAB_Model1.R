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

#function to calculate thiosulfate fraction J
CalcJ <-function(ThiosulfS, TrithioS)
{
  J=(ThiosulfS/(ThiosulfS+TrithioS)) 
  return(J)
}

#function to calculate isotope ratio of the total product (Rp)
CalcRp <-function(J, Rred, Rox)
{
  Rp=((1/2)*J+(1/3)*(1-J))*Rred + ((1/2)*J+(2/3)*(1-J))*Rox
  return(Rp)
}

#function to calculate alpha-total (alphaT)
Calc.alphaT <-function(Rp, Rao, f)
{
  alphaT = log(((Rp/Rao)*(f-1))+1)/(log(f))
  return(alphaT)
}

#function to calculate alpha-x (alpha-x where x is reduced or oxidized)
Calc.alphax <-function(Rx, Rao, alphaT, f)
{
  alphax = (Rx/Rao)*(alphaT/((f^alphaT)-1))*(f-1)
  return(alphax)
}

###### Start of program  ###### 
Dsr<-read.table("/Users/abradley/Documents/Rdata/WilDataFrame.txt",header=TRUE)
attach(Dsr)
names(Dsr)
cdt = .045005  		#CDT ratio for 34/32. 
cdt3x = 0.007379		#CDT ratio for 3x/32. 

#calculate j,f
j = CalcJ(Dsr$Thiosulf,Dsr$Trithio)
f = Dsr$SO3n/Dsr$SO30

#Pull out delta values & convert to Rs
Dsr.34deltas <- Dsr[, c(10:16)]
Dsr.33deltas <- Dsr[,c(17:20)]
Dsr.36deltas <- Dsr[,c(21:24)]
Dsr.34Rs <- Delta2R(Dsr.34deltas,cdt)
colnames(Dsr.34Rs) <- c("SO30", "SO3D","SO3","oxD","ox","redD","red")
Dsr.33Rs <- Delta2R(Dsr.33detlas,cdt3x)
colnames(Dsr.33Rs) <- c("SO30", "SO3", "ox", "red")
Dsr.36Rs <- Delta2R(Dsr.36deltas,cdt3x)
colnames(Dsr.36Rs) <- c("SO30", "SO3", "ox", "red")

############ 34 S ############ 
#Calculate total R of product - using measurements on 253 (default) and on Delta (labeled Delta)
Rp.34 = CalcRp(j,Dsr.34Rs$red,Dsr.34Rs$ox)
Rp.34.Delta = CalcRp(j,Dsr.34Rs$redD,Dsr.34Rs$oxD)

#calculate alpha total in each case
alphaT.34 = Calc.alphaT(Rp.34, Dsr.34Rs$SO30, f)
alphaT.34.Delta = Calc.alphaT(Rp.34.Delta,Dsr.34Rs$SO30,f)

#calculate alphas for oxidized and reduced moieties
alpha.red.34 = Calc.alphax(Dsr.34Rs$red, Dsr.34Rs$SO30, alphaT.34, f)
alpha.red.34.Delta = Calc.alphax(Dsr.34Rs$redD, Dsr.34Rs$SO30, alphaT.34.Delta, f)

alpha.ox.34 = Calc.alphax(Dsr.34Rs$ox, Dsr.34Rs$SO30, alphaT.34, f)
alpha.ox.34.Delta = Calc.alphax(Dsr.34Rs$oxD, Dsr.34Rs$SO30, alphaT.34.Delta, f)

############ 33 S ############ 
#Calculate total R of product - using measurements on 253 (default) and on Delta (labeled Delta)
Rp.33 = CalcRp(j,Dsr.33Rs$red,Dsr.33Rs$ox)
Rp.33.Delta = CalcRp(j,Dsr.33Rs$redD,Dsr.33Rs$oxD)

#calculate alpha total in each case
alphaT.33 = Calc.alphaT(Rp.33,Dsr.33Rs$SO30,f)
alphaT.33.Delta = Calc.alphaT(Rp.33.Delta,Dsr.33Rs$SO30,f)

#calculate alphas for oxidized and reduced moieties
alpha.red.33 = Calc.alphax(Dsr.33Rs$red, Dsr.33Rs$SO30, alphaT.33, f)
alpha.red.33.Delta = Calc.alphax(Dsr.33Rs$redD, Dsr.33Rs$SO30, alphaT.33.Delta, f)

alpha.ox.33 = Calc.alphax(Dsr.33Rs$ox, Dsr.33Rs$SO30, alphaT.33, f)
alpha.ox.33.Delta = Calc.alphax(Dsr.33Rs$oxD, Dsr.33Rs$SO30, alphaT.33.Delta, f)

############ 36 S ############ 
#Calculate total R of product - using measurements on 253 (default) and on Delta (labeled Delta)
Rp.36 = CalcRp(j,Dsr.36Rs$red,Dsr.36Rs$ox)
Rp.36.Delta = CalcRp(j,Dsr.36Rs$redD,Dsr.36Rs$oxD)

#calculate alpha total in each case
alphaT.36 = Calc.alphaT(Rp.36,Dsr.36Rs$SO30,f)
alphaT.36.Delta = Calc.alphaT(Rp.34.Delta,Dsr.34Rs$SO30,f)

#calculate alphas for oxidized and reduced moieties
alpha.red.34 = Calc.alphax(Dsr.34Rs$red, Dsr.34Rs$SO30, alphaT.34, f)
alpha.red.34.Delta = Calc.alphax(Dsr.34Rs$redD, Dsr.34Rs$SO30, alphaT.34.Delta, f)

alpha.ox.34 = Calc.alphax(Dsr.34Rs$ox, Dsr.34Rs$SO30, alphaT.34, f)
alpha.ox.34.Delta = Calc.alphax(Dsr.34Rs$oxD, Dsr.34Rs$SO30, alphaT.34.Delta, f)


############################# 


