############tests to run

aT = Reliable(alphaT.34)
aT.a = Reliable(alphaT.34.a)
qqnorm(alphaT.34)
qqlines(alphaT.34)
shapiro.test(alphaT.34)

red.253 = Reliable(alpha.red.34.a)
red.Delta = Reliable(alpha.red.34.Delta.a)
mean(red.253, na.rm=TRUE)

ks.test(red.253, red.Delta)

xlab.name = "-log(f)"
#ylab.name=expression(paste(delta^ 34,"S sulfonate"))
ylab.name = "delta 34S sulfonate"
  
  
poslnf = -log(f)
oxd34 = Reliable(Dsr.34deltas$d34ox)
redd34 = Reliable(Dsr.34deltas$d34ox)
plot(poslnf, oxd34, xlab = xlab.name, ylab=ylab.name)
oxmod <- lm(oxd34 ~ poslnf)
redmod <- lm(redd34 ~ poslnf)
