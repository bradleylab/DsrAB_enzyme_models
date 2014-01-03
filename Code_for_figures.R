
###Figure2A
RedPlot<-data.frame(Reliable(alpha.red.34))
RedPlot<-cbind(RedPlot,Reliable(lambda.red.33))
RedPlot<-cbind(RedPlot,f)
RedPlot.complete = RedPlot[complete.cases(RedPlot),]
plot(RedPlot.complete[,1],RedPlot.complete[,2],ylim=c(0.490,0.525), xlim=c(0.982,0.992),xlab="alpha.red",ylab="lambda",pch=16, col=rgb(1,RedPlot.complete$f,RedPlot.complete$f))

###Figure2B
OxPlot<-data.frame(Reliable(alpha.ox.34))
OxPlot<-cbind(OxPlot,Reliable(lambda.ox.33))
OxPlot<-cbind(OxPlot,f)
OxPlot.complete = OxPlot[complete.cases(OxPlot),]
plot(OxPlot.complete[,1],OxPlot.sorted[,2],ylim=c(0.490,0.525), xlim=c(1.000,1.015),xlab="alpha.ox",ylab="lambda",pch=16, col=rgb(1,OxPlot.complete$f,OxPlot.complete$f))


###Figure 2C
###reverse plot with grayscale
alpha.range = seq(0.980, 1.1, length=99)
plot(alpha.range, Xs[1,], type="n", xlim=c(0.976,1.1),xlab="alpha.secondary",ylab="X")
lines(MixMod.complete[1,4:n],Xs[1,], col=rgb(1,f.MM.complete[1],f.MM.complete[1]))
lines(MixMod.complete[2,4:n],Xs[1,], col=rgb(1,f.MM.complete[2],f.MM.complete[2]))
lines(MixMod.complete[3,4:n],Xs[1,], col=rgb(1,f.MM.complete[3],f.MM.complete[3]))
lines(MixMod.complete[4,4:n],Xs[1,], col=rgb(1,f.MM.complete[4],f.MM.complete[4]))
lines(MixMod.complete[5,4:n],Xs[1,], col=rgb(1,f.MM.complete[5],f.MM.complete[5]))
lines(MixMod.complete[6,4:n],Xs[1,], col=rgb(1,f.MM.complete[6],f.MM.complete[6]))
#lines 7 and 8 excluded, because these are included in calc of alpha for f near 1
#lines(MixMod.complete[7,4:n],Xs[1,], col=rgb(1,f.MM.complete[7],f.MM.complete[7]))
#lines(MixMod.complete[8,4:n],Xs[1,], col=rgb(1,f.MM.complete[8],f.MM.complete[8]))
lines(MixMod.complete[9,4:n],Xs[1,], col=rgb(1,f.MM.complete[9],f.MM.complete[9]))


