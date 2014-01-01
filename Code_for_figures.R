

###Figure 2A
RedPlot<-data.frame(Reliable(alpha.red.34))
RedPlot<-cbind(RedPlot,Reliable(lambda.red.33))
RedPlot.complete = RedPlot[complete.cases(RedPlot),]
RedPlot.complete = cbind(RedPlot.complete, f.rank)
RedPlot.sorted = RedPlot.complete[order(RedPlot.complete[,3]),]
plot(RedPlot.sorted[,1],RedPlot.sorted[,2],ylim=c(0.490,0.525), xlim=c(0.982,0.992),xlab="alpha.red",ylab="lambda", pch=16, col=gray.colors(9,start=0.3,end=0.9))


###Figure2A
RedPlot<-data.frame(Reliable(alpha.red.34))
RedPlot<-cbind(RedPlot,Reliable(lambda.red.33))
RedPlot<-cbind(RedPlot,f)
RedPlot.complete = RedPlot[complete.cases(RedPlot),]
RedPlot.sorted<-RedPlot.complete[order(RedPlot.complete$f), ]
plot(RedPlot.sorted[,1],RedPlot.sorted[,2],ylim=c(0.490,0.525), xlim=c(0.982,0.992),xlab="alpha.ox",ylab="lambda",pch=16, col=gray.colors(9,start=0.3,end=0.9))



row.names(Dsr)
attr(Dsr,"row.names")

###Figure2B
OxPlot<-data.frame(Reliable(alpha.ox.34))
OxPlot<-cbind(OxPlot,Reliable(lambda.ox.33))
OxPlot<-cbind(OxPlot,f)
OxPlot.complete = OxPlot[complete.cases(OxPlot),]
OxPlot.sorted<-OxPlot.complete[order(OxPlot.complete$f), ]
plot(OxPlot.sorted[,1],OxPlot.sorted[,2],ylim=c(0.490,0.525), xlim=c(1.000,1.015),xlab="alpha.ox",ylab="lambda",pch=16, col=gray.colors(9,start=0.3,end=0.9))


###Figure 2C
###reverse plot with grayscale
alpha.range = seq(0.980, 1.1, length=99)
plot(alpha.range, Xs[1,], type="n", xlim=c(0.976,1.1),xlab="alpha.secondary",ylab="X")
lines(MixMod.complete[1,4:n],Xs[1,], col=gray.colors(9,start=0.3, end=0.9)[f.rank[1]])
lines(MixMod.complete[2,4:n],Xs[1,], col=gray.colors(9,start=0.3, end=0.9)[f.rank[2]])
lines(MixMod.complete[3,4:n],Xs[1,], col=gray.colors(9,start=0.3, end=0.9)[f.rank[3]])
lines(MixMod.complete[4,4:n],Xs[1,], col=gray.colors(9,start=0.3, end=0.9)[f.rank[4]])
lines(MixMod.complete[5,4:n],Xs[1,], col=gray.colors(9,start=0.3, end=0.9)[f.rank[5]])
lines(MixMod.complete[6,4:n],Xs[1,], col=gray.colors(9,start=0.3, end=0.9)[f.rank[6]])
#lines 7 and 8 excluded, because these are included in calc of alpha for f near 1
#lines(MixMod.complete[7,4:n],Xs[1,], col=gray.colors(9,start=0.3, end=0.9)[f.rank[7]])
#lines(MixMod.complete[8,4:n],Xs[1,], col=gray.colors(9,start=0.3, end=0.9)[f.rank[8]])
lines(MixMod.complete[9,4:n],Xs[1,], col=gray.colors(9,start=0.3, end=0.9)[f.rank[9]])
lines(MixMod.complete[10,4:n],Xs[1,], col=gray.colors(9,start=0.3, end=0.9)[f.rank[10]])
lines(MixMod.complete[11,4:n],Xs[1,], col=gray.colors(9,start=0.3, end=0.9)[f.rank[11]])

