# From the intermediate results to the entries in Figure B.15.

figureB15 = read.csv("InterResult//FigureB15//figureB15.csv", skip = 1, header = FALSE)


pv=as.matrix(figureB15[1])
fp100=as.matrix(figureB15[2])
fp200=as.matrix(figureB15[3])

postscript("results/figureB15.eps", width = 12, height = 5,  
           horizontal = FALSE, paper = "special")
par(mfrow=c(1,2))
# Left panel of Figure B.15
plot(pv,fp100,type="l",xlab="Nominal Significance Level (n=100)",ylab="Rejection Rate")
lines(pv,pv,lty=3)
# Right panel of Figure B.15
plot(pv,fp200,type="l",xlab="Nominal Significance Level (n=200)",ylab="Rejection Rate")
lines(pv,pv,lty=3)


dev.off()