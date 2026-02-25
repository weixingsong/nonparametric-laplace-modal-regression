# From the intermediate results to the entries in Figure B.16.

figureB16 = read.csv("InterResult//FigureB16//figureB16.csv", skip = 1, header = FALSE)



pv=as.matrix(figureB16[1])
fp100=as.matrix(figureB16[2])
fp200=as.matrix(figureB16[3])

postscript("results/figureB16.eps", width = 12, height = 5,  
           horizontal = FALSE, paper = "special")

par(mfrow=c(1,2))
# Left panel of Figure B.16
plot(pv,fp100,type="l",xlab="Nominal Significance Level (n=200)",ylab="Rejection Rate")
lines(pv,pv,lty=3)
# Right panel of Figure B.16
plot(pv,fp200,type="l",xlab="Nominal Significance Level (n=400)",ylab="Rejection Rate")
lines(pv,pv,lty=3)

dev.off()
