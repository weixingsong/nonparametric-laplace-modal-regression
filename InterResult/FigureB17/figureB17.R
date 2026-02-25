# From the intermediate results to the entries in Figure B.17.

figureB17 = read.csv("InterResult//FigureB17//figureB17.csv", skip = 1, header = FALSE)

pv=as.matrix(figureB17[1])
cvm100=as.matrix(figureB17[2])
ks100=as.matrix(figureB17[3])
cvm200=as.matrix(figureB17[4])
ks200=as.matrix(figureB17[5])

ymaxi=max(figureB17)
ymini=min(figureB17)

postscript("results/figureB17.eps", width = 12, height = 5,  
           horizontal = FALSE, paper = "special")

par(mfrow=c(1,2))
# Left panel of Figure B.17
plot(pv,cvm100,type="l",xlim=c(0.01,0.5),ylim=c(ymini,ymaxi),xlab="Nominal Significance Level (n=100)",ylab="Rejection Rate")
lines(pv,ks100,lty=2)
lines(pv,pv,lty=3)
# Right panel of Figure B.17
plot(pv,cvm200,type="l",xlim=c(0.01,0.5),ylim=c(ymini,ymaxi),xlab="Nominal Significance Level (n=200)",ylab="Rejection Rate")
lines(pv,ks200,lty=2)
lines(pv,pv,lty=3)

dev.off()
