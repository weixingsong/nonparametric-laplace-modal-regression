# From the intermediate results to the entries in Figure B.18.

figureB18 = read.csv("InterResult//FigureB18//figureB18.csv", skip = 1, header = FALSE)

pfun=function(x){mean(x<=0.05)}

Pvaluel=as.vector(apply(figureB18[1:12],2,pfun))

Pmccldc=as.vector(apply(figureB18[10:18],2,pfun))
Pcvmks=as.vector(apply(figureB18[19:24],2,mean))
Pvaluer=c(Pmccldc,Pcvmks)

postscript("results/figureB18.eps", width = 12, height = 5,  
           horizontal = FALSE, paper = "special")

# Figure B.18(left panel): Powers of MCCM1 in Model B.1, B.2, B.3 and B.4

par(mfrow=c(1,2))

nsize=c(100,200,300)
modell=matrix(Pvaluel,nrow=4,byrow=T)

plot(nsize,modell[1,],type="l",ylim=c(min(modell)-0.1,max(modell)),
     xlab="n",ylab="Rejection Rate",lty=1,lwd=2,xaxt="n")
axis(side=1,at=c(100,200,300))
lines(nsize,modell[2,],lty=2,lwd=2)
lines(nsize,modell[3,],lty=3,lwd=2)
lines(nsize,modell[4,],lty=4,lwd=2)
legend("bottomright",                     # position
       legend = c("Model B.1","Model B.2","Model B.3","Model B.4"), # legend text
       lwd = 2,                        # line width
       lty = c(1,2,3,4),
       y.intersp = 1.2,
       bty = "n")   

# Figure B.18(right panel): Powers of MCCM1, MCCM2, DC, CvM and KS in Model B.4

modelr=matrix(Pvaluer,nrow=5,byrow=T)
nsize=c(100,200,300)
plot(nsize,modelr[1,],type="l",ylim=c(min(modelr),max(modelr)),
     xlab="n",ylab="Rejection Rate",lty=1,lwd=2,xaxt="n")
axis(side=1,at=c(100,200,300))
lines(nsize,modelr[2,],lty=2,lwd=2)
lines(nsize,modelr[4,],lty=3,lwd=2)
lines(nsize,modelr[5,],lty=4,lwd=2)
lines(nsize,modelr[3,],lty=5,lwd=2)

legend("bottomright",                     # position
       legend = c("MCCM1","MCCM2","CvM","KS","DC"), # legend text
       lwd = 2,                        # line width
       lty = c(1,2,3,4,5),
       y.intersp = 1.2,
       bty = "n") 


dev.off()
