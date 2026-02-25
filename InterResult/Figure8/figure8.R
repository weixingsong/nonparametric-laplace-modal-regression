# From the intermediate results to the entries in Figure 8.

figure8 =read.csv("InterResult//Figure8//figure8.csv", skip = 1, header = FALSE)

pfun=function(x){mean(x<=0.05)}


Pcvmks=as.vector(apply(figure8[10:15],2,mean))
Pmccldc=as.vector(apply(figure8[1:9],2,pfun))
Pvalue=c(Pmccldc,Pcvmks)


# Figure 8: Powers of MCCM1, MCCM2, DC, CvM and KS in Model 2.3

postscript("results/figure8.eps", width = 12, height = 5,  
           horizontal = FALSE, paper = "special")

par(mfrow=c(1,1))
model=matrix(Pvalue,nrow=5,byrow=T)
nsize=c(100,200,300)
plot(nsize,model[1,],type="l",ylim=c(min(model),max(model)),
     xlab="n (Model 2.3)",ylab="Rejection Rate",lty=1,lwd=2,xaxt="n")
axis(side=1,at=c(100,200,300))
lines(nsize,model[2,],lty=2,lwd=2)
lines(nsize,model[3,],lty=3,lwd=2)
lines(nsize,model[4,],lty=4,lwd=2)
lines(nsize,model[5,],lty=5,lwd=2)

legend("bottomright",                     # position
       legend = c("MCCM1","MCCM2","DC","CvM","KS"), # legend text
       lwd = 2,                        # line width
       lty = c(1,2,3,4,5),
       y.intersp = 1.2,
       bty = "n") 

dev.off()
