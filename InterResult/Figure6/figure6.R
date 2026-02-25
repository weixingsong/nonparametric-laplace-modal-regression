# From the intermediate results to the entries in Figure 6.

figure6 = read.csv("InterResult//Figure6//figure6.csv", skip = 1, header = FALSE)

pv=as.matrix(figure6[1])
cvm100=as.matrix(figure6[2])
ks100=as.matrix(figure6[3])
cvm200=as.matrix(figure6[4])
ks200=as.matrix(figure6[5])

ymaxi=max(figure6)
ymini=min(figure6)

postscript("results/figure6.eps", width = 12, height = 5,  
           horizontal = FALSE, paper = "special")

par(mfrow=c(1,2))
# Left panel of Figure 6
plot(pv,cvm100,type="l",xlim=c(0.01,0.5),ylim=c(ymini,ymaxi),xlab="Nominal Significance Level (n=100)",ylab="Rejection Rate")
lines(pv,ks100,lty=2)
lines(pv,pv,lty=3)
# Right panel of Figure 6
plot(pv,cvm200,type="l",xlim=c(0.01,0.5),ylim=c(ymini,ymaxi),xlab="Nominal Significance Level (n=200)",ylab="Rejection Rate")
lines(pv,ks200,lty=2)
lines(pv,pv,lty=3)

dev.off()

