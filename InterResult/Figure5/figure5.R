# From the intermediate results to the entries in Figure 5.

figure5 = read.csv("InterResult//Figure5//figure5.csv", skip = 1, header = FALSE)


pv=as.matrix(figure5[1])
fp100=as.matrix(figure5[2])
fp200=as.matrix(figure5[3])

postscript("results/figure5.eps", width = 12, height = 5,  
           horizontal = FALSE, paper = "special")
par(mfrow=c(1,2))
# Left panel of Figure 5
plot(pv,fp100,type="l",xlab="Nominal Significance Level (n=200)",ylab="Rejection Rate")
lines(pv,pv,lty=3)
# Right panel of Figure 5
plot(pv,fp200,type="l",xlab="Nominal Significance Level (n=400)",ylab="Rejection Rate")
lines(pv,pv,lty=3)

dev.off()
