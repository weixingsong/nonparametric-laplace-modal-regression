# From the intermediate results to the entries in Figure 3.

figure3 = read.csv("InterResult//Figure3//figure3.csv", skip = 1, header = FALSE)


pv=as.matrix(figure3[1])
fp100=as.matrix(figure3[2])
fp200=as.matrix(figure3[3])

postscript("results/figure3.eps", width = 12, height = 5,  
           horizontal = FALSE, paper = "special")

par(mfrow=c(1,2))
# Left panel of Figure 3
plot(pv,fp100,type="l",xlab="Nominal Significance Level (n=100)",ylab="Rejection Rate")
lines(pv,pv,lty=3)
# Right panel of Figure 3
plot(pv,fp200,type="l",xlab="Nominal Significance Level (n=200)",ylab="Rejection Rate")
lines(pv,pv,lty=3)

dev.off()
