# From the intermediate results to the entries in Figure 4.

figure4 = read.csv("InterResult//Figure4//figure4.csv", skip = 1, header = FALSE)


pv=as.matrix(figure4[1])
fp100=as.matrix(figure4[2])
fp200=as.matrix(figure4[3])

postscript("results/figure4.eps", width = 12, height = 5, 
           horizontal = FALSE, paper = "special")

par(mfrow=c(1,2))
# Left panel of Figure 4
plot(pv,fp100,type="l",xlab="Nominal Significance Level (n=100)",ylab="Rejection Rate")
lines(pv,pv,lty=3)
# Right panel of Figure 4
plot(pv,fp200,type="l",xlab="Nominal Significance Level (n=200)",ylab="Rejection Rate")
lines(pv,pv,lty=3)


dev.off()