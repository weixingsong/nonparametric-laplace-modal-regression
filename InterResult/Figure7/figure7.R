# From the intermediate results to the entries in Figure 7.

figure7 = read.csv("InterResult//Figure7//figure7.csv", skip = 1, header = FALSE)

pfun=function(x){mean(x<=0.05)}
Pvalue=as.vector(apply(figure7,2,pfun))

# Figure 7: Powers of MCCM1 in Model 2.1, 2.2, and 2.3

postscript("results/figure7.eps", width = 12, height = 5,  
           horizontal = FALSE, paper = "special")

par(mfrow=c(1,3))

nsize=c(100,200,300)

plot(nsize,Pvalue[1:3],type="l",ylim=c(min(Pvalue[1:9]),max(Pvalue[1:9])),
     xlab="n (Model 2.1)",ylab="Rejection Rate",lty=3,lwd=2,xaxt="n")
axis(side=1,at=c(100,200,300))
lines(nsize,Pvalue[4:6],lty=2,lwd=2)
lines(nsize,Pvalue[7:9],lty=1,lwd=2)
legend("bottomright",                     # position
       legend = c(expression(beta[4] == -0.25), expression(beta[4]==-0.08),expression(beta[4]==-0.05)), # legend text
       lwd = 2,                        # line width
       lty = c(1,2,3),
       y.intersp = 1.2,
       bty = "n")   

plot(nsize,Pvalue[10:12],type="l",ylim=c(min(Pvalue[10:21]),max(Pvalue[10:21])),
     xlab="n (Model 2.2)",ylab="Rejection Rate",lty=4,lwd=2,xaxt="n")
axis(side=1,at=c(100,200,300))
lines(nsize,Pvalue[13:15],lty=3,lwd=2)
lines(nsize,Pvalue[16:18],lty=2,lwd=2)
lines(nsize,Pvalue[19:21],lty=1,lwd=2)
legend("bottomright",                     # position
       legend = c(expression(beta[4] == 1), expression(beta[4]==0.27),expression(beta[4]==0.25),expression(beta[4]==0.2)), # legend text
       lwd = 2,                        # line width
       lty = c(1,2,3,4),
       y.intersp = 1.2,
       bty = "n") 

plot(nsize,Pvalue[22:24],type="l",ylim=c(min(Pvalue[22:33]),max(Pvalue[22:33])),
     xlab="n (Model 2.3)",ylab="Rejection Rate",lty=4,lwd=2,xaxt="n")
axis(side=1,at=c(100,200,300))
lines(nsize,Pvalue[25:27],lty=3,lwd=2)
lines(nsize,Pvalue[28:30],lty=2,lwd=2)
lines(nsize,Pvalue[31:33],lty=1,lwd=2)
legend("bottomright",                     # position
       legend = c(expression(beta[4] == 2), expression(beta[4]==1.4),expression(beta[4]==1.3),expression(beta[4]==1.2)), # legend text
       lwd = 2,                        # line width
       lty = c(1,2,3,4),
       y.intersp = 1.2,
       bty = "n") 

dev.off()