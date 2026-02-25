# From the intermediate results to the entries in Figure 2.

library(latex2exp)

figure2 = read.csv("InterResult//Figure2//figure2.csv", skip = 3, header = FALSE)

postscript("results/figure2.eps", width = 12, height = 5,  
           horizontal = FALSE, paper = "special")

par(mfrow = c(1, 2))


# Left panel of Figure 2
result=as.matrix(figure2[1:4])
Nresult=as.matrix(figure2[5:8])


resultNM=cbind(result[,1],Nresult[,1],result[,2],Nresult[,2],result[,3],Nresult[,3],result[,4],Nresult[,4])
boxplot.matrix(resultNM,xlab = TeX(r"($\textbf{X}_1$ and $\textbf{X}_2$ are dependent)"),col=c("gray61","gray94","gray61","gray94","gray61","gray94","gray61","gray94"),  at=c(1,2,4,5,7,8,10,11),xaxt="n")
xtick = c(1.5,4.5, 7.5, 10.5)
names=c(TeX(r"( $\beta_{0}$ )"),TeX(r"( $\beta_{1}$ )"),TeX(r"( $\beta_{2}$ )"),TeX(r"( $\beta_{3}$ )"))
axis(side=1, at=xtick, labels = names)
abline(h=-0.25,lty=3)

# Right panel of Figure 2
result=as.matrix(figure2[9:12])
Nresult=as.matrix(figure2[13:16])


resultNM=cbind(result[,1],Nresult[,1],result[,2],Nresult[,2],result[,3],Nresult[,3],result[,4],Nresult[,4])
boxplot.matrix(resultNM,xlab = TeX(r"($\textbf{X}_1$ and $\textbf{X}_2$ are independent)"),col=c("gray61","gray94","gray61","gray94","gray61","gray94","gray61","gray94"),  at=c(1,2,4,5,7,8,10,11),xaxt="n")
xtick = c(1.5,4.5, 7.5, 10.5)
names=c(TeX(r"( $\beta_{0}$ )"),TeX(r"( $\beta_{1}$ )"),TeX(r"( $\beta_{2}$ )"),TeX(r"( $\beta_{3}$ )"))
axis(side=1, at=xtick, labels = names)
abline(h=-0.25,lty=3)

dev.off()