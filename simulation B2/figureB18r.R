# Simulation study for the paper "Parametric Modal Regression with Error Contaminated Covariates"

figureB18r=function(nsize=c(100,200,300),B=100,total=100,M=100)
 {
   path=getwd()
   source(paste0(path,"\\simulation B2\\FigureB18 Model B4.R"))
   source(paste0(path,"\\simulation B2\\FigureB18 MCCM2.R"))
   source(paste0(path,"\\simulation B2\\FigureB18 DC.R"))
   source(paste0(path,"\\simulation B2\\FigureB18 CK.R"))

   modelB4 <- matrix(NA, nrow = length(nsize), ncol = 5)
   rownames(modelB4) <- paste0("n=", nsize)
   colnames(modelB4) <- c("figure8M1", "figure8M2", "figure8DC", "figure8CK1", "figure8CK2")

   for(i in seq_along(nsize)) 
    {
      n <- nsize[i]
  
      modelB4[i, 1] <- figureB18B4(n=n,B,1,total,M)
  
      modelB4[i, 2] <- figureB18B4M2(n=n,B,1,total,M)
  
      modelB4[i, 3] <- figureB18B4DC(n=2*n,B,1,total,M)
  
      resultCK <- figureB18B4CK(n=n,B,1,total)
      modelB4[i, 4] <- resultCK[1]  
      modelB4[i, 5] <- resultCK[2]  
    }

   plot(nsize,modelB4[,1],type="l",ylim=c(min(modelB4),max(modelB4)),
     xlab="n",ylab="Rejection Rate",lty=1,lwd=2,xaxt="n")
   axis(side=1,at=c(100,200,300))
   lines(nsize,modelB4[,2],lty=2,lwd=2)
   lines(nsize,modelB4[,3],lty=3,lwd=2)
   lines(nsize,modelB4[,4],lty=4,lwd=2)
   lines(nsize,modelB4[,5],lty=5,lwd=2)

   legend("bottomright",                     # position
       legend = c("MCCM1","MCCM2","DC","CvM","KS"), # legend text
       lwd = 2,                        # line width
       lty = c(1,2,3,4,5),
       y.intersp = 1.2,
       bty = "n") 

   cat("Simulation Results in Figure B.18(right)","\n")
   cat("MCCM1,MCCM2,DC,Cvm and KS test in Model B.4","\n\n") 
 }

btime=Sys.time()
postscript("results/figureB18r.eps", width = 12, height = 8,  
           horizontal = FALSE, paper = "special")

figureB18r(nsize=c(100,200,300),B=100,total=100,M=100)

dev.off()
Sys.time()-btime

