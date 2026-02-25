# Simulation study for the paper "Parametric Modal Regression with Error Contaminated Covariates"

modelB3=function(nsize=c(100,200,300),B=100,total=100,M=100)
 {
   path=getwd()
   source(paste0(path,"\\simulation B2\\FigureB18 Model B3.R"))
   source(paste0(path,"\\simulation B2\\Model B3M2.R"))
   source(paste0(path,"\\simulation B2\\Model B3DC.R"))
   source(paste0(path,"\\simulation B2\\Model B3CK.R"))

   modelB3 <- matrix(NA, nrow = length(nsize), ncol = 5)
   rownames(modelB3) <- paste0("n=", nsize)
   colnames(modelB3) <- c("MCCM1", "MCCM2", "DC", "CvM", "KS")

   for(i in seq_along(nsize)) 
    {
      n <- nsize[i]
  
      modelB3[i, 1] <- figureB18B3(n=n,B,1,total,M)
  
      modelB3[i, 2]  <- modelB3M2(n=n,B,1,total,M)
  
      modelB3[i, 3]<- modelB3DC(n = 2*n,B,1,total,M)
  
      resultCK <- modelB3CK(n = n,B,1,total)
      modelB3[i, 4] <- resultCK[1]  
      modelB3[i, 5] <- resultCK[2]  
    }

   plot(nsize,modelB3[,1],type="l",ylim=c(min(modelB3),max(modelB3)),
     xlab="n (Model B.3)",ylab="Rejection Rate",lty=1,lwd=2,xaxt="n")
   axis(side=1,at=c(100,200,300))
   lines(nsize,modelB3[,2],lty=2,lwd=2)
   lines(nsize,modelB3[,3],lty=3,lwd=2)
   lines(nsize,modelB3[,4],lty=4,lwd=2)
   lines(nsize,modelB3[,5],lty=5,lwd=2)

   legend("bottomright",                     # position
       legend = c("MCCM1","MCCM2","DC","CvM","KS"), # legend text
       lwd = 2,                        # line width
       lty = c(1,2,3,4,5),
       y.intersp = 1.2,
       bty = "n") 

   modelB3
   cat("Simulation Results in Model B.3","\n")
   cat("Rejection rate of MCCM1,MCCM2,DC,Cvm and KS in Model B.3","\n\n") 
 }

btime=Sys.time()
sink("results/modelB3_results.txt")
modelB3(nsize=c(100,200,300),B=100,total=100,M=100)
sink()
Sys.time()-btime

