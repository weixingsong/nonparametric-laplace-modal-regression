# Simulation study for the paper "Parametric Modal Regression with Error Contaminated Covariates"

figure8=function(nsize=c(100,200,300),B=100,total=100,M=100)
 {
   path=getwd()
   source(paste0(path,"\\simulation 2\\Figure8 MCCM1.R"))
   source(paste0(path,"\\simulation 2\\Figure8 MCCM2.R"))
   source(paste0(path,"\\simulation 2\\Figure8 DC.R"))
   source(paste0(path,"\\simulation 2\\Figure8 CK.R"))

   model <- matrix(NA, nrow = length(nsize), ncol = 5)
   rownames(model) <- paste0("n=", nsize)
   colnames(model) <- c("figure8M1", "figure8M2", "figure8DC", "figure8CK1", "figure8CK2")

   for(i in seq_along(nsize)) 
    {
      n <- nsize[i]
  
      model[i, 1] <- figure8M1(n = n,B,1,1.3,total,M)
  
      model[i, 2] <- figure8M2(n = n,B,1,1.3,total,M)
  
      model[i, 3] <- figure8DC(n = 2*n,B,1,1.3,total,M)
  
      resultCK <- figure8CK(n = n,B,1,1.3,total)
      model[i, 4] <- resultCK[1]  
      model[i, 5] <- resultCK[2]  
    }

   plot(nsize,model[,1],type="l",ylim=c(min(model),max(model)),
     xlab="n (Model 2.3)",ylab="Rejection Rate",lty=1,lwd=2,xaxt="n")
   axis(side=1,at=c(100,200,300))
   lines(nsize,model[,2],lty=2,lwd=2)
   lines(nsize,model[,3],lty=3,lwd=2)
   lines(nsize,model[,4],lty=4,lwd=2)
   lines(nsize,model[,5],lty=5,lwd=2)

   legend("bottomright",                     # position
       legend = c("MCCM1","MCCM2","DC","CvM","KS"), # legend text
       lwd = 2,                        # line width
       lty = c(1,2,3,4,5),
       y.intersp = 1.2,
       bty = "n") 

   cat("Simulation Results in Figure 8","\n")
   cat("Power of MCCM1,MCCM2,DC,Cvm and KS in Model 2.3","\n\n") 
 }

btime=Sys.time()
postscript("results/figure8.eps", width = 12, height = 8,  
           horizontal = FALSE, paper = "special")

figure8(nsize = c(100, 200, 300), B = 100, total = 100, M = 100)

dev.off()
Sys.time()-btime

