# Simulation study for the paper "Parametric Modal Regression with Error Contaminated Covariates"

figureB18l=function(nsize=c(100,200,300),B=100,total=100,M=100)
 {
   path=getwd()
   source(paste0(path,"\\simulation B2\\FigureB18 Model B1.R"))
   source(paste0(path,"\\simulation B2\\FigureB18 Model B2.R"))
   source(paste0(path,"\\simulation B2\\FigureB18 Model B3.R"))
   source(paste0(path,"\\simulation B2\\FigureB18 Model B4.R"))

   modelB <- matrix(NA, nrow = length(nsize), ncol = 4)
   rownames(modelB) <- paste0("n=", nsize)
   colnames(modelB) <- c("figure18B1", "figure18B2", "figure18B3", "figure18B4")

   for(i in seq_along(nsize)) 
    {
      n <- nsize[i]
  
      modelB[i, 1] <- figureB18B1(n,B,1,total,M)
  
      modelB[i, 2] <- figureB18B2(n,B,1,total,M)
  
      modelB[i, 3] <- figureB18B3(n,B,1,total,M)
  
      modelB[i, 4] <- figureB18B4(n,B,1,total,M) 
    }

   plot(nsize,modelB[,1],type="l",ylim=c(min(modelB)-0.1,max(modelB)),
     xlab="n",ylab="Rejection Rate",lty=1,lwd=2,xaxt="n")
   axis(side=1,at=c(100,200,300))
   lines(nsize,modelB[,2],lty=2,lwd=2)
   lines(nsize,modelB[,3],lty=3,lwd=2)
   lines(nsize,modelB[,4],lty=4,lwd=2)

   legend("bottomright",                     # position
       legend = c("Model B.1","Model B.2","Model B.3","Model B.4"), # legend text
       lwd = 2,                        # line width
       lty = c(1,2,3,4),
       y.intersp = 1.2,
       bty = "n")   

   cat("Simulation Results in Figure B.18(left)","\n")
   cat("MCCM1 test","\n\n") 
 }

btime=Sys.time()
postscript("results/figureB18l.eps", width = 12, height = 8,  
           horizontal = FALSE, paper = "special")

figureB18l(nsize=c(100,200,300),B=100,total=100,M=100)

dev.off()

Sys.time()-btime

