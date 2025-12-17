# Simulation study for the paper "Nonparametric Modal Regression with Laplace Measurement Error"
#  Dietary data analysis 
# Estimate modal curve and var.

rm(list=ls())
ts=Sys.time()

# Load packages
 library(MASS)
 library(KernSmooth)
 library(readxl)  # 载入readxl包

  source("bandsel_huang.R")
  source("meanshift_huang.R")
  source("optimcor.R")
  source("CLy.bandwidth.R") 
  source("bandwidth_selection.R")

 data<- read.csv("wishreg.csv", header=T)

  y=data$ffq
  z=data[,2:7]
  Y=log(y)
  z=log(z)

  n=length(Y)   # sample size
  gridlen=200

 set.seed(123)
 
 #   # Define the objective function
   su=sqrt(0.022)
   Z=apply(z,1,mean)
   xseq=seq(min(Z),max(Z),length=gridlen)
   su2=su^2 

   h0 = n^(-1/7)    
   bandn = band.sel(Y,Z,su,h0,method="Naive")
   h1n = bandn[1]
   h2n = bandn[2]
   print(c(h1n,h2n))
   
   bandc = CLy_band_width(Z,Y,su,n) 
   h1c =  bandc[1]
   h2c =  bandc[2]
   print(c(h1c,h2c))
   
   bandh = bandsel_huang(Y,Z,su)
   h1h = bandh[1]
   h2h = bandh[2]
   print(c(h1h,h2h))
   m0cor = m0naive = m0huang = rep(0,gridlen)
   m1cor = m1naive = rep(0,gridlen)
      
   k=1
   for(x0 in xseq)
     {
     
         # The Proposed Local Linear Method

         y1_in=7.2
         y2_in=0.2
         myoptcor = optimcor(Z,Y,y1_in,y2_in,x0,h1c,h2c)
         
      
        # The Naive Method
      
        objfun_naive=function(bt)
          {
               mean(dnorm(Z,x0,h1n)*dnorm(Y,bt[1]+bt[2]*(Z-x0),h2n))
           }
            
        myoptnaive=optim(c(7.2,0.2),objfun_naive, control=list(fnscale=-1))$par
      
    
      # Huang's Local Linear Method
      
      # Mean-Shift Estimate   
      
        myopthuang = meanshift(Z,Y,7.1,x0,h1h,h2h)[1]
      
        m0huang[k]=myopthuang
        m0cor[k]=myoptcor[1]
        m0naive[k]=myoptnaive[1]
        k=k+1
        print(k)
     }

     corhuang=(1/2)*(m0cor+m0huang)

     postscript("output.eps", 
           width = 5,   
           height = 5,  
           horizontal = FALSE,  
           paper = "special",   
           onefile = TRUE      
      )

 
   par(mfrow = c(1, 1), mar = c(6, 4, 2, 2) + 0.1, oma = c(1, 0, 0, 0), bty = "l")

    library(latex2exp)    


   x_min <- quantile(xseq,0.1)
   x_max <- quantile(xseq,0.99)
   idx <- xseq >= x_min & xseq <= x_max
   plot(Z,Y,xlab=TeX(r"(log long-term usual intake)"),ylab="log food frequency intake",col = "gray")
   lines(xseq[idx], m0naive[idx], type = "l", lwd = 2, col = "green")
   lines(xseq[idx], m0huang[idx], lwd = 2, col = "blue")
   lines(xseq[idx], m0cor[idx], lwd = 2, col = "red")
   lines(xseq[idx], corhuang[idx], lwd = 2, col = "cyan")

   legend(x=6.3,y=8.9,
       legend = c("LocPoly", "LocLinear", "LocPL", "Naive"),
       col = c("red", "blue", "cyan", "green"),
       lwd = c(2, 2, 2, 2),
       lty = c(1, 1, 1, 1),
       horiz = FALSE,
       xpd = NA,
       cex = 0.8,  
       y.intersp = 1.0,   
       text.width = NULL,
       bty = "n"
       )

dev.off()

te=Sys.time()
time=te-ts

