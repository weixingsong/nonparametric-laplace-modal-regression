# Simulation study for the paper "Nonparametric Modal Regression with Laplace Measurement Error"
#  Dietary data analysis 
# Estimate modal curve.
 
 rm(list=ls())
 btime = Sys.time()



# Load packages
 library(MASS)
 library(KernSmooth)
 library(readxl)  # 载入readxl包

  source("bandsel_huang.R")
  source("meanshift_huang.R")
  source("optimcor.R")
  source("CLy.bandwidth.R") 
  source("bandwidth_selection.R")

 set.seed(12345)

 gridlen=200
 reptime=500

 #original data

 data<- read.csv("wishreg.csv", header=T)

 y=data$ffq
 z=data[,2:7]
 Yo=log(y)
 z=log(z)
 Zo=apply(z,1,mean)
 n=length(Yo)
 su=sqrt(0.022) 
 su2=su^2

 xseq=seq(min(Zo),max(Zo),length=gridlen)
 
 
 m0cor = m0naive = m0huang = matrix(0,nrow=reptime,ncol=gridlen)
 m1cor = m1naive = matrix(0,nrow=reptime,ncol=gridlen)

 for(j in seq(reptime))
    {
   
           ##resample
         bootstrap_sample <- function(Za, Ya) 
            {
                 idx <- sample(n, replace = TRUE)  
                 list(Za = Za[idx], Ya = Ya[idx])  
            }

         resampled_data <- bootstrap_sample(Zo, Yo)
         Z=resampled_data$Za
         Y=resampled_data$Ya

   
         # Define the objective function
   
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

               m0huang[j,k] = myopthuang
               m0cor[j,k] = myoptcor[1]
               m0naive[j,k] = myoptnaive[1]
               m1cor[j,k] = myoptcor[2]
               m1naive[j,k] = myoptnaive[2]
               k=k+1
           }
        cat(j,"\n")
    } 
   
    corhuang=(1/2)*(m0cor+m0huang)

    mod0cormean = apply(m0cor,2,mean)
    mod0naivemean = apply(m0naive,2,mean)
    mod0huangmean = apply(m0huang,2,mean)
    mod0corhuangmean = apply(corhuang,2,mean)


    mod0cortmean = apply(m0cor,2,mean,trim=0.05)
    mod0naivetmean = apply(m0naive,2,mean,trim=0.05)
    mod0huangtmean = apply(m0huang,2,mean,trim=0.05)
    mod0corhuangtmean = apply(corhuang,2,mean,trim=0.05)

    mod0cormedian = apply(m0cor,2,median)
    mod0naivemedian = apply(m0naive,2,median)
    mod0huangmedian = apply(m0huang,2,median)
    mod0corhuangmedian = apply(corhuang,2,median)

    misenaivemean = mean(rowMeans((sweep(m0naive, 2, mod0naivemean,"-"))^2))
    misecormean = mean(rowMeans((sweep(m0cor, 2, mod0cormean,"-"))^2))
    misehuangmean = mean(rowMeans((sweep(m0huang, 2, mod0huangmean))^2))
    misecorhuangmean = mean(rowMeans((sweep(corhuang, 2, mod0corhuangmean,"-"))^2))

    Var <- rbind(misenaivemean, misecormean,  misehuangmean, misecorhuangmean)
    rownames(Var) <- c("Naive", "Corrected", "Huang", "Corhuang")
    Var


tm=Sys.time()-btime  
save(file="Diet_boot.re=500.RData",list=ls())