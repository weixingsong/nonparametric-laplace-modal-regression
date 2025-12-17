# Simulation study for the paper
#  Nonparametric Modal Regression with Laplace Measurement Error
# 
# Simulation 1

 #install.packages("KernSmooth")
 rm(list=ls())
 btime = Sys.time()
 library(KernSmooth)
 library(ks)

 # Empirical Optimal Bandwidth Selection Based on Theory 

  source("bandsel_huang.R")
  source("meanshift_huang.R")
  source("optimcor.R")
  source("hss.option4.R")
  source("naive.option4.R")
 

 # Data

  n=200
  gridlen=100
  reptime=500
   
  xseq=seq(0,1,length=gridlen)
 
  mod0true = 2*sin(pi*xseq)+1+2*xseq
  mod1true = 2*pi*cos(pi*xseq)+2
 
 
  mod0cor = mod0naive = mod0huang = rep(0,gridlen)
  mod1cor = mod1naive = rep(0,gridlen)
 
  m0cor = m0naive = m0huang = matrix(0,nrow=reptime,ncol=gridlen)
  m1cor = m1naive = matrix(0,nrow=reptime,ncol=gridlen)

  set.seed(12345)
  for(j in seq(reptime))
   {
     #data1  Yao
     utemp = runif(n,0,1)
     e = (utemp<=0.5)*rnorm(n,-1,2.5)+(utemp>0.5)*rnorm(n,1,0.5)
     X = runif(n,0,1)
     su = sqrt((1/0.75-1)*(1/12))
     #su = sqrt((1/0.90-1)*(1/12))
     su2 = su^2
     Z = X+rexp(n,sqrt(2)/su)-rexp(n,sqrt(2)/su)
     Y = 2*sin(pi*X)+(1+2*X)*e
   
     # Define the objective function
   
      h0 = n^(-1/7)    
      bandn = Naive.band.sel(Y,Z,su,h0,method="Naive")
      h1n = bandn[1]
      h2n = bandn[2]
      print(c(h1n,h2n))
   
      bandc = HSS.band.sel(Y,Z,su,h0,method="lap_cor")
      h1c = bandc[1]
      h2c = bandc[2]
      print(c(h1c,h2c))
   
      bandh = bandsel_huang(Y,Z,su)
      h1h = bandh[1]
      h2h = bandh[2]
      print(c(h1h,h2h))
     
      k = 1
      for(x0 in xseq)
         {
            # The Proposed Local Linear Method

             y1_in=mod0true[k]+0.1
             y2_in=mod1true[k]+0.1
             myoptcor = optimcor(Z,Y,y1_in,y2_in,x0,h1c,h2c)
    
     
              # The Naive Method
      
             objfun_naive = function(bt)
              {
                 mean(dnorm(Z,x0,h1n)*dnorm(Y,bt[1]+bt[2]*(Z-x0),h2n))
              }
      
            myoptnaive = optim(c(mod0true[k]+0.1,mod1true[k]+0.1),objfun_naive, control=list(fnscale=-1))$par
      
           # Huang's Local Linear Method
           
           # Mean-Shift Estimate   
      
            myopthuang = meanshift(Z,Y,mod0true[k],x0,h1h,h2h)[1]
      
            m0huang[j,k] = myopthuang
            m0cor[j,k] = myoptcor[1]
            m0naive[j,k] = myoptnaive[1]
            m1cor[j,k] = myoptcor[2]
            m1naive[j,k] = myoptnaive[2]
            k = k+1
         }
      cat(j,"\n")
   }

    
 tm=Sys.time()-btime
 
 save(file="data1.n=200.100.500.0.75.RData",list=ls())