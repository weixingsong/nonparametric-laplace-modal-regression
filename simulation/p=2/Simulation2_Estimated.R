# Simulation study for the paper
#  Nonparametric Modal Regression with Laplace Measurement Error
# 
# Simulation 2
#Choose n=200, 500. su=sqrt((1/0.75-1)*(4/3)), sqrt((1/0.90-1)*(4/3)).

 rm(list=ls())
 btime = Sys.time()

  source("p=2_bandwidth_selection.R")
  #source("bandsel_huang.R")
  #source("meanshift_huang.R")
  source("optimcor.R")
  source("p=2_CLy.bandwidth.R")
 

 # Data

  n=200
  gridlen=100
  reptime=500
  
  xseq=seq(-2,2,length=gridlen)  #huang

 
  mod0true = xseq+xseq^2     #huang
  mod1true = 1+2*xseq        #huang
  mod2true = rep(2,gridlen)        #huang
 
  mod0cor = mod0naive = mod0huang = rep(0,gridlen)
  mod1cor = mod1naive = rep(0,gridlen)
 
  m0cor = m0naive = m0huang = matrix(0,nrow=reptime,ncol=gridlen)
  m1cor = m1naive = matrix(0,nrow=reptime,ncol=gridlen)

  set.seed(12345)
  for(j in seq(reptime))
    {
      #data2  huang
      X = runif(n,-2,2)
      mx=X+X^2
      sigx=0.5+exp(-X^2)
      utemp = runif(n, 0, 1)  
      Y = (utemp<=0.5)*rnorm(n, mx - 2 * sigx, 2.5 * sigx)+(utemp>0.5)*rnorm(n, mx, 0.5 * sigx) 
      su = sqrt((1/0.75-1)*(4/3))
      #su = sqrt((1/0.90-1)*(4/3))
      su2 = su^2
      Z = X+rexp(n,sqrt(2)/su)-rexp(n,sqrt(2)/su)
   
      # Define the objective function

   
       h0 = n^(-1/7)    
       bandn = band.sel(Y,Z,su,h0,method="Naive")
       h1n = bandn[1]
       h2n = bandn[2]
       print(c(h1n,h2n))
   
       bandc = CLy_band_width(Z,Y,su,n) 

       h1c = bandc[1]
       h2c = bandc[2]
       print(c(h1c,h2c))   
      #bandh = bandsel_huang(Y,Z,su)
      #h1h = bandh[1]
      #h2h = bandh[2]
      #print(c(h1h,h2h))
     
      k = 1
      for(x0 in xseq)
         {
            # The Proposed Local Linear Method

             y1_in=mod0true[k]+0.1
             y2_in=mod1true[k]+0.1
             y3_in=mod2true[k]+0.1

             myoptcor = optimcor(Z,Y,y1_in,y2_in,y3_in,x0,h1c,h2c)
    
     
              # The Naive Method
      
             objfun_naive = function(bt)
              {
                 mean(dnorm(Z,x0,h1n)*dnorm(Y,bt[1]+bt[2]*(Z-x0)+bt[3]*(Z-x0)^2,h2n))
              }
      
            myoptnaive = optim(c(mod0true[k]+0.1,mod1true[k]+0.1,mod2true[k]+0.1),objfun_naive, control=list(fnscale=-1))$par
      
           # Huang's Local Linear Method
           
           # Mean-Shift Estimate   
      
            #myopthuang = meanshift(Z,Y,mod0true[k],x0,h1h,h2h)[1]
      
            #m0huang[j,k] = myopthuang
            m0cor[j,k] = myoptcor[1]
            m0naive[j,k] = myoptnaive[1]
            m1cor[j,k] = myoptcor[2]
            m1naive[j,k] = myoptnaive[2]
            k = k+1
         }
      cat(j,"\n")
   }




    
 tm=Sys.time()-btime 
 #save(file="data2.n=200.100.500.0.75.RData",list=ls())