# Bandwidth Selection for
# Nonparametric Modal Regression With Laplace Measurement Error.

# Huang's Method
 
  bandsel_huang = function(Y,Z,su)
     {
        # Initial Residual
  
        n=length(Y)
        h2=1.06*sd(Y)*n^(-1/5)
  
        B=15
  
        Ws=matrix(0,nrow=n,ncol=B)
        Wss=matrix(0,nrow=n,ncol=B)
  
        U=matrix(rexp(n*B,sqrt(2)/su)-rexp(n*B,sqrt(2)/su),nrow=n)
        V=matrix(rexp(n*B,sqrt(2)/su)-rexp(n*B,sqrt(2)/su),nrow=n)
        Ws=U+Z
        Wss=U+V+Z
  
        meanCp1=function(h1)
             {  
                 temp = 0
                     for(b in seq(B))
                           {
                               W=Ws[,b]
                               wdiff = kronecker(W,W,"-")
                                ydiff = kronecker(Y,Y,"-")
                                Kx=matrix(dnorm(wdiff,0,h1),nrow=n)
                                Ky1=matrix(dnorm(ydiff,0,h2),nrow=n)
                                Ky2=matrix(dnorm(ydiff,0,h2*sqrt(2)),nrow=n)
      
                                wL=quantile(W,0.025)
                                wU=quantile(W,0.975)
                                w=(W>=wL)*(W<=wU)
  
                                T1=rowSums((Kx %*% Ky2) * t(Kx))
                                T2=2*dnorm(0,0,h1)*rowSums(Kx * t(Ky2))
                                T3=(dnorm(0,0,h1))^2*dnorm(0,0,su*sqrt(2))
                                T4 = apply(Kx,2,sum)-dnorm(0,0,h1)
                                T5 =rowSums(Kx * t(Ky1))-dnorm(0,0,h1)*dnorm(0,0,h2)
  
                                temp=temp+mean(w*((T1-T2+T3)/T4^2)-2*T5/T4)
                            }
                  temp/B
              }
  
  
             h1s= optim(c(0.1),meanCp1,lower=0,upper=4,method="Brent")$par
  
             meanCp2=function(h1)
                  {  
                      temp = 0
                      for(b in seq(B))
                           {
                                W=Wss[,b]
                                wdiff = kronecker(W,W,"-")
                                 ydiff = kronecker(Y,Y,"-")
                                 Kx=matrix(dnorm(wdiff,0,h1),nrow=n)
                                 Ky1=matrix(dnorm(ydiff,0,h2),nrow=n)
                                 Ky2=matrix(dnorm(ydiff,0,h2*sqrt(2)),nrow=n)
      
                                 wL=quantile(W,0.025)
                                 wU=quantile(W,0.975)
                                 w=(W>=wL)*(W<=wU)
      
                                 T1=rowSums((Kx %*% Ky2) * t(Kx))
                                 T2=2*dnorm(0,0,h1)*rowSums(Kx * t(Ky2))
                                 T3=(dnorm(0,0,h1))^2*dnorm(0,0,su*sqrt(2))
                                 T4 = apply(Kx,2,sum)-dnorm(0,0,h1)
                                 T5 =rowSums(Kx * t(Ky1))-dnorm(0,0,h1)*dnorm(0,0,h2)
      
                                 temp=temp+mean(w*((T1-T2+T3)/T4^2)-2*T5/T4)
                            }
                       temp/B
                   }
  
            h1ss= optim(c(0.1),meanCp2,lower=0,upper=4,method="Brent")$par
            h1=h1s^2/h1ss  
  
            return(c(h1,h2))
      }

