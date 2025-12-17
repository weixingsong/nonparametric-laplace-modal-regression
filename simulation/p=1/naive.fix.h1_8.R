# Bandwidth Selection for
# Nonparametric Modal Regression With Laplace Measurement Error.

# Naive Method


Naive.band.sel = function(Y,Z,su,h0,method="Naive")
     {
         # Initial Residual
         su2=su^2
         n=length(Y)
         if(method == "Naive")
             {
                xh1 = Z
                zdiff=kronecker(Z,Z,"-")
                fhat=apply(matrix(dnorm(zdiff,0,h0),nrow=n),1,mean)
              } else if(method=="lap_cor")
                   {
                       zdiff=kronecker(Z,Z,"-")
                       fhat=apply(matrix(dnorm(zdiff,0,h0),nrow=n),1,mean)
                       amatrix=t(matrix((zdiff>=0)-(zdiff<0),nrow=n))
                       ediff = matrix(exp(-abs(zdiff)*sqrt(2)/su),nrow=n)
                       ediff = ediff*amatrix
                       ed = apply(ediff,2,mean)
                       xh1= Z+ed/fhat
                       xdiff=kronecker(xh1,xh1,"-")
                       fhat=apply(matrix(dnorm(xdiff,0,h0),nrow=n),1,mean)
                   } else if(method=="nor_cor")
                         {
                             xh1 = mean(Z)*su2/var(Z)+(var(Z)-su2)*Z/var(Z)
                             zdiff=kronecker(Z,Z,"-")
                             fhat=apply(matrix(dnorm(zdiff,0,h0),nrow=n),1,mean)
                          }
          xh2 = xh1^2
          xh3 = xh1^3
          xh = cbind(1,xh1,xh2,xh3)

          #step 1

          myreg = lm(Y~xh1+xh2+xh3)
          iniires = resid(myreg)

          #step 2  
          h0=dpik(iniires)
  
          # Constant mode
  
          inif=function(bt)
             {
                mean(dnorm(iniires-bt,0,h0))
             }
  
          mypar = optim(c(1.2),inif,control=list(fnscale=-1),lower=0,upper=4,method="Brent")$par

          #step 3

          res1=iniires-mypar
          g01=mean(dnorm(res1,0,h0))
          g21=mean(((res1/h0)^2-1)*dnorm(res1,0,h0)/h0^2)
          g31=mean((3*(res1/h0)-(res1/h0)^3)*dnorm(res1,0,h0)/h0^3)
          hatJ=g21*crossprod(xh)/nrow(xh) 
          hatK=apply(g31*xh,2,mean)
          hatL=g01*crossprod(xh)/nrow(xh)
          h_opt=as.numeric(((3/(2*sqrt(pi)))/(hatK%*%solve(hatL)%*%hatK))^(1/7)*n^(-1/7))
   
          #step 4  
          # 3-order modal
   
          secf=function(bt)
               {
                        mean(dnorm(Y-bt[1]-bt[2]*xh1-bt[3]*xh2-bt[4]*xh3,0,h_opt))
                }
          bt0=c(1,2,0.1,0.1)
          myparsec=optim(bt0,inif,control=list(fnscale=-1))$par
          res=Y-xh%*%myparsec

  
          g0=mean(dnorm(res,0,h_opt))
          g2=mean(((res/h_opt)^2-1)*dnorm(res,0,h_opt)/h_opt^2)
          g3=mean((3*(res/h_opt)-(res/h_opt)^3)*dnorm(res,0,h_opt)/h_opt^3)
  
          m2=2*myparsec[3]+6*myparsec[4]*xh1 

          vtild=1/(4*sqrt(pi))
          v0=1/(2*sqrt(pi))
          K = vtild*v0*mean(g0/(g2^2*fhat))
          M = 0.5^2*mean(m2^2)
          N = 0.5^2*mean((g3/g2)^2)
          L = -0.5^2*mean(m2*g3/g2)


          h1=n^(-1/8)
          yfun=function(h)
             {
                b=h
                K/(n*h1*b^3)+M*h1^4+N*b^4+2*L*h1^2*b^2
             }
          inva=n^(-1/8)
          sol=optim(inva,yfun,lower=0,upper=4,method="Brent")$par
          h2 = sol
          return(c(h1,h2))
      }