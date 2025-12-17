# Bandwidth Selection for
# Nonparametric Modal Regression With Laplace Measurement Error.

# Naive and Correction Method


   band.sel = function(Y,Z,su,h0,method="Naive")
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

          myreg = lm(Y~xh1+xh2+xh3)
          ires = resid(myreg)
          lmcof = myreg$coef

          # Constant mode
  
          inif=function(bt)
             {
                mean(dnorm(ires-bt,0,h0))
             }
  
           mypar = optim(c(1.2),inif,control=list(fnscale=-1),lower=0,upper=4,method="Brent")$par
           res=ires-mypar
   
  
           g0=mean(dnorm(res,0,h0))
           g2=mean(((res/h0)^2-1)*dnorm(res,0,h0)/h0^2)
           g3=mean((3*(res/h0)-(res/h0)^3)*dnorm(res,0,h0)/h0^3)
 
           m2=2*lmcof[3]+6*lmcof[4]*xh1 

           vtild=1/(4*sqrt(pi))
           v0=1/(2*sqrt(pi))
           K = vtild*v0*mean(g0/(g2^2*fhat))
           M = 0.5^2*mean(m2^2)
           N = 0.5^2*mean((g3/g2)^2)
           L = -0.5^2*mean(m2*g3/g2)

           delt = sqrt((sqrt(L^2+3*M*N)+L)/N)
           h1 = (3*K/(4*n*delt^5*(L+N*delt^2)))^(1/8)
           h2 = delt*h1
           return(c(h1,h2))
       }