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
          xh4 = xh1^4
          xh5 = xh1^5
          xh = cbind(1,xh1,xh2,xh3,xh4,xh5)

          myreg = lm(Y~xh1+xh2+xh3+xh4+xh5)
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
 
           m4=24*lmcof[5]+120*lmcof[6]*xh1 

           fhat=apply(matrix(dnorm(zdiff,0,h0),nrow=n),1,mean)
           fhat1=apply(matrix((-zdiff/h0^2)*dnorm(zdiff,0,h0),nrow=n),1,mean)


           vtild=1/(4*sqrt(pi))
           PPK<-function(t) (1.5^2-1.5*t^2+0.5^2*t^4)*((1/(sqrt(2*pi)))*exp(-(t^2)/2))^2    
           v0=integrate(PPK, lower = -Inf, upper = Inf)$value  ##(K'')^2

           K = vtild*v0*mean(g0/(g2^2*fhat))
           M = (1/8)^2*mean((m4+fhat1/fhat)^2)
           N = mean((g3/g2)^2)
           L = -(1/8)*mean((m4+fhat1/fhat)*g3/g2)

          yfun=function(h)
             {
                a=exp(h[1])
                b=exp(h[2])
                K/(n*a*b^3)+M*a^8+N*b^4+2*L*a^4*b^2
             }
          inva=n^(-1/12)
          sol=optim(c(inva,inva),yfun)$par
          h1 = exp(sol[1]) 
          h2 = exp(sol[2])
          return(c(h1,h2))
       }