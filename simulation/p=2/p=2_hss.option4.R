# Bandwidth Selection for
# Nonparametric Modal Regression With Laplace Measurement Error.

#  Correction-I Method


HSS.band.sel = function(Y,Z,su,h0,method="lap_cor")
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
                    fhat1=apply(matrix((-xdiff/h0^2)*dnorm(xdiff,0,h0),nrow=n),1,mean)
                 } else if(method=="nor_cor")
                       {
                          xh1 = mean(Z)*su2/var(Z)+(var(Z)-su2)*Z/var(Z)                      
                       }
        xh2 = xh1^2
        xh3 = xh1^3
        xh4 = xh1^4
        xh5 = xh1^5
        xh = cbind(1,xh1,xh2,xh3,xh4,xh5)

        #step 1

         myreg = lm(Y~xh1+xh2+xh3+xh4+xh5)
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
         hatJ=g21*crossprod(xh)/nrow(xh)  #  t(xh) %*% xh
         hatK=apply(g31*xh,2,mean)
         hatL=g01*crossprod(xh)/nrow(xh)
         h_opt=as.numeric(((3/(2*sqrt(pi)))/(hatK%*%solve(hatL)%*%hatK))^(1/7)*n^(-1/7))
   
         #step 4  

         # 3-order modal
   
         secf=function(bt)
            {
               mean(dnorm(Y-bt[1]-bt[2]*xh1-bt[3]*xh2-bt[4]*xh3-bt[5]*xh4-bt[6]*xh5,0,h_opt))
            }
         bt0=c(1,2,0.1,0.1,0.1,0.1)
         myparsec=optim(bt0,inif,control=list(fnscale=-1))$par
         res=Y-xh%*%myparsec

         #step 5  

         g0=mean(dnorm(res,0,h_opt))
         g2=mean(((res/h_opt)^2-1)*dnorm(res,0,h_opt)/h_opt^2)
         g3=mean((3*(res/h_opt)-(res/h_opt)^3)*dnorm(res,0,h_opt)/h_opt^3)

         m1=myparsec[2]+2*myparsec[3]*xh1+3*myparsec[4]*xh1^2+4*myparsec[5]*xh1^3+5*myparsec[6]*xh1^4    
         m4=24*myparsec[5]+120*myparsec[6]*xh1 

         zdiff=kronecker(xh1,xh1,"-")
         fhat=apply(matrix(dnorm(zdiff,0,h_opt),nrow=n),1,mean)
         fhat1=apply(matrix((-xdiff/h0^2)*dnorm(xdiff,0,h0),nrow=n),1,mean)

         #f(x-u)f(u)####
         zxdiff=kronecker(Z,xh1,"-")
         gz_x=apply(matrix(dnorm(zxdiff,0,h_opt),nrow=n),1,mean)


         vtild=1/(4*sqrt(pi))
         PPK1<-function(t) (0.5^2*t^8-4*t^6+15.5*t^4-17*t^2+6.25)*((1/(sqrt(2*pi)))*exp(-(t^2)/2))^2    
          v1=integrate(PPK1, lower = -Inf, upper = Inf)$value 

          
          eta = su^4*vtild*(1/4) 
          K1 = eta*v1*mean(gz_x*g0/(g2*fhat)^2)


          M = (1/8)^2*mean((m4+fhat1/fhat)^2)
          N = mean((g3/g2)^2)
          L = -(1/8)*mean((m4+fhat1/fhat)*g3/g2)

          yfun=function(h)
             {
                a=exp(h[1])
                b=exp(h[2])
                K1/(n*a^5*b^3)+M*a^8+N*b^4+2*L*a^4*b^2
             }
          inva=n^(-1/12)
          #sol=optim(c(inva,inva),yfun,method="BFGS")$par
          sol=optim(c(inva,inva),yfun)$par
          h1 = exp(sol[1]) 
          h2 = exp(sol[2])
          return(c(h1,h2))
       }