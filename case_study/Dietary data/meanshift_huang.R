# Estimate for
# Nonparametric Modal Regression With Laplace Measurement Error.

# Huang's Method


Ku0=function(xv)
  {
     dnorm(xv)*(1-su2*(xv^2-1)/(2*h1h^2))
   }

Ku1=function(xv)
   {
       dnorm(xv)*((1+3*su2/(2*h1h^2))*xv-su2*xv^3/(2*h1h^2)) 
    }  

Ku2=function(xv)
   {
      dnorm(xv)*(-su2*xv^4/(2*h1h^2)+(1+5*su2/(2*h1h^2))*xv^2-su2/h1h^2) 
   }
 

meanshift=function(Z,Y,mod0true_k,x0,h1h,h2h)
   {
       iter=0
       tol=1e-6
       y_in=mod0true_k+0.1
       delta=1
       while(delta>tol&&iter<200)
          {  
             Sn1=mean(Ku1((Z-x0)/h1h)/h1h)
             Sn2=mean(Ku2((Z-x0)/h1h)/h1h)
             numerator=sum((Ku0((Z-x0)/h1h)*Sn2-Ku1((Z-x0)/h1h)*Sn1)*dnorm((Y-y_in)/h2h)*Y)
             denominator=sum((Ku0((Z-x0)/h1h)*Sn2-Ku1((Z-x0)/h1h)*Sn1)*dnorm((Y-y_in)/h2h))+1e-6
             y_k=numerator/denominator
             delta=abs(y_k-y_in)
              y_in=y_k
             iter=iter+1
          }
     return(c(y_k,iter))
   }
