# Estimate for
# Nonparametric Modal Regression With Laplace Measurement Error.

# this paper's Method

optimcor=function(Z,Y,y1_in,y2_in,y3_in,x0,h1c,h2c)
   {
      iter=0
      tol=1e-6
     delta=1
     while(delta>tol&&iter<200)
        {  

            objfun2 = function(bt)
                {
                    W = Z-x0 
                    V = Y-(y1_in)-y2_in*(Z-x0)-bt*(Z-x0)^2
                    H = y2_in+2*bt*(Z-x0)
                    temp = 1-0.5*su2*(W^2/h1c^2-1)/h1c^2+su2*H*W*V/(h1c^2*h2c^2)
                            -su2*bt*V/h2c^2-0.5*su2*H^2*(V^2/h2c^2-1)/h2c^2
                    mean(temp*dnorm(W,0,h1c)*dnorm(V,0,h2c))
                 }
       
             myopt2 = optim(c(y3_in),objfun2, control=list(fnscale=-1))$par

            objfun1 = function(bt)
                {
                    W = Z-x0 
                    V = Y-(y1_in)-bt*(Z-x0)-myopt2*(Z-x0)^2
                    H = bt+2*myopt2*(Z-x0)
                    temp = 1-0.5*su2*(W^2/h1c^2-1)/h1c^2+su2*H*W*V/(h1c^2*h2c^2)
                            -su2*myopt2*V/h2c^2-0.5*su2*H^2*(V^2/h2c^2-1)/h2c^2
                    mean(temp*dnorm(W,0,h1c)*dnorm(V,0,h2c))
                 }
       
             myopt1 = optim(c(y2_in),objfun1, control=list(fnscale=-1))$par


             objfun0 = function(bt)
                  {
                    W = Z-x0 
                    V = Y-bt-myopt1*(Z-x0)-myopt2*(Z-x0)^2
                    H = myopt1+2*myopt2*(Z-x0)
                    temp = 1-0.5*su2*(W^2/h1c^2-1)/h1c^2+su2*H*W*V/(h1c^2*h2c^2)
                            -su2*myopt2*V/h2c^2-0.5*su2*H^2*(V^2/h2c^2-1)/h2c^2
                    mean(temp*dnorm(W,0,h1c)*dnorm(V,0,h2c))
                  }
       
             myopt0 = optim(c(y1_in),objfun0, control=list(fnscale=-1))$par
       
             delta=abs(myopt2-y3_in)
             y1_in=myopt0
             y2_in=myopt1
             y3_in=myopt2

             iter=iter+1
          }
      return(c(y1_in,y2_in,y3_in,iter))
   }
