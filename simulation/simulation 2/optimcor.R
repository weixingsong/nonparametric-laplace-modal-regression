# Estimate for
# Nonparametric Modal Regression With Laplace Measurement Error.

# this paper's Method

optimcor=function(Z,Y,y1_in,y2_in,x0,h1c,h2c)
   {
      iter=0
      tol=1e-6
     delta=1
     while(delta>tol&&iter<200)
        {  

            objfun0 = function(bt)
                {
                    W = Z-x0 
                     V = Y-(y1_in)-bt*(Z-x0)
                     temp = 1+0.5*su2*(h1c^(-2)+bt^2*h2c^(-2))-0.5*su2*(W*h1c^(-2)-bt*V*h2c^(-2))^2
                     mean(temp*dnorm(W,0,h1c)*dnorm(V,0,h2c))
                 }
       
             myopt1 = optim(c(y2_in),objfun0, control=list(fnscale=-1))$par

             objfun1 = function(bt)
                  {
                     W = Z-x0 
                     V = Y-bt-myopt1*(Z-x0)
                     temp = 1+0.5*su2*(h1c^(-2)+myopt1^2*h2c^(-2))-0.5*su2*(W*h1c^(-2)-myopt1*V*h2c^(-2))^2
                     mean(temp*dnorm(W,0,h1c)*dnorm(V,0,h2c))
                  }
       
             myopt0 = optim(c(y1_in),objfun1, control=list(fnscale=-1))$par
       
             delta=abs(myopt0-y1_in)
             y1_in=myopt0
             y2_in=myopt1
             iter=iter+1
          }
      return(c(y1_in,y2_in,iter))
   }
