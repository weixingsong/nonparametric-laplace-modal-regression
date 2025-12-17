# Bandwidth Selection for
# Nonparametric Modal Regression With Laplace Measurement Error.

# Correction-II Method

 CLy_band_width<-function(z,y,sigmau,n) 
     {
        n=n
        ##############hatX##########
        hatX<-function(z,sigmau,n) 
               {
                  h=n^{-1/7}
                  hatX<-rep(0,n)
                  for(i in 1:n)
                       {
                           hatgz_i=sum((1/h)*dnorm((z - z[i])/h, sd = 1))
                           b=sigmau/sqrt(2)
                           hatl_z=(1/hatgz_i)*sum((exp(-abs(z-z[i])/b))*((z>=z[i])-(z<z[i])))
                           hatX[i]=z[i]+hatl_z
                       }   
                  return(hatX)
               }
        hatx=hatX(z,sigmau,n)
        hatx <- as.numeric(hatx)
      
        #m(x)=a0+a1*X+a2*X^2+a3*X^3

        least_squares_reg <- function(x, y, n) 
             {                
                  X <- cbind(1, x, x^2, x^3,x^4,x^5) 
                  theta <- solve(t(X) %*% X) %*% t(X) %*% y
                  return(theta)
              }
        hattheta=least_squares_reg(hatx,y,n)
        haterr=y-hattheta[1]-hattheta[2]*hatx-hattheta[3]*hatx^2
                  -hattheta[4]*hatx^3-hattheta[5]*hatx^4-hattheta[6]*hatx^5

 

        modal_linear_reg <- function(haterr, n)
             {
                   f2<-function(theta)
                            {
                                h=n^{-1/7}
                                K_h=dnorm(haterr-theta,0,h)
                                L=sum(K_h)
                                -L
                             }
                  theta=optim(c(1.2),f2,lower=0,upper=4,method="Brent")$par
                  return(theta)
             }
        hatmerr=modal_linear_reg(haterr,n)

        haterror=haterr-hatmerr

        ############f(x),f(epsilon),m(x)#####
        h=n^{-1/7}
        hat_f_error=(1/h)*mean(dnorm(haterror/h, sd = 1))
        phat_f_error=-(1/(h^3))*mean(haterror*dnorm(haterror/h, sd = 1))
        pphat_f_error=(1/(h^3))*mean(((haterror/h)^2-1)*dnorm(haterror/h, sd = 1))
        ppphat_f_error=(1/(h^4))*mean((3*(haterror/h)-(haterror/h)^3)*dnorm(haterror/h, sd = 1))

        f_x <- function(hatx,n)   
             {
                 fx <- sapply(hatx, function(x) 
                      {
                          h=n^{-1/7}
                          hat_f_x=(1/h)*mean(dnorm((hatx-x)/h, sd = 1))
                          hat_f_x
                       })
                 return(fx)
              }
        hat_f_x=f_x(hatx,n)

        f_x1 <- function(hatx,n)   
             {
                 fx1 <- sapply(hatx, function(x) 
                      {
                          h=n^{-1/7}
                          hat_f_x1=mean((-(hatx-x)/h^2)*(1/h)*dnorm((hatx-x)/h, sd = 1))
                          hat_f_x1
                       })
                 return(fx1)
              }
        hat_f_x1=f_x1(hatx,n)


        ######### m''''(hatx)
         m1=hattheta[2]+2*hattheta[3]*hatx+3*hattheta[4]*hatx^2+4*hattheta[5]*hatx^3+5*hattheta[6]*hatx^4    
         m4=24*hattheta[5]+120*hattheta[6]*hatx 


          ###f(x-u)f(u)####
        gz_x <- function(z,hatx,n)
           {
                 gx <- sapply(hatx, function(x) 
                          {
                              h=n^{-1/5}
                              phat_f_x=mean((1/h)*dnorm((z - x)/h, sd = 1))
                           })
                 return(gx)
            }

        hatgz_x=gz_x(z,hatx,n)

         vtild=1/(4*sqrt(pi))
         PPK1<-function(t) (0.5^2*t^8-4*t^6+15.5*t^4-17*t^2+6.25)*((1/(sqrt(2*pi)))*exp(-(t^2)/2))^2    
          v1=integrate(PPK1, lower = -Inf, upper = Inf)$value  ##(K'')^2
          
          eta = su^4*vtild*(1/4) 
          K1 = eta*v1*mean(hatgz_x*hat_f_error/(pphat_f_error*hat_f_x)^2)


          M = (1/8)^2*mean((m4+hat_f_x1/hat_f_x)^2)
          N = mean((ppphat_f_error/pphat_f_error)^2)
          L = -(1/8)*mean((m4+hat_f_x1/hat_f_x)*ppphat_f_error/pphat_f_error)

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
 
