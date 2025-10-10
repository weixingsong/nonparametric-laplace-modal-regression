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
                  X <- cbind(1, x, x^2, x^3) 
                  theta <- solve(t(X) %*% X) %*% t(X) %*% y
                  return(theta)
              }
        hattheta=least_squares_reg(hatx,y,n)
        haterr=y-hattheta[1]-hattheta[2]*hatx-hattheta[3]*hatx^2-hattheta[4]*hatx^3
 

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

        ######### m'(hatx),m''(hatx)
         phatm=hattheta[2]+2*hattheta[3]*hatx+3*hattheta[4]*hatx^2
         pphatm=2*hattheta[3]+6*hattheta[4]*hatx


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

        PL<-function(s) (s*(1/(sqrt(2*pi)))*exp(-(s^2)/2))^2 ##(L')^2
        tv=integrate(PL, lower = -Inf, upper = Inf)$value

        PPK<-function(t) ((t^2-1)*(1/(sqrt(2*pi)))*exp(-(t^2)/2))^2    ##PPK^2 
        v1=(1/4)*integrate(PPK, lower = -Inf, upper = Inf)$value  ##(K'')^2

        ttPPK<-function(t) t^2*((t^2-1)*(1/(sqrt(2*pi)))*exp(-(t^2)/2))^2    ##t^2*PPK^2
        ttPPKPK<-function(t) t*((t^2-1)*(1/(sqrt(2*pi)))*exp(-(t^2)/2))*(-t*(1/(sqrt(2*pi)))*exp(-(t^2)/2))    ##t*PPK^2*PK 
        v2=(1/2)*integrate(ttPPK, lower = -Inf, upper = Inf)$value+tv+integrate(ttPPKPK, lower = -Inf, upper = Inf)$value

        ttPPL<-function(t) t^2*((t^2-1)*(1/(sqrt(2*pi)))*exp(-(t^2)/2))    ##t^2*PPL
         v3=(1/2)*integrate(ttPPL, lower = -Inf, upper = Inf)$value
 
         tPL<-function(t) t*(-t*(1/(sqrt(2*pi)))*exp(-(t^2)/2)) ##tL'
         v4=(1/2)*integrate(tPL, lower = -Inf, upper = Inf)$value

          tttPL<-function(t) (t^3)*(-t*(1/(sqrt(2*pi)))*exp(-(t^2)/2)) ##t^3L'
          v5=(1/6)*integrate(tttPL, lower = -Inf, upper = Inf)$value

          tK<-function(t) (t^2)*(1/(sqrt(2*pi)))*exp(-(t^2)/2) ##t^2K
          mu2=integrate(tK, lower = -Inf, upper = Inf)$value


          V_2=(1/(v3*pphat_f_error))*ppphat_f_error*v5
          V_1=-(1/v3)*pphatm*v4*mu2
          K=(1/(hat_f_x*pphat_f_error*v3)^2)*(sigmau^4)*tv*v1*hatgz_x*hat_f_error

          M=mean(V_1^2)
          N=mean(V_2^2)
          L=mean(V_1*V_2)
          K=mean(K)

          delta=sqrt((sqrt(L^2+15*M*N)-L)/(5*N))
          h_1=((3*K)/(4*n*delta^5*(L+N*delta^2)))^(1/12)
          h_2=delta*h_1

          return(c(h_1,h_2))
      }
 
