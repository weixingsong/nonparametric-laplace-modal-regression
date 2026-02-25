# Simulation study for the paper "Parametric Modal Regression with Error Contaminated Covariates"
#  AD data analysis 
# Estimate parameters and standard (result：estimates of parameters，Bcov：standard errors).
#Note: mode-link function have three type as below.   The covariance of U have four type.

tableAD=function(B=50, link_type=1, covu_type=1) 
 {
  # mode-link function 
  
  if (link_type == 1) 
   {
    g <- function(x) { exp(x) }
    d1g <- function(x) { exp(x) }
    
   } else if (link_type == 2) 
      {
       g <- function(x) { 1/(1+exp(-x)) }  
       d1g <- function(x) { exp(-x)/(1+exp(-x))^2 }
    
      } else if (link_type == 3) 
         {
          g <- function(x) { exp(-exp(-x)) }
          d1g <- function(x) { exp(-exp(-x)) * exp(-x) }
    
         }   
  # Coviance
  if (covu_type == 1) 
   {
    covu <- diag(c(0, 0))
    covu_name <- "diag(0, 0)"
    
   } else if (covu_type == 2)
      {
       covu <- diag(c(0.03, 0))
       covu_name <- "diag(0.03, 0)"
    
      } else if (covu_type == 3) 
         {
           covu <- diag(c(0, 0.16))
           covu_name <- "diag(0, 0.16)"
    
         } else if (covu_type == 4) 
            {
              covu <- diag(c(0.03, 0.16))
              covu_name <- "diag(0.03, 0.16)"
    
            } 

  set.seed(6666)
  
  path=getwd()
  data<- read.csv(paste0(path,"\\case_study\\AD data\\BJadni.csv"), header=T)

  y=data$y
  X1=data[2:3]

  n=length(y)    # sample size
  d=2                 # X1 are contaminated with normal error
  nj=4                # replication at each (X1,X2)
  mu=rep(0,d)   # mean of measurement error U


  U=mvrnorm(nj*n,mu,covu)  
  X1rep=X1[rep(seq_len(nrow(X1)),each =nj), ]  # repeat each row of X1 nj times
  Zrep=X1rep+U   # nj replicated observations at each row of X1

  Zmean=matrix(rep(0,n*d),nrow=n)                 # container of mean of Z's at each X1
  Zsd=matrix(rep(0,d*d*n),nrow=n,ncol=d*d)   # container of covariance of Z's at each X1
 
  for(j in seq(n))
   {
    Zmean[j,]=apply(Zrep[(nj*(j-1)+1):(nj*j),],2,mean)
    # Spectral decomposition of covariance
    Temp1=eigen(cov(Zrep[(nj*(j-1)+1):(nj*j),]))
    Temp2=(Temp1$vectors)%*%diag(sqrt(Temp1$values))%*%t((Temp1$vectors))
    Zsd[j,]=c(Temp2)
   }

  #######################MCCL###########

  ZmB=Zmean[rep(seq_len(nrow(Zmean)),each=B),] # repeat each row of Zmean B times  
  ZvB=Zsd[rep(seq_len(nrow(Zsd)),each =B), ] # repeat each row of Zcov B times 
 
  Tpm=mvrnorm(n*B,rep(0,nj-1),diag(nj-1))      # Generating Tp-values 
  Tpf=function(vv){vv[1:d]/sqrt(sum(vv^2))}
  Tp=t(apply(Tpm,1,Tpf))        
   
  Tp1=apply(ZvB[,1:d]*Tp,1,sum)          # Compute S*Tp
  Tp2=apply(ZvB[,(d+1):(2*d)]*Tp,1,sum)
 
  impart=sqrt((nj-1)/nj)*cbind(Tp1,Tp2)  # compute the imaginary part 
  realpart=ZmB                           # extract the real part
 
 
 
  #  Monte-Carlo Corrected Log-Likelihood 
 
  mccl=function(bet)
  {
   rp=ZmB%*%bet[2:(d+1)]+bet[1]
   ip=impart%*%bet[2:(d+1)]
   Wb=complex(real=rp,imaginary=ip)
   logg=Re(apply(matrix(log(g(Wb)),nrow=B),2,mean))
   invg=Re(apply(matrix(1/g(Wb),nrow=B),2,mean))
   phi=bet[d+2]
   out1=n*(1+phi)*log(phi)-n*log(gamma(phi+1))+phi*sum(log(y))
   out2=-(1+phi)*sum(logg)-phi*sum(y*invg)
   MCCL=out1+out2
   return(-MCCL/n)
  }
 

  
 
  mcov=function(bet)
  {
    Wb01=complex(real=realpart[,1],imaginary=impart[,1])
    Wb02=complex(real=realpart[,2],imaginary=impart[,2])
    rp=ZmB%*%bet[2:(d+1)]+bet[1]
    ip=impart%*%bet[2:(d+1)]
    Wb=complex(real=rp,imaginary=ip)
    phi=bet[d+2]
    yrep=kronecker(y,rep(1,B),"*")
    
    P0seq= -(1+phi)*d1g(Wb)/g(Wb)+phi*yrep*d1g(Wb)/(g(Wb))^2
    P1seq1= -(1+phi)*Wb01*d1g(Wb)/g(Wb)+phi*yrep*Wb01*d1g(Wb)/(g(Wb))^2
    P1seq2= -(1+phi)*Wb02*d1g(Wb)/g(Wb)+phi*yrep*Wb02*d1g(Wb)/(g(Wb))^2
    P3seq= 1+log(phi)-digamma(phi)+log(yrep)-log(g(Wb))-yrep/g(Wb)
    
    Psi0=Re(apply(matrix(P0seq,nrow=B),2,mean))
    Psi11=Re(apply(matrix(P1seq1,nrow=B),2,mean))
    Psi12=Re(apply(matrix(P1seq2,nrow=B),2,mean))
    Psi3=Re(apply(matrix(P3seq,nrow=B),2,mean))
    Psi=rbind(Psi0,Psi11,Psi12,Psi3)
    return((Psi%*%t(Psi))/n)
   }
   
 
  result=round(optim(c(-0.71,-0.15,-0.08,2.93),mccl)$par,3)
 
  An=hessian(mccl,result)
  Bn=mcov(result)
  invAn=solve(An)
  Bcov=round(sqrt(diag(invAn%*%Bn%*%invAn/n)),3)

  pvalue1=2*(1-pnorm(abs(result[2]/Bcov[2])))
  pvalue2=2*(1-pnorm(abs(result[3]/Bcov[3])))
  pvalue=round(c(pvalue1,pvalue2),3)

  outA <- paste0(result, "(", Bcov, ")")
  tempM <- rep("", 4)
  tempM[2:3] <- paste0("[", pvalue, "]")
  output=noquote(rbind(outA,tempM))
  dimnames(output)=list(c(covu_name," "),
                        c("beta0","beta1","beta2","phi"))
  output
}  


  table8=table9=table10 <- matrix(NA, nrow = 8, ncol = 4)
  final_rownames <- character(8)
  
  for(i in 1:4) 
    {  
      table8[(2*i-1):(2*i), ] <- as.matrix(tableAD(50,1,i))
  
      table9[(2*i-1):(2*i), ] <- as.matrix(tableAD(50,2,i))
  
      table10[(2*i-1):(2*i), ] <- as.matrix(tableAD(50,3,i))

      final_rownames[2*i-1] <- rownames(tableAD(50,1,i))[1]
      final_rownames[2*i] <- "" 
    }
  table8=noquote(table8)
  table9=noquote(table9)
  table10=noquote(table10)
  colnames(table8) <- colnames(table9) <- colnames(table10) <- c("beta0", "beta1", "beta2", "phi")
  rownames(table8) <- rownames(table9) <- rownames(table10) <- final_rownames

sink("results/tables_8_9_10.txt")

cat("========================================\n")
cat("Simulation Results - Tables 8, 9, and 10\n")
cat("Generated on:", date(), "\n")
cat("========================================\n\n")

# Table 8
cat("Simulation Results in Table 8\n")
cat("============================\n")
print(table8)
cat("\n\n")

# Table 9
cat("Simulation Results in Table 9\n")
cat("============================\n")
print(table9)
cat("\n\n")

# Table 10
cat("Simulation Results in Table 10\n")
cat("=============================\n")
print(table10)
cat("\n\n")

cat("========================================\n")
cat("End of Tables\n")
cat("========================================\n")

sink()



###############################results####################
#Simulation Results in Table 8 
#                 beta0         beta1         beta2         phi         
#diag(0, 0)       -2.037(0.06)  -0.34(0.158)  -0.199(0.076) 2.965(0.563)
#                               [0.031]       [0.009]                   
#diag(0.03, 0)    -2.039(0.061) -0.353(0.172) -0.198(0.077) 2.969(0.562)
#                               [0.04]        [0.01]                    
#diag(0, 0.16)    -2.036(0.06)  -0.338(0.156) -0.204(0.085) 2.965(0.57) 
#                               [0.03]        [0.016]                   
#diag(0.03, 0.16) -2.035(0.061) -0.33(0.163)  -0.201(0.085) 2.958(0.567)
#                               [0.043]       [0.018]                   
#
#Simulation Results in Table 9 
#                 beta0         beta1         beta2         phi         
#diag(0, 0)       -1.897(0.069) -0.392(0.181) -0.237(0.091) 2.968(0.566)
#                               [0.03]        [0.009]                   
#diag(0.03, 0)    -1.899(0.07)  -0.411(0.2)   -0.236(0.092) 2.973(0.564)
#                               [0.04]        [0.01]                    
#diag(0, 0.16)    -1.895(0.069) -0.388(0.18)  -0.244(0.103) 2.969(0.571)
#                               [0.031]       [0.018]                   
#diag(0.03, 0.16) -1.894(0.07)  -0.381(0.188) -0.241(0.103) 2.961(0.568)
#                               [0.043]       [0.019]                   
#
#Simulation Results in Table 10 
#                 beta0        beta1         beta2         phi         
#diag(0, 0)       -0.711(0.03) -0.168(0.077) -0.11(0.042)  2.975(0.571)
#                              [0.029]       [0.009]                   
#diag(0.03, 0)    -0.712(0.03) -0.179(0.086) -0.109(0.043) 2.979(0.569)
#                              [0.037]       [0.011]                   
#diag(0, 0.16)    -0.71(0.03)  -0.165(0.075) -0.116(0.05)  2.976(0.576)
#                              [0.028]       [0.02]                    
#diag(0.03, 0.16) -0.709(0.03) -0.164(0.079) -0.114(0.05)  2.97(0.569) 
#                              [0.038]       [0.023]            
##################################################################
