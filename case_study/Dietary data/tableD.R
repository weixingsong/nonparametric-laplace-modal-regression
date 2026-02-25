# Simulation study for the paper "Parametric Modal Regression with Error Contaminated Covariates"
#  Dietary data analysis 
# Estimate parameters and standard errors 
     #Method Notation in Code
    #MCCL : result(estimate of parameters)  and Bcov(estimate of standard errors),
    #Naive : Nresult(estimate of parameters)  and NBcov(estimate of standard errors)
# Full data(D1), corresponding  Z=data[,2:7]  y=y
   #D2,  Z=data[y<5000,2:7]  y=y[y<5000]
   #D3,  Z=data[y<4000,2:7] y=y[y<4000]


tableD=function(ytype=1,B=50)
 {
  set.seed(6666)

  # mode-link function

  g=function(x){exp(x)}
  d1g=function(x){exp(x)}

  ##Data retrieval
  path=getwd()
  data<- read.csv(paste0(path,"\\case_study\\Dietary data\\wishreg.csv"), header=T)

  if (ytype == 1) 
   {
    Z <- data[, 2:7]
    y <- data$ffq 
   } else if (ytype == 2) 
      {
        Z <- data[data$ffq < 5000, 2:7]
        y <- data$ffq[data$ffq < 5000]
       } else if (ytype == 3)  
          {
            Z <- data[data$ffq < 4000, 2:7]
            y <- data$ffq[data$ffq < 4000]  
           }
  y=scale(y,center=F,scale=T)
  Z=scale(Z,center=T,scale=T)
  Zrep=c(t(Z))

  n=length(y)   # sample size
  nj=ncol(Z)    # replication at each (X1,X2)

  Zmean=rep(0,n)  # container of mean of Z's at each X1
  Zsd=rep(0,n) # container of sd of Z's at each X1
  Zcov=rep(0,n) # container of sd of Z's at each X1
 
  for(j in seq(n))
   {
    Zmean[j]=mean(Zrep[(nj*(j-1)+1):(nj*j)])
    Zcov[j]=var(Zrep[(nj*(j-1)+1):(nj*j)])

    # Spectral decomposition of covariance
    Zsd[j]=sqrt(var(Zrep[(nj*(j-1)+1):(nj*j)]))
   }


  ########################Naive#############

  Naivmccl=function(bet)
   {
    Wb=Zmean*bet[2]+bet[1]
    logg=log(g(Wb))
    invg=1/(g(Wb))
    phi=bet[3]
    out1=n*(1+phi)*log(phi)-n*log(gamma(phi+1))+phi*sum(log(y))
    out2=-(1+phi)*sum(logg)-phi*sum(y*invg)
    MCCL=out1+out2
    return(-MCCL/n)
   }

  Naivmcov=function(bet)
   {
    Wb01=Zmean
    Wb=Zmean*bet[2]+bet[1]
    phi=bet[3]
    
    Psi0= -(1+phi)*d1g(Wb)/g(Wb)+phi*y*d1g(Wb)/(g(Wb))^2
    Psi11= -(1+phi)*Wb01*d1g(Wb)/g(Wb)+phi*y*Wb01*d1g(Wb)/(g(Wb))^2
    Psi3= 1+log(phi)-digamma(phi)+log(y)-log(g(Wb))-y/g(Wb)
   
    Psi=cbind(Psi0,Psi11,Psi3)
    return((t(Psi)%*%Psi)/n)
   }

  Nresult=optim(c(-1.58,0.27,5),Naivmccl)$par
 
  An=hessian(Naivmccl,Nresult)
  Bn=Naivmcov(Nresult)
  invAn=solve(An)
  NBcov=round(sqrt(diag(invAn%*%Bn%*%invAn/n)),3)



  #######################MCCL###########

  ZmB=rep(Zmean,each=B) # repeat each row of Zmean B times  
  ZvB=rep(Zsd,each =B)  # repeat each row of Zcov B times 
  
  Tpm=mvrnorm(n*B,rep(0,nj-1),diag(nj-1))      # Generating Tp-values 
  Tpf=function(vv){vv[1:1]/sqrt(sum(vv^2))}
  Tp=t(apply(Tpm,1,Tpf))        
   
  Tp=ZvB*Tp          # Compute S*Tp
 
  impart=sqrt((nj-1)/nj)*Tp   # compute the imaginary part 
  realpart=ZmB                           # extract the real part
 
 
 
  #  Monte-Carlo Corrected Log-Likelihood 
 
  mccl=function(bet)
  {
   rp=ZmB*bet[2]+bet[1]
   ip=impart*bet[2]
   Wb=complex(real=rp,imaginary=ip)
   logg=Re(apply(matrix(log(g(Wb)),nrow=B),2,mean))
   invg=Re(apply(matrix(1/(g(Wb)),nrow=B),2,mean))
   phi=bet[3]
   out1=n*(1+phi)*log(phi)-n*log(gamma(phi+1))+phi*sum(log(y))
   out2=-(1+phi)*sum(logg)-phi*sum(y*invg)
   MCCL=out1+out2
   return(-MCCL/n)
  }
 
  mcov=function(bet)
  {
    Wb01=complex(real=realpart,imaginary=impart)
    rp=ZmB*bet[2]+bet[1]
    ip=impart*bet[2]
    Wb=complex(real=rp,imaginary=ip)
    phi=bet[3]
    yrep=kronecker(y,rep(1,B),"*")
    
    P0seq= -(1+phi)*d1g(Wb)/g(Wb)+phi*yrep*d1g(Wb)/(g(Wb))^2
    P1seq1= -(1+phi)*Wb01*d1g(Wb)/g(Wb)+phi*yrep*Wb01*d1g(Wb)/(g(Wb))^2
    P3seq= 1+log(phi)-digamma(phi)+log(yrep)-log(g(Wb))-yrep/g(Wb)
    
    Psi0=Re(apply(matrix(P0seq,nrow=B),2,mean))
    Psi11=Re(apply(matrix(P1seq1,nrow=B),2,mean))
    Psi3=Re(apply(matrix(P3seq,nrow=B),2,mean))
    Psi=rbind(Psi0,Psi11,Psi3)
    return((Psi%*%t(Psi))/n)
  }
   
 
  result=optim(c(-1.58,0.27,5),mccl)$par
 
  An=hessian(mccl,result)
  Bn=mcov(result)
  invAn=solve(An)
  Bcov=round(sqrt(diag(invAn%*%Bn%*%invAn/n)),3)

  ################################HatX#################
  # Generating hatX1
  
  barZm=mean(Zmean)      ##sample mean of barZ
  barZmrep=rep(barZm,each=n) ##each row of barZ n times   
  biasZmbarZ=Zmean-barZmrep     ## sample bias  barZi-barZ
  hatsigmaz=var(Zmean)          ###hatsigmaz
  
  # Spectral decomposition of covariance hatsigmaz
   hatsigmazz=sqrt(solve(hatsigmaz))
  
  hatsigmau=mean(Zcov)/nj
  
  # choice sigmahat:(hatsigmaz-hatsigmau)+#
  hatsigmazu=round(hatsigmaz-hatsigmau,8)
  chatsigma=hatsigmazu*(hatsigmazu>0)
  # Spectral decomposition of covariance chatsigma
  choicsigmax=sqrt(chatsigma)
  
  ##sigmax/sigmaz
  choicsigma=choicsigmax*hatsigmazz
  
  # hatX1
  hatX1=barZmrep+c(choicsigma)*biasZmbarZ



  result=round(result,3)
  Nresult=round(Nresult,3)

  outB=paste0(result, "(", Bcov, ")")
  outN=paste0(Nresult, "(", NBcov, ")")

  output=noquote(rbind(outB,outN))
  dimnames(output)=list(c("MCCL","Naive"),
                        c("beta0","beta1","phi"))

  par(mfrow = c(1, 2))

  library(latex2exp)
  plot(Zmean,y,type="p",xlab=TeX(r"(Long Term Intake $\{\bar{Z}_{j}\}_{j=1}^{271}$)"),ylab="Scaled FFQ Intake")
  Zmean=sort(Zmean)
  Mode.mccl=(g(result[1]+Zmean*result[2]))
  Mode.naive=(g(Nresult[1]+Zmean*Nresult[2]))
  lines(Zmean,Mode.mccl,lty=2,lwd=2)
  lines(Zmean,Mode.naive,lty=3,lwd=2)

  plot(hatX1,y,type="p",xlab=TeX(r"(Long Term Intake $\{\hat{X}_{j}\}_{j=1}^{271}$)"),ylab="Scaled FFQ Intake")
  hatX1=sort(hatX1)
  hMode.mccl=(g(result[1]+hatX1*result[2]))
  hMode.naive=(g(Nresult[1]+hatX1*Nresult[2]))
  lines(hatX1,hMode.mccl,lty=2,lwd=2)
  lines(hatX1,hMode.naive,lty=3,lwd=2)

  print(output)
 }



sink("results/table4.txt")
cat("# Simulation Results in Table 4 \n")
cat("# \n")
tableD(1,50)
sink()
postscript("results/figure9.eps", width = 8, height = 6, 
           horizontal = FALSE, paper = "special")
tableD(1,50)  
dev.off()


sink("results/table5.txt")
cat("# Simulation Results in Table 5 \n")
cat("# \n")
tableD(2,50)
sink()
postscript("results/figure10.eps", width = 8, height = 6, 
           horizontal = FALSE, paper = "special")
tableD(2,50)  
dev.off()

sink("results/table6.txt")
cat("# Simulation Results in Table 6 \n")
cat("# \n")
tableD(3,50)
sink()
postscript("results/figure11.eps", width = 8, height = 6, 
           horizontal = FALSE, paper = "special")
tableD(3,50)  
dev.off()


