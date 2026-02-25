# Simulation study for the paper "Parametric Modal Regression with Error Contaminated Covariates"
#Btnumsum TBtnumsum and NBtnumsum correspond to the proposed method MCCL_1 MCCL_2 and the Naive method. 

tableB13=function(n=200,B=100,su2=0.25)
 {
  set.seed(1234)

  # mode-link function

   g=function(x){1/(1+exp(x))}
   d1g=function(x){-exp(x)/(1+exp(x))^2}

  # Simulated Data

   total=1000
   d=2         # X1,X2 are contaminated with normal error
   nj=4        # replication at each (X1,X2)
   mx1=rep(0,d)    # mean of (X1,X2)
   mu=rep(0,d)     # mean of measurement error U
   covx1=diag(d)   # covariance of (X1,X2)
   covu=diag(rep(su2,d))   # covariance of u

########################Naive#############

  Nresult=matrix(0,nrow=total,ncol=5)
  Mresult=matrix(0,nrow=total,ncol=5)
  Tresult=matrix(0,nrow=total,ncol=5)


  for(kk in seq(total))
   {
     X2<-rbinom (n,1,0.5)   # generate  n sample from X3
     X1=matrix(0,nrow=n,ncol=2)  # generate n sample from (X1,X2)
     for(i in seq(n))
      {
        X1[i,]<-mvrnorm(1,rep((X2[i]==1)-(X2[i]==0),d),covx1)
       }
     bt0=-0.25
     bt1=c(-0.25,-0.25)   # true beta-values for X1, X2
     bt2=-0.25            # true beta-value for X3
     phi=5                # true phi-value

    # generating sample from the gamma mode regression

     alph=1+phi
     bt=phi/g(bt0+X1%*%bt1+X2*bt2)
     y=rgamma(n,shape=alph,rate=bt)

    #Generating normal measurement error

     U=draw.multivariate.laplace(nj*n,d=2,gamma=2,mu=mu,Sigma=covu)
     X1rep=X1[rep(seq_len(nrow(X1)),each =nj), ]  # repeat each row of X1 nj times
     Zrep=X1rep+U   # nj replicated observations at each row of X1

     Zmean=matrix(rep(0,n*d),nrow=n)             # container of mean of Z's at each X1
     Zcov=matrix(rep(0,d*d*n),nrow=n,ncol=d*d)   # container of covariance of Z's at each X1


     for(j in seq(n))
      {
       Zmean[j,]=apply(Zrep[(nj*(j-1)+1):(nj*j),],2,mean)
       # Spectral decomposition of covariance
       Temp1=eigen(cov(Zrep[(nj*(j-1)+1):(nj*j),]))
       Temp2=(Temp1$vectors)%*%diag(sqrt(Temp1$values))%*%t((Temp1$vectors))
       Zcov[j,]=c(Temp2)
      }


    # Naive Likelihood

     Naivmccl=function(bet)
      {
       Wb=Zmean%*%bet[2:(d+1)]+as.vector(bet[d+2]*X2)+bet[1]
       logg=log(g(Wb))
       invg=1/g(Wb)
       phi=bet[d+3]
       out1=n*(1+phi)*log(phi)-n*log(gamma(phi+1))+phi*sum(log(y))
       out2=-(1+phi)*sum(logg)-phi*sum(y*invg)
       MCCL=out1+out2
       return(-MCCL/n)
      }
     Nresult[kk,]=optim(c(-0.2,-0.2,-0.2,-0.3,4),Naivmccl)$par

    # MCCL-2

     Zt=matrix(rep(0,n*d),nrow=n)
     for(j in seq(n)) 
       {
         Zt[j,]=X1[j,]+U[nj*(j-1)+1,]}
         Ztmean=Zt
         Temp1t=eigen(covu)
         Temp2t=(Temp1t$vectors)%*%diag(sqrt(Temp1t$values))%*%t((Temp1t$vectors))
         Ztcov=c(Temp2t)

         ZmBt=Ztmean[rep(seq_len(nrow(Ztmean)),each=B),]
         ZvBt=t(matrix(rep(Ztcov,n*B),nrow=d*d))

         Tpt=mvrnorm(n*B,rep(0,d),diag(d))      # Generating Tp-values 
   
         Tp1t=apply(ZvBt[,1:d]*Tpt,1,sum)          # Compute S*Tp
         Tp2t=apply(ZvBt[,(d+1):(2*d)]*Tpt,1,sum)
 
         impartt=cbind(Tp1t,Tp2t)                       # compute the imaginary part 
         realpartt=ZmBt                                      # extract the real part
 
         X2rep=kronecker(X2,rep(1,B),"*")       # repeat each row of X2 B times
 
 
    #  Monte-Carlo Corrected Log-Likelihood 
 
     Tmccl=function(bet)
       {
         rp=ZmBt%*%bet[2:(d+1)]+as.vector(bet[d+2]*X2rep)+bet[1]
         ip=impartt%*%bet[2:(d+1)]
         Wb=complex(real=rp,imaginary=ip)
         logg=Re(apply(matrix(log(g(Wb)),nrow=B),2,mean))
         invg=Re(apply(matrix(1/g(Wb),nrow=B),2,mean))
         phi=bet[d+3]
        out1=n*(1+phi)*log(phi)-n*log(gamma(phi+1))+phi*sum(log(y))
        out2=-(1+phi)*sum(logg)-phi*sum(y*invg)
        MCCL=out1+out2
        return(-MCCL/n)
      }

     Tresult[kk,]=optim(c(-0.3,-0.2,-0.2,-0.3,4),Tmccl)$par


    # MCCL-1

     ZmB=Zmean[rep(seq_len(nrow(Zmean)),each=B),] # repeat each row of Zmean B times
     ZvB=Zcov[rep(seq_len(nrow(Zcov)),each =B), ] # repeat each row of Zcov B times

     Tpm=mvrnorm(n*B,rep(0,nj-1),diag(nj-1))      # Generating Tp-values
     Tpf=function(vv){vv[1:d]/sqrt(sum(vv^2))}
     Tp=t(apply(Tpm,1,Tpf))

     Tp1=apply(ZvB[,1:d]*Tp,1,sum)         # Compute S*Tp
     Tp2=apply(ZvB[,(d+1):(2*d)]*Tp,1,sum)

     impart=sqrt((nj-1)/nj)*cbind(Tp1,Tp2)  # compute the imaginary part
     realpart=ZmB                           # extract the real part

    # Monte-Carlo Corrected Log-Likelihood

     mccl=function(bet)
     {
       rp=ZmB%*%bet[2:(d+1)]+as.vector(bet[d+2]*X2rep)+bet[1]
       ip=impart%*%bet[2:(d+1)]
       Wb=complex(real=rp,imaginary=ip)
       logg=Re(apply(matrix(log(g(Wb)),nrow=B),2,mean))
       invg=Re(apply(matrix(1/g(Wb),nrow=B),2,mean))
       phi=bet[d+3]
       out1=n*(1+phi)*log(phi)-n*log(gamma(phi+1))+phi*sum(log(y))
       out2=-(1+phi)*sum(logg)-phi*sum(y*invg)
       MCCL=out1+out2
       return(-MCCL/n)
     }
  
     Mresult[kk,]=optim(c(-0.2,-0.2,-0.2,-0.3,4),mccl)$par
   }


  # Naive Numerical summary of beta and phi

   NBtM=apply(Nresult,2,mean)
   NBtMD=apply(Nresult,2,median)
   NBtS=apply(Nresult,2,sd)
   NBtIQ=apply(Nresult,2,IQR)

   NBtnumsum=round(rbind(NBtMD,NBtIQ),3)
   dimnames(NBtnumsum)=list(c("Median","IQR"),
                              c("beta0","beta1","beta2","beta3","phi"))


   TBtM=apply(Tresult,2,mean)
   TBtMD=apply(Tresult,2,median)
   TBtS=apply(Tresult,2,sd)
   TBtIQ=apply(Tresult,2,IQR)

   TBtnumsum=round(rbind(TBtMD,TBtIQ),3)
   dimnames(TBtnumsum)=list(c("Median",
                         "IQR"),
                       c("beta0","beta1","beta2","beta3","phi"))


   BtM=apply(Mresult,2,mean)
   BtMD=apply(Mresult,2,median)
   BtS=apply(Mresult,2,sd)
   BtIQ=apply(Mresult,2,IQR)

   Btnumsum=round(rbind(BtMD,BtIQ),3)
   dimnames(Btnumsum)=list(c("Median", "IQR"),
                             c("beta0","beta1","beta2","beta3","phi"))

   outB=paste0(Btnumsum[1, ], "(", Btnumsum[2, ], ")")
   outT=paste0(TBtnumsum[1, ], "(", TBtnumsum[2, ], ")")
   outN=paste0(NBtnumsum[1, ], "(", NBtnumsum[2, ], ")")

   output=noquote(rbind(outB,outT,outN))
   dimnames(output)=list(c("MCCL1","MCCL2","Naive"),
                        c("beta0","beta1","beta2","beta3","phi"))
   cat("Simulation Results in Table B.13","\n")
   print(output)
 }

btime=Sys.time()
sink("results/tableB13.txt")
tableB13(200,100,0.25)
sink()
Sys.time()-btime


####################results#########################
#Simulation Results in Table B.13 
#      beta0         beta1         beta2         beta3         phi         
#MCCL1 -0.26(0.203)  -0.252(0.095) -0.252(0.094) -0.249(0.34)  5.132(0.757)
#MCCL2 -0.274(0.236) -0.26(0.114)  -0.25(0.107)  -0.229(0.365) 5.169(0.809)
#Naive -0.224(0.192) -0.236(0.089) -0.235(0.087) -0.305(0.322) 5.076(0.75) 
#> Sys.time()-btime
#Time difference of 2.641566 hours
####################################################