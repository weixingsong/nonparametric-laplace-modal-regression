# Simulation study for the paper "Parametric Modal Regression with Error Contaminated Covariates"
# Btnumsum and NBtnumsum correspond to the proposed method and the Naive method. 

# Syntax:
#   table1(n=100,B=100,su2=0.25)

table1=function(n=100,B=100,su2=0.25)
 {
  set.seed(1234) 

  # mode-link function 
 
   g=function(x){exp(x)}
   d1g=function(x){exp(x)}
  
  # Simulated Data
   
   total=1000 
   d=2         # X1,X2 are contaminated with normal error
   nj=4        # replication at each (X1,X2)
   mx1=rep(0,d)   # mean of (X1,X2)
   mu=rep(0,d)     # mean of measurement error U
   covx1=diag(d)  # covariance of (X1,X2)
   covu=diag(rep(su2,d))  # covariance of u

########################Naive#############

  Nresult=matrix(0,nrow=total,ncol=5) 
  Mresult=matrix(0,nrow=total,ncol=5) 
  
  for(kk in seq(total))
   {
     X2<-rbinom (n,1,0.5)  # generate  n sample from X3
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
 
     U=mvrnorm(nj*n,mu,covu)  
     X1rep=X1[rep(seq_len(nrow(X1)),each =nj), ]  # repeat each row of X1 nj times
     Zrep=X1rep+U   # nj replicated observations at each row of X1

     Zmean=matrix(rep(0,n*d),nrow=n)                  # container of mean of Z's at each X1
     Zcov=matrix(rep(0,d*d*n),nrow=n,ncol=d*d)  # container of covariance of Z's at each X1

     
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
  
  # MCCL-1   
    
     ZmB=Zmean[rep(seq_len(nrow(Zmean)),each=B),] # repeat each row of Zmean B times  
     ZvB=Zcov[rep(seq_len(nrow(Zcov)),each =B), ] # repeat each row of Zcov B times 
     
     Tpm=mvrnorm(n*B,rep(0,nj-1),diag(nj-1))      # Generating Tp-values 
     Tpf=function(vv){vv[1:d]/sqrt(sum(vv^2))}
     Tp=t(apply(Tpm,1,Tpf))        
     
     Tp1=apply(ZvB[,1:d]*Tp,1,sum)          # Compute S*Tp
     Tp2=apply(ZvB[,(d+1):(2*d)]*Tp,1,sum)
     
     impart=sqrt((nj-1)/nj)*cbind(Tp1,Tp2)  # compute the imaginary part 
     realpart=ZmB                           # extract the real part
     
     X2rep=kronecker(X2,rep(1,B),"*")       # repeat each row of X2 B times

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

   BtM=apply(Mresult,2,mean)
   BtMD=apply(Mresult,2,median)
   BtS=apply(Mresult,2,sd)
   BtIQ=apply(Mresult,2,IQR)

   Btnumsum=round(rbind(BtMD,BtIQ),3)
   dimnames(Btnumsum)=list(c("Median", "IQR"),
                       c("beta0","beta1","beta2","beta3","phi"))
   

   outB=paste0(Btnumsum[1, ], "(", Btnumsum[2, ], ")")
   outN=paste0(NBtnumsum[1, ], "(", NBtnumsum[2, ], ")")
   rbind(outB,outN)
 }

# Using the following codes to generate Table 1

final=as.data.frame(matrix(NA,nrow=16,ncol=5))

jj=1
for(su2 in c(0.25,1))
 {
  for(n in c(100,200))  
   {
    for(B in c(100,200))
     {
      final[(1+2*(jj-1)):(2*jj),]=table1(n=n,B=B,su2=su2) 
      jj=jj+1
     }
   }
}

final0=final[-c(2,6,10,14),]
final1=final0[1:3,]
dimnames(final1)=list(c("100","200","Naive"),c("beta1","beta2","beta3","beta4","phi"))
final2=final0[4:6,]
dimnames(final2)=list(c("100","200","Naive"),c("beta1","beta2","beta3","beta4","phi"))
final3=final0[7:9,]
dimnames(final3)=list(c("100","200","Naive"),c("beta1","beta2","beta3","beta4","phi"))
final4=final0[10:12,]
dimnames(final4)=list(c("100","200","Naive"),c("beta1","beta2","beta3","beta4","phi"))

sink("results/table1.txt")

cat("##### Simulation Results \n\n")
cat("# Simulation Results in Table 1 \n")
cat("# \n")
cat("# n=100, su2=0.25 \n")
cat("# beta0         beta1         beta2         beta3          phi\n")
print(final1, quote = FALSE)
cat("# \n")
cat("# \n")
cat("# n=200, su2=0.25 \n")
cat("# beta0         beta1         beta2         beta3          phi\n")
print(final2, quote = FALSE)
cat("# \n")
cat("# \n")
cat("# n=100, su2=1 \n")
cat("# beta0         beta1         beta2         beta3          phi\n")
print(final3, quote = FALSE)
cat("# \n")
cat("# \n")
cat("# n=200, su2=1 \n")
cat("# beta0         beta1         beta2         beta3          phi\n")
print(final4, quote = FALSE)

sink()

#####  Simulation Results 

# Simulation Results in Table 1 
# 
# n=100, su2=0.25 
# beta0         beta1         beta2         beta3          phi
# 100   -0.241(0.118) -0.248(0.059) -0.249(0.057) -0.261(0.193) 5.248(1.155)
# 200   -0.251(0.122) -0.252(0.059) -0.251(0.063) -0.237(0.218) 5.341(1.306)
# Naive -0.226(0.122) -0.237(0.055) -0.236(0.059) -0.301(0.211) 5.059(1.208)
# 
# 
# n=200, su2=0.25 
# beta0         beta1         beta2         beta3          phi
# 100   -0.245(0.087) -0.249(0.039)  -0.25(0.043) -0.251(0.146) 5.171(0.886)
# 200   -0.251(0.084) -0.251(0.038) -0.253(0.041)  -0.24(0.137) 5.198(0.846)
# Naive -0.226(0.082) -0.236(0.036) -0.237(0.038)   -0.3(0.129) 4.943(0.797)
# 
# 
# n=100, su2=1 
# beta0         beta1         beta2         beta3          phi
# 100   -0.248(0.127)  -0.249(0.07)  -0.25(0.068) -0.251(0.238) 5.303(1.354)
# 200   -0.259(0.138) -0.256(0.074) -0.256(0.073) -0.229(0.253) 5.402(1.555)
# Naive -0.166(0.124) -0.201(0.055) -0.202(0.056) -0.444(0.208) 4.536(1.105)
# 
# 
# n=200, su2=1 
# beta0         beta1         beta2         beta3          phi
# 100   -0.249(0.098) -0.251(0.048)  -0.251(0.05) -0.246(0.169) 5.224(0.971)
# 200   -0.256(0.094) -0.252(0.044) -0.255(0.051) -0.231(0.162) 5.214(0.971)
# Naive -0.168(0.085)   -0.2(0.035) -0.202(0.038)  -0.44(0.128) 4.416(0.692)

