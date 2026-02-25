# Simulation study for the paper "Parametric Modal Regression with Error Contaminated Covariates"

# Testing hypothesis: 
# Y~Gamma given link function and linear part are correctely specified.
#
# Test methods: MCCM1, CvM, KS    
# Note: mode-link function have three type as below.  The covariance of U have four type.

testAD=function(B=50,M=100) 
 { 
   path=getwd()
   source(paste0(path,"\\case_study\\AD data\\testMCCM1.R"))
   source(paste0(path,"\\case_study\\AD data\\testCK.R"))

  tableMCCM1=tableCvM=tableKS <- matrix(NA, nrow = 4, ncol = 3)
  
  for(j in 1:4) 
    {  
      for(i in 1:3)
       {
         tableMCCM1[j, i] <- testMCCM1(B,M,i,j)
         
         testC = testCK(B,i,j)
         tableCvM[j, i] <- testC[1]
  
         tableKS[j, i] <- testC[2]
       }
    }
  colnames(tableMCCM1) <- colnames(tableCvM) <- colnames(tableKS) <- c("link1", "link2", "link3")
  rownames(tableMCCM1) <- rownames(tableCvM) <- rownames(tableKS) <- c("diag(0,0)","diag(0.03,0)","diag(0,0.16)","diag(0.03,0.16)")

   cat("\n")

  cat("MCCM1 test","\n")
  print(tableMCCM1)
  cat("\n")

  cat("CvM test","\n")
  print(tableCvM)
  cat("\n")

  cat("KS test","\n")
  print(tableKS)

}

btime=Sys.time()
sink("results/testAD_results.txt")
testAD(50,100)
sink()
Sys.time()-btime

###########  results  #######
#MCCM1 test 
#                link1 link2 link3
#diag(0,0)        0.02  0.02  0.02
#diag(0.03,0)     0.02  0.01  0.01
#diag(0,0.16)     0.04  0.04  0.03
#diag(0.03,0.16)  0.06  0.06  0.06
#
#CvM test 
#                link1 link2 link3
#diag(0,0)        0.41  0.40  0.37
#diag(0.03,0)     0.39  0.38  0.37
#diag(0,0.16)     0.36  0.38  0.39
#diag(0.03,0.16)  0.44  0.44  0.43
#
#KS test 
#                link1 link2 link3
#diag(0,0)        0.50  0.40  0.36
#diag(0.03,0)     0.47  0.49  0.45
#diag(0,0.16)     0.29  0.42  0.30
#diag(0.03,0.16)  0.43  0.48  0.56
###############################