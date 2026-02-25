# Simulation stuy for the paper "Parametric Modal Regression with Error Contaminated Covariates"
# Testing hypothesis: 
# Y~Gamma given link function and linear part are correctely specified.
#
# Test methods: MCCM1   
# Full data(D1), corresponding  Z=data[,2:7]  y=y
   #D2,  Z=data[y<5000,2:7]  y=y[y<5000]
   #D3,  Z=data[y<4000,2:7] y=y[y<4000]

tFFQD=function(B=50,M=100) 
 { 
   path=getwd()
   source(paste0(path,"\\case_study\\Dietary data\\testDMCCM.R"))
   source(paste0(path,"\\case_study\\Dietary data\\testDCK.R"))


  testFFQD<- matrix(NA, nrow = 3, ncol = 3)
 
  for(i in 1:3)
       {
         testFFQD[i, 1] <- testDMCCM(i,B,M)
         testCK <-testDCK(i,B)
         testFFQD[i, 2] <- testCK[1] 
         testFFQD[i, 3] <- testCK[2]
       }

  colnames(testFFQD) <- c("MCCM","CvM","KS")
  rownames(testFFQD) <- c("D1", "D2", "D3")

   cat("\n")

  cat("FFQ: D1-D3 test","\n")
  print(testFFQD)
  cat("\n")

}

btime=Sys.time()
sink("results/FFQDtest_results.txt")
tFFQD(50,100)
sink()
Sys.time()-btime


############## results #############
#FFQ: D1-D3 test 
#   MCCM  CvM   KS
#D1 0.01 0.00 0.01
#D2 0.06 0.00 0.02
#D3 0.53 0.02 0.01
##################################