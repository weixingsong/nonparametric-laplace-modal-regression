# Simulation stuy for the paper "Parametric Modal Regression with Error Contaminated Covariates"
# Testing hypothesis: 
# Y~Gamma given link function and linear part are correctely specified.
#
# Test methods: MCCM1, CvM and KS  

testDT=function(B=50,M=100) 
 { 
   path=getwd()
   source(paste0(path,"\\case_study\\Dietary data\\testTwoMCCM.R"))
   source(paste0(path,"\\case_study\\Dietary data\\testTwoCK.R"))


  DTtest <- rep(0,3)
  DTtest[1] <- testDTMCCM(B,M) 
  testCK <- testDTCK(B)
  DTtest[2] <- testCK[1]  
  DTtest[3] <- testCK[2]   
  names(DTtest) <- c("MCCM","CvM","KS")
  
  cat("\n")
  cat("test: FFQ with two covariates ","\n")
  print(DTtest)
}

btime=Sys.time()
sink("results/FFQTtest_results.txt")
testDT(50,100)
sink()
Sys.time()-btime

#######  results  ##########
#test: FFQ with two covariates  
#MCCM  CvM   KS 
#0.54 0.10 0.19 
#######################