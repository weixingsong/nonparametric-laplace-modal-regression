This README file provides supplementary information / reproducible research files for the manuscript 
Title: "Nonparametric Modal Regression with Laplace Measurement Error"

Authors: Yanfei He, Jianhong Shi, Weixing Song
Code was written by Yanfei He, Jianhong Shi and Weixing Song
In case of questions or comments please contact  weixing@ksu.edu.

The code was written/evaluated in R with the following software versions:
R version 4.4.0
Processor: Inter(R) Core(TM) i5-8500 CPU @ 3.00GHz 3.00GHZ 
Computer RAM: 8.00GB

Matrix products: default

attached  packages: MASS   KernSmooth readxl          

This repository contains the following  files that can be used to reproduce all table and figures of the manuscript.
It contains two subfolders containing the following files:

./case_study/:
This folder contains the following  files：
 
     ./Dietary data/
     This subfolder folder contains the following code(.R) and Dietary data(.xls)
         estimate.full.R
         This R program generates the left-panel of Figure 9

         estimate.boot.R
         This R program generates the right-panel of Figure 9

         CLy.bandwidth.R  bandsel_huang.R and  bandwidth_selection.R
         Subroutine required by estimate.full.R and estimate.boot.R. In the code, bandwidth_selection.R bandsel_huang.R 
         and CLy.bandwidth.R correspond to the proposed bandwidth selection method(Corrected-II) ,  CV-SIMEX bandwidth 
         selector and the Naive_h method in the manuscript, respectively.


         optimcor.R and meanshift_huang.R
         Subroutine required by estimate.full.R and estimate.boot.R.  In the code, optimcor.R and meanshift_huang.R correspond 
         to the proposed estimate method  and  meanshift method in Zhou and Huang (2016).

         wishreg.xls
         An XLS sheet containing Dietary data.
    
    
./simulation/
This folder contains the following two files, 

   ./simulation 1/
    This subfolder folder contains the following code(.R).
        Simulation1_Estimated.R
        All results of Table 1-2 and Figure 1-4 can be generated using the following parameter combinations: n = {200, 500},  RR = {0.75, 0.90}.  
 
       hss.option4.R  and naive.option4.R
       Subroutine required by Simulation1_Estimated.R. hss.option4.R  and naive.option4.R correspond to the 
       proposed bandwidth selection method(Corrected-I) and Naive-h method.

      bandsel_huang.R, optimcor.R and meanshift_huang.R
      As well as ./Dietary data/.

   ./simulation 2/
    This subfolder folder contains the following code(.R).
       Simulation2_Estimated.R
       All results of Table 3-4  and Figure 5-8 can be generated using the following parameter combinations: n = {200, 500},  RR = {0.75, 0.90}.
      

      CLy.bandwidth.R  bandsel_huang.R and  bandwidth_selection.R, optimcor.R and meanshift_huang.R
      As well as ./Dietary data/.
    
  results
  Each folder's result file contains the results generated from within that folder. Specifically, these results correspond to 
  the tables and figures presented in the manuscript. In the code, the Corrected , Huang, Corhuang and Naive correspond to the LocPoly, 
  LocLinear LocPL and Naive presented in the manuscript. 
