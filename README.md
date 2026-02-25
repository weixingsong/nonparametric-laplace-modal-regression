# Parametric Modal Regression with Error-Contaminated Covariates

This repository contains all R codes and data to reproduce the numerical results from the manuscript "Parametric Modal Regression with Error Contaminated Covariates"
by Yanfei He, Jianhong Shi, and Weixing Song.

##  Overview

This study develops parametric modal regression methods for datasets with error-contaminated covariates. The repository includes:

- R codes to conduct all the simulation studies and real data analysis validating the proposed methods.
- A master R-file to source all the R codes.
- A folder containing all the intermediate results to quickly reproduces all the tables and figures in the paper.
- A folder contains all the tables and figures in the paper.
- An excel file reporting the running time for simulation studies. Note that the running time for Tables 4, 5, 6, 7, 8, 9 ,10 and Figures 9, 10, 11, 12 is short. 

Please refer to the repository structure described below for details.
 
##  Environment

R version 4.4.0 (2024-04-24 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 10 x64 (build 19045)
Matrix products: default

locale:

[1] LC_COLLATE=Chinese (Simplified)_China.utf8 
[2] LC_CTYPE=Chinese (Simplified)_China.utf8  
[3] LC_MONETARY=Chinese (Simplified)_China.utf8
[4] LC_NUMERIC=C                               
[5] LC_TIME=Chinese (Simplified)_China.utf8    

time zone: Asia/Shanghai
tzcode source: internal

attached base packages:

[1] splines   stats4    stats     graphics  grDevices utils     datasets 
[8] methods   base     

other attached packages:

 [1] plot3Drgl_1.0.4      plot3D_1.4.1         rgl_1.3.1           
 [4] scatterplot3d_0.3-44 VGAM_1.1-11          dcov_0.1.1          
 [7] MultiRNG_1.2.4       pracma_2.4.4         latex2exp_0.9.6     
[10] gofgamma_1.0         readxl_1.4.3         matrixcalc_1.0-6   
[13] smoothmest_0.1-3     MASS_7.3-60.2       

loaded via a namespace (and not attached):

 [1] knitr_1.47        cli_3.6.2         xfun_0.45         rlang_1.1.3      
 [5] stringi_1.8.4     misc3d_0.9-1      jsonlite_1.8.8    glue_1.7.0       
 [9] htmltools_0.5.8.1 cellranger_1.1.0  base64enc_0.1-3   fastmap_1.2.0    
[13] lifecycle_1.0.4   stringr_1.5.1     compiler_4.4.0    htmlwidgets_1.6.4
[17] Rcpp_1.0.12       tcltk_4.4.0       digest_0.6.36     R6_2.5.1         
[21] magrittr_2.0.3    tools_4.4.0      

##  Repository Structure
├── results/    # File containing simulation results (case study generate from master.R, other from InterResult file)\
├── InterResult/    # File containing simulation Inter-results \
├── result_time.xlsx    # Directory containing simulation run time results\
├── master.R              # Master script for repeating the entire simulation study\
├── case_study/         # Real-world data case studies\
│ ├── AD data/         # Alzheimer's Disease (AD) data analysis section\
│ │ ├── tableAD.R    # Main estimation program (generates Tables 8-10)\
│ │ ├── testAD.R      # Main P-value calculations (primary hypothesis tests for AD data)\
│ │ ├── testCK.R      # P-value calculations for Cramer-von Mises (CvM) and Kolmogorov-Smirnov (KS) tests\
│ │ ├── testMCCM1.R   # P-value calculations for MCCM1 test\
│ │ ├── gammaTest.R   # Subroutine: Implementation of testCK\
│ │ ├── resampling.R    # Subroutine: Implementation of testCK\
│ │ └── BJadni.xls         # Original dataset for AD study (Excel format)\
│ └── Dietary data/      # Dietary intake data analysis section\
│ ├── tableD.R             # Main analysis program (generates Tables 4-6, Figures 9-11 from the paper)\
│ ├── table7.R            # Supplementary analysis (generates Table 7, Figure 12 from the paper)\
│ ├── testD.R             # Main P-value calculations (tests for hypotheses D1-D3)\
│ ├── testDCK.R           # P-value calculations for CvM and KS tests\
│ ├── testDMCCM.R     # P-value calculations for MCCM1 test\
│ ├── testTwo.R        # Main P-value calculations \
│ ├── testTwoCK.R         # P-value calculations for CvM and KS tests\
│ ├── testTwoMCCM.R   # P-value calculations for MCCM1 test\
│ ├── gammaTest.R         # Subroutine: Implementation of CvM and KS tests\
│ ├── resampling.R         # Subroutine: Implementation of CvM and KS tests\
│ └── wishreg.xls           # Original dataset for dietary study (Excel format)\
├── simulation 1/      # First simulation study (corresponds to Simulation 1 in the paper)\
│ ├── table1.R           # Code to generate results for Table 1\
│ ├── table2.R           # Code to generate results for Table 2\
│ ├── table3.R          # Code to generate results for Table 3\
│ └── figure2.R          # Code to generate results for Figure 2\
├── simulation 2/     # Second simulation study (corresponds to Simulation 2 in the paper)\
│ ├── figure3.R         # Code to generate Figure 3\
│ ├── figure4.R        # Code to generate Figure 4\
│ ├── figure5.R         # Code to generate Figure 5\
│ ├── figure6.R         # Code to generate Figure 6\
│ ├── figure71.R       # Code to generate the Model 2.1 part of Figure 7\
│ ├── figure72.R        # Code to generate the Model 2.2 part of Figure 7\
│ ├── figure73.R        # Code to generate the Model 2.3 part of Figure 7\
│ ├── figure8.R          # Code to generate the complete Figure 8\
│ ├── Figure8 CK.R     # Code to generate the CvM and KS test parts of Figure 8\
│ ├── Figure8 DC.R     # Code to generate the DC test part of Figure 8\
│ ├── Figure8 MCCM1.R  # Code to generate the MCCM1 test part of Figure 8\
│ └── Figure8 MCCM2.R  # Code to generate the MCCM2 test part of Figure 8\
├── simulation B1/        # Appendix B - First simulation study\
│ ├── tableB11.R          # Code to generate results for Table B.11\
│ ├── tableB12.R          # Code to generate results for Table B.12\
│ ├── tableB13.R         # Code to generate results for Table B.13\
│ └── figureB13.R      # Code to generate Figure B.13\
├── simulation B2/     # Appendix B - Second simulation study\
│ ├── figureB14.R      # Code to generate Figure B.14\
│ ├── figureB15.R        # Code to generate Figure B.15\
│ ├── figureB16.R       # Code to generate Figure B.16\
│ ├── figureB17.R         # Code to generate Figure B.17\
│ ├── figureB18l.R        # Code to generate the left half of Figure B.18\
│ ├── FigureB18 Model B1.R   # Code to generate the Model B.1 part of Figure B.18\
│ ├── FigureB18 Model B2.R    # Code to generate the Model B.2 part of Figure B.18\
│ ├── FigureB18 Model B3.R    # Code to generate the Model B.3 part of Figure B.18\
│ ├── FigureB18 Model B4.R     # Code to generate the Model B.4 part of Figure B.18\
│ ├── figureB18r.R         # Code to generate the right half of Figure B.18\
│ ├── FigureB18 CK.R     # Code to generate the CvM and KS test parts of Figure B.18\
│ ├── FigureB18 DC.R      # Code to generate the DC test part of Figure B.18\
│ ├── FigureB18 MCCM2.R    # Code to generate the MCCM2 test part of Figure B.18\
│ ├── Model B3.R         # Main test code for Model B.3\
│ ├── Model B3CK.R        # CvM and KS test code for Model B.3\
│ ├── Model B3DC.R       # DC test code for Model B.3\
│ └── Model B3M2.R       # MCCM2 test code for Model B.3\

##  Support

For questions or comments about this code, please contact Weixing Song via email: weixing@ksu.edu
