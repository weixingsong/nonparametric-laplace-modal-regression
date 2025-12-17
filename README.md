# Nonparametric Modal Regression with Laplace Measurement Error

This repository contains all code and data to reproduce the results from the manuscript "Nonparametric Modal Regression with Laplace Measurement Error"
by Yanfei He, Jianhong Shi, and Weixing Song.


#Overview
This study develops nonparametric modal regression methods for datasets with Laplace measurement errors in covariates. The repository includes:

Two simulation scenarios covering different estimator configurations

A real-data case study using Dietary data

Complete reproducible research files for all tables and figures

# Environment
R version: 4.4.0

Processor: Intel(R) Core(TM) i5-8500 CPU @ 3.00GHz

RAM: 8.00 GB

Matrix products: Default

# Required R Packages
 
 MASS, KernSmooth, readxl

#Repository Structure

simulation/

├── p=1/

│   ├── Simulation1_Estimated.R         # Simulation 1 results (Local linear modal regression)

│   ├── Simulation2_Estimated.R         # Simulation 2 results (Local linear modal regression)

│   ├── optimcor.R                      # Proposed estimation method (subroutine)

│   ├── meanshift_huang.R               # Zhou & Huang (2016) meanshift method (subroutine)

│   ├── hss.option4.R                        # Proposed bandwidth selection  (Corrected-I) (subroutine)

│   ├── bandsel_huang.R                  # CV-SIMEX bandwidth selector (subroutine)

│   ├── naive.option4.R                    # Naive-h method (subroutine)

│   ├── CLy.bandwidth.R                  # Proposed bandwidth selection (Corrected-II) (subroutine)

│   ├── bandwidth_selection.R         # Naive bandwidth selection (subroutine)

│   ├── hss.fix_h1_5.R              # Corrected-I with fixed h1= n^{-1/5} (subroutine)

│   ├── hss.fix_h1_8.R              # Corrected-I with fixed h1= n^{-1/8} (subroutine)

│   ├── naive.fix.h1_5.R           # Naive-h with fixed h1= n^{-1/5} (subroutine)

│   ├── naive.fix.h1_8.R           # Naive-h with fixed h1= n^{-1/8} (subroutine)

│   ├── data1.results                   # Results for Simulation 1

│   ├── data2.results                   # Results for Simulation 2

│   └── fix_h1.results                   # Results for Simulation 1 with fixed h1 variants 

│

└── p=2/

    ├── Simulation1_Estimated.R         # Simulation 1 results (Local quadratic modal regression)
    
    ├── Simulation2_Estimated.R         # Simulation 2 results (Local quadratic modal regression)
    
    ├── optimcor.R                      # Proposed estimation method (subroutine)
    
    ├── p=2_hss.option4.R               # Proposed bandwidth selection (Corrected-I) (subroutine)
    
    ├── p=1_naive.option4.R             # Naive-h method (subroutine)
    
    ├── p=2_CLy.bandwidth.R             # Proposed bandwidth selection (Corrected-II) (subroutine)
    
    ├── p=2_bandwidth_selection.R       #  Naive bandwidth selection 
    
    └── p=2.results                           # Results for  quadratic simulations


case_study/

└── Dietary_data/

    ├── estimate.full.R                 # Full estimation for Dietary data
    
    ├── estimate.boot.R                 # Bootstrap estimation for Dietary data
    
    ├── CLy.bandwidth.R                 # Proposed bandwidth selection (Corrected-II) (subroutine)
    
    ├── bandsel_huang.R                 # CV-SIMEX bandwidth selector (subroutine)
    
    ├── bandwidth_selection.R           # Naive_h bandwidth method (subroutine)
    
    ├── optimcor.R                      # Proposed estimation method (subroutine)
    
    ├── meanshift_huang.R               # Zhou & Huang (2016) meanshift method (subroutine)
    
    └── wishreg.xls                     # Dietary dataset
    


# Simulation Scenarios
p=1: Local linear modal regression 

p=2: Local quadratic modal regression

Sample sizes: n = {200, 500}

Reliability ratios: RR = {0.75, 0.90}

# Method Notation in Code
Corrected: Proposed method (LocPoly in manuscript)

Huang: Zhou & Huang (2016) LocLinear method

Corhuang: Combined approach (LocPL in manuscript)

Naive: Naive method ignoring measurement error

#Output

data*.results files containing simulation outputs

#Support
For questions or comments about this code, please contact:

Weixing Song - weixing@ksu.edu
