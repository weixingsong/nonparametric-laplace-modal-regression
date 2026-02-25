# Simulation Study

# Note: To repeat the simulation, please set the the folder 
# where the master.R file is located as the working directory
  

# Install and load packages

load_or_install_packages <- function(packages) {
  # 1. Identify packages that are NOT yet installed
  new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  
  # 2. Install missing packages
  if(length(new_packages)) {
    message("Installing missing packages: ", paste(new_packages, collapse = ", "))
    install.packages(new_packages, dependencies = TRUE)
  }
  
  # 3. Load all packages using library()
  lapply(packages, library, character.only = TRUE)
  
  invisible(NULL)
}

load_or_install_packages (c("MASS", "smoothmest","matrixcalc","readxl","gofgamma",                      "latex2exp","pracma","MultiRNG","dcov","VGAM",
                            "scatterplot3d","plot3Drgl"))

path=getwd()

## Simulation 1

### Table 1
 source(paste0(path,"\\simulation 1\\table1.R"))

### Table 2
 source(paste0(path,"\\simulation 1\\table2.R"))
 
### Table 3
 source(paste0(path,"\\simulation 1\\table3.R"))
 
### Figure 2 
 source(paste0(path,"\\simulation 1\\figure2.R"))
 
 
## Simulation 2
 
### Figure 3 
 source(paste0(path,"\\simulation 2\\figure3.R"))
 
### Figure 4 
 source(paste0(path,"\\simulation 2\\figure4.R"))
 
### Figure 5 
 source(paste0(path,"\\simulation 2\\figure5.R"))
 
### Figure 6 
 source(paste0(path,"\\simulation 2\\figure6.R"))
 
### Figure 7
 source(paste0(path,"\\simulation 2\\figure71.R"))
 source(paste0(path,"\\simulation 2\\figure72.R"))
 source(paste0(path,"\\simulation 2\\figure73.R"))
 
### Figure 8 
 source(paste0(path,"\\simulation 2\\figure8.R"))

 
 
## Simulation Studies in Appendix B1
 
### Table B11
 source(paste0(path,"\\simulation B1\\tableB11.R"))
 
### Table B12
 source(paste0(path,"\\simulation B1\\tableB12.R"))
 
### Table B13
 source(paste0(path,"\\simulation B1\\tableB13.R"))

### Figure B13
 source(paste0(path,"\\simulation B1\\figureB13.R"))
 
 
## Simulation Study in Appendix B2
 
### Figure B14
 source(paste0(path,"\\simulation B2\\figureB14.R"))
 
### Figure B15
 source(paste0(path,"\\simulation B2\\figureB15.R"))
 
### Figure B16
 source(paste0(path,"\\simulation B2\\figureB16.R"))
 
### Figure B17
 source(paste0(path,"\\simulation B2\\figureB17.R"))
 
### Figure B18
 source(paste0(path,"\\simulation B2\\figureB18l.R"))
 source(paste0(path,"\\simulation B2\\figureB18r.R"))



## Case Study  

### Dietary data 

 #### Table 4-6, Figure 9-11
  source(paste0(path,"\\case_study\\Dietary data\\tableD.R"))

 #### Table 7, Figure 12
  source(paste0(path,"\\case_study\\Dietary data\\table7.R"))
 
### AD data  

#### Table 8-10
  source(paste0(path,"\\case_study\\AD data\\tableAD.R"))


#############   Supplement  ###############################

#  In application 1, the gamma assumption is tested by
#  the MCCM1, CvM, and KS procedure. The p-values can be
#  recovered by the following R-codes. 
#  (See the subsections 
#     Full data set (D1)
#     Data set with FFQ value less than 5000 (D2)
#     Data set with FFQ value less than 4000 (D3)
#   and the last paragraph in Application 1.)

# D1-D3
source(paste0(path,"\\case_study\\Dietary data\\testD.R"))

# Two covariates
source(paste0(path,"\\case_study\\Dietary data\\testTwo.R"))

#  In application 2, the gamma assumption is tested by
#  the MCCM1, CvM, and KS procedure. The p-values can be
#  recovered by the following R-codes. (See the third paragraph in
#  Application 2)

source(paste0(path,"\\case_study\\AD data\\testAD.R"))

#  Simulation Study in Appendix B2. The rejection rate of MCCM1,MCCM2,
#  DC,Cvm and KS in Model B.3 can be recovered by the following R-codes.

source(paste0(path,"\\simulation B2\\Model B3.R"))
