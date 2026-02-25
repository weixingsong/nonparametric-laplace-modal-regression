-- Readme --

 1. Intermediate Results of Table 2 are contained in table2.csv

    The structure of the table2.csv file

    1st row:  Estimation Methods
    2nd row:  hat_s.d. and s.d.
    3rd row:  The modal regression parameters
    4th row -- 1003th row: The 1000 estimates of parameters

 2. The R-script to generate Table 2 from the intermediate results table2.csv

    The output includes:

       mcclhatsd:  hat_sd from MCCL method
       mcclsd:     sd from MCCL method
       naivehatsd: hat_sd from Naive method
       naivesd:    sd from Naive method