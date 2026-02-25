-- Readme --

 1. Intermediate Results of Table B.11 are contained in tableB11.csv

    The structure of the tableB11.csv file

    1st row:  Covariance Structure
    2nd row:  The sample size (n) and bootstrap sample size (B)
    3rd row:  The modal regression parameters
    4th row -- 1003th row: The 1000 estimates of parameters

 2. The R-script to generate Table B.11 from the intermediate results tableB11.csv

    The output includes:

       mean: the average of 1000 estimates of parameters for each scenario
        std: the standard deviation of 1000 estimates of parameters for each scenario, which is the number within the parentheses.
