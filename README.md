# Self-tuning Rank Selection In Multivariate Response Regression

This package implements a self-tuning procedure of selecting the rank of the regression coefficient matrix in multivariate response regression models. The procedure can be used to select the rank of the factor components under a factor model. 

Example:  
 	Let Y and X be a n by m response matrix and a n by p feature design matrix generated as Y = XA + E. 
  SRS(Y, X) returns the estimated rank of the coefficient matrix A via regressing Y on X.
     
The following commands install the `STRS` package in R (require the installation of `devtools` package).
  ```
  library(devtools)
  install_github("bingx1990/STRS")
  ```

