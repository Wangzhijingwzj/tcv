# Thining Cross-Validation
Thinning Cross-Validation for Determining the Number of Factors in Poisson Factor Models.

# Installation
```r
install.packages("devtools")  
devtools::install_github("Wangzhijingwzj/tcv")  
library(tcv) 
```

# Usage
```r
TCV <- multiDT(x, K = 5, rmax = 8)
```
* x: n by p matrix containing count responses, where n is the number of samples, p is the number of features
* K: The number of folds for conducting thinning cross-validation. Default as 5.
* rmax: The maximum number of factors for conducting thinning cross-validation. Default as 8.

# Value
* The TCV error and TICV error. The number of factors corresponding to the smallest result is TCV estimator.

