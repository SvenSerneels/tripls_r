# tripls_r - Univariate Trilinear Partial Least Squares Regression
R package for univariate trilinear partial least squares regression, including case specific prediction intervals

Description
-----------
This package provides R code for univariate trilinear partial least squares regression (tri-PLS1)[1]. The method is a provides a combination of dimension reduction and prediction and is suited to predict a univariate dependent variabe from three-way tensor input data. Besides the estimator itself, the package also delivers inference for the estimates and the predictions based on them. Inference in this package is based on a first order linear approximation of the estimator, including an efficient algorithm to calculate the Jacobian matrices of the tri-PLS1 estimators theselves [2]. 

Overview
--------
The package provides: 
- a cross-validation function to determine the optimal number of latent components 
- the tri-PLS1 estimator 
- a summary function 
- a predict function for new data
- a plot function with several options, e.g. y vs. y predicted including prediction intervals.  

To do
-----
This package was written back in 2014 and was intended for distribution on CRAN. While the package is functional when installed locally, it does not comply with all CRAN standards. At some point, the corresponding adaptations may be included and the package may be distributed. 

References
----------
1. [Regression coefficients in multilinear PLS](https://onlinelibrary.wiley.com/doi/10.1002/%28SICI%291099-128X%28199801/02%2912%3A1%3C77%3A%3AAID-CEM496%3E3.0.CO%3B2-7), Sijmen de Jong, Journal of Chemometrics, 12 (1998), 77–81.
