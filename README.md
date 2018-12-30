# tripls_r - Univariate Trilinear Partial Least Squares Regression
R package for univariate trilinear partial least squares regression, including case specific prediction intervals

Description
-----------
This package provides R code for univariate trilinear partial least squares regression (tri-PLS1). The method is a provides a combination of dimension reduction and prediction. It is suited to predict a univariate dependent variable from three-way tensor input data. The method was first introduced by Ståhle [1]. A more efficient algorithm was later introduced for the univariate case [2], which is implemented here. 

Besides the estimator itself, the package also delivers inference for the estimates and the predictions based on them. Inference in this package is based on a first order linear approximation of the estimator, including an efficient algorithm to calculate the Jacobian matrices of the tri-PLS1 estimators theselves as proposed in [3]. 

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
This package was written back in 2014 and was intended for distribution on CRAN. While the package is functional when installed locally, it does not comply with all CRAN standards. At some point, the corresponding adaptations may be included and the package may be distributed through the official channels. 

References
----------
1. [Aspects of the analysis of three-way data](https://doi.org/10.1016/0169-7439(89)80114-5),Lars Ståhle, Chemometrics and Intelligent Laboratory Systems, 7 (1989), 95–100.
2. [Regression coefficients in multilinear PLS](https://onlinelibrary.wiley.com/doi/10.1002/%28SICI%291099-128X%28199801/02%2912%3A1%3C77%3A%3AAID-CEM496%3E3.0.CO%3B2-7), Sijmen de Jong, Journal of Chemometrics, 12 (1998), 77–81.
3. [Case specific prediction intervals for tri-PLS1: The full local linearisation](https://doi.org/10.1016/j.chemolab.2011.05.002), Sven Serneels, Klaas Faber, Tim Verdonck, Pierre J. Van Espen, Chemometrics and Intelligent Laboratory Systems, 108 (2011), 93–99.
