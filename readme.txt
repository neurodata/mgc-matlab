The code is the Matlab version for fast, powerful, and universally consistent independence testing.
It implements Multiscale Graph Correlation and Distance Correlation for a number of metric and kernel choices (in particular, it also computes the Hilbert-Schmidt Independence Criterion)
It is able to test independence, do k-sample testing, do partial test, and do feature screening. 

The R version is available on CRAN and https://github.com/neurodata/r-mgc.
The Python version is available on PyPi and https://github.com/neurodata/hyppo.

Author: Cencheng Shen, shenc@udel.edu



Main reference papers:

Multiscale Graph Correlations: 
J. T. Vogelstein and Q. Wang and E. Bridgeford and C. E. Priebe and M. Maggioni and C. Shen, "Discovering and Deciphering Relationships Across Disparate Data Modalities", eLife, 2019.
C. Shen and C. E. Priebe and J. T. Vogelstein, "From Distance Correlation to Multiscale Graph Correlation", Journal of the American Statistical Association, 2020.

Distance Correlation (Biased, Unbiased, and Partial):
G. Szekely and M. Rizzo and N. Bakirov, "Measuring and Testing Independence by Correlation of Distances", Annals of Statistics, 2007.
G. Szekely and M. Rizzo, "Partial Distance Correlation with Methods for Dissimilarities", Annals of Statistics, 2014.

Fast Correlation Computation and Fast Testing:
A. Chaudhuri and W. Hu, "A Fast Algorithm for Computing Distance Correlation", 2018.
C. Shen and J. T. Vogelstein, "The Chi-Square Test of Distance Correlation", 2020.

Various Metric / Kernal Choices:
R. Lyons, "Distance Covariance in Metric Spaces", Annals of Probability, 2013.
A. Gretton and L. Gyorfi, "Consistent Nonparametric Tests of Independence", Journal of Machine Learning Research, 2010.
C. Shen and  J. T. Vogelstein, "The Exact Equivalence of Distance and Kernel Methods in Hypothesis Testing", 2020.

K-Sample Testing: 
M. Rizzo and G. Szekely, "Energy distance", Wiley Interdisciplinary Reviews: Computational Statistics, 2016.
C. Shen and  C. E. Priebe and J. T. Vogelstein, "The Exact Equivalence of Independence Testing and Two-Sample Testing", 2020.
