# MissAdapt

This repository contains the code to replicate all numerical results in ["Adapting to Misspecification"](https://arxiv.org/pdf/2305.14265.pdf) by Armstrong, Kline and Sun (2023).  It also provides MATLAB and R code implementing the adaptive estimator, its soft-thresholding approximation, and their risk limited variants that we propose in the paper.  The vignette included here is also discussed in the Chamberlain online seminar recorded [here](https://youtu.be/JrDsCW-1h6A).


## Introduction
The outputs and the plots of the numerical examples are included in their respective folders, along with replication code for the original result.

All Matlab scripts for implementing the adaptive estimator are included in the folder `Matlab/`.   

All R scripts for implementing the adaptive estimator are included in the folder `R/`.  

## A vignette for example usage
Here we provide an example usage of the adaptive estimator in the context of:  

Dobkin, Carlos, Amy Finkelstein, Raymond Kluender, and Matthew J. Notowidigdo. 2018. "The Economic Consequences of Hospital Admissions." American Economic Review, 108 (2): 308-52.

In this example, the parameter of interest is the contemporaneous impact of hospitalization on out-of-pocket medical spending (dollar per year). The restricted estimator assumes away a confounding trend, and the robust estimator allows for a linear confounding trend. If there is indeed no confounding trend, the restricted estimator is more precise. However, we can not be sure whether there is no confounding trend or not. Instead, the adaptive estimator pools information across the restricted and and the robust estimates to arrive at a single estimate that balances efficiency and robustness.

This example can be implemented using either the Matlab or R script.  Please note that the script assumes that the `/lookup_tables/` directory, which contains the pre-tabulated adaptive estimators, is correctly downloaded and referenced in the provided paths.
	
### 1. Load data
We first load the typical data from robustness checks: the robust and restricted estimates for the impact of hospitalization on out-of-pocket medical spending. We also load their variance-covariance matrix. Note that the `VUR` is similar to `VR` because this setting is close to the Hausman setting where where `YR` is efficient under the restriction of no confounding bias. 
```r
YR <- 2409; VR <- 221^2; # the restricted estimator and its squared standard error
YU <- 2217; VU <- 257^2; # the robust estimator and its squared standard error
VUR <- 211^2; # the covariance between the restricted and robust esitmators
```
### 2. Calculate the Over-ID Test Statistic, Efficient estimator, and correlation coefficient
We calculate the over-identification (over-ID) test statistic to be `tO=1.2`:
```r
YO <- YR - YU; VO <- VR - 2*VUR + VU;
tO <- YO / sqrt(VO);
```
We also calculate the efficient estimate to be `GMM=2379`:
```r
GMM <- YU - VUO /VO * YO;
```
Finally, we compute the correlation coefficient to be `corr=-0.52`:
```r
VUO <- (VUR - VU); corr <- VUO/sqrt(VO)/sqrt(VU);
```

### 3. Calculate the Adaptive Estimate based on Interpolation
Based on the correlation coefficient, we interpolate adaptive estimator `psi.grid.extrap` based on pre-tabulated results in the `/lookup_tables/` directory. We then input `psi.grid.extrap` with the over-ID test statistic `tO`, which returns the adaptive estimate:
```r
t_tilde <- psi.grid.extrap(tO) 
adaptive_nonlinear <- VUO/sqrt(VO) * t_tilde + GMM
```
The adaptive estimate is `adaptive_nonlinear=2302`. The script also includes codes for approximating the adaptive estimator based on soft-thresholding. The various estimates are summarized in the following table. The metric of "max regret" describes how much worse the mean squared error of the adaptive estimator can be than an oracle estimator that uses prior knowledge of the magnitude of the bias in the restricted estimator $Y_{R}$.  

In this context, the high value of the max regret points out that the pre-test estimator, which alternates between using the restricted estimator $Y_{R}$ and robust estimator $Y_{U}$ based on the over-identification statistic, doesn't exhibit an optimal efficiency-robustness tradeoff. In other words, while it may perform well for some values of bias, there exist values of bias at which the pre-test estimator doesn't achieve good efficiency-robustness tradeoff.

In contrast, the adaptive estimator provides a sensible summary of $Y_{U}$ and $Y_{R}$. Its max regret is low, indicating a more favorable efficiency-robustness tradeoff in the estimator's performance. Essentially, when an estimator's max regret is close to zero, its performance mirrors that of having complete knowledge about the exact extent of confounding bias (oracle performance).


| Hospitalization Year=0 | $Y_{U}$    | $Y_{R}$ | $Y_O$  |   GMM   | Adaptive | Soft-threshold | Pre-test  |
|-----------|------------|---------|--------|---------|----------|-----------|-------|
| Estimate   | 2,217   | 2,409  | 192     | 2,379    | 2,302     | 2,287 |  2,409  |
|Std Error  | (257)   | (221)  | (160)   | (219)    |           |       |       |
|Max Regret | 38%     | ∞      |       |  ∞       | 15%        | 15%       | 68%   |
| Threshold  |        |         |        |         |          | 0.52      | 1.96  |



## Matlab scripts

`/lookup_tables/` contains pre-tabulated solutions to the adaptation problem over a grid of correlation coefficients.

`adaptive_estimate.m` : a function that takes the input of the pre-tabulated solutions and outputs the adaptive estimates via an instantaneous lookup table.

`const_estimate.m` : a function that takes the input of the pre-tabulated solutions and outputs the adaptive estimates via an instantaneous lookup table where the worst-case risk is constrained to be no more than {5%,15%,20%} of Y_U.  The default is set to 20%.

`minimax_locus_plot.m` : a function that takes the input of the pre-tabulated solutions and outputs the locus of B-minimax estimates via an instantaneous lookup table.  The oracle risk function (risk of each B-minimax estimator) as well as the risk function of the adaptive estimator are also plotted.

`example.m` : an example script that calls `adaptive_estimate.m` and `const_estimate.m` and output the adaptive results in a formatted csv file; it also calls `minimax_locus_plot.m` and output the B-minimax locus plot along with risk functions.

## R scripts

`example.R` : an example script that reads the pre-tabulated solutions and outputs the adaptive estimates via an instantaneous lookup table.
 
