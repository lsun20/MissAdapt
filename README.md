# MissAdapt

This repository contains the code to replicate all numerical results in ["Adapting to Misspecification"](https://arxiv.org/pdf/2305.14265.pdf) by Armstrong, Kline, and Sun (2023).  It also provides MATLAB and R code implementing the adaptive estimator, its soft-thresholding approximation, and their risk limited variants proposed in the paper.  The vignette included here is also discussed in the Chamberlain online seminar recorded [here](https://youtu.be/JrDsCW-1h6A).


## Introduction
The output and the plots of the numerical examples are included in their respective folders, along with replication code for the original result.

All Matlab scripts for implementing the adaptive estimator are included in the folder `Matlab/`.   

All R scripts for implementing the adaptive estimator are included in the folder `R/`.  

## A vignette for example usage
Here we provide an example usage of the adaptive estimator to results from the paper:  

Dobkin, Carlos, Amy Finkelstein, Raymond Kluender, and Matthew J. Notowidigdo. 2018. "The Economic Consequences of Hospital Admissions." American Economic Review, 108 (2): 308-52.

In this example, the parameter of interest is the contemporaneous impact of hospitalization on out-of-pocket medical spending (in dollars per year). Two linear regression specifications are considered: an *unrestricted* specification that allows for a linear trend and a *restricted* specification omits this trend. These two specifications are estimated by ordinary least squares. Interpreting these results is difficult because we can't be sure whether a trend is actually present or not. If no trend is present, the restricted estimator will be more precise. However, if a linear trend is present, the restricted estimator will be biased, while the unrestricted estimator will remain unbiased. 

The adaptive estimator pools information across the unrestricted and restricted estimators to arrive at a single estimate that yields a nearly optimal balance between precision and worst case bias. By "nearly optimal" we mean that the adaptive estimator achieves a worst case mean squared error (MSE) as close as possible to that of an *oracle* who knows the magnitude (but not the sign) of the bias faced by the restricted estimator. We call the excess worst case MSE of the adaptive estimator over the oracle the *adaptation regret*. A small adaptation regret indicates a nearly optimal balance is being struck between precision and robustness. As discussed in the paper -- and illustrated below -- the adaptive estimator will always exhibit lower adaptation regret than than selecting a model based upon a pre-test, with especially large differences resulting when such tests exhibit low power.

While the vignette is written in `R`, this example can be implemented using either the `example.m` or `example.R` script.  Please note that the script assumes that the `/lookup_tables/` directory, which contains the pre-tabulated adaptive estimators, is correctly downloaded and referenced in the provided paths.
	
### 1. Load data
We first load the typical data from robustness checks: the unrestricted and restricted estimates for the impact of hospitalization on out-of-pocket medical spending. We also load their variance-covariance matrix. Note that the `VUR` is similar to `VR` because this setting is close to the Hausman setting where where `YR` is efficient under the restriction of no confounding bias. 
```r
YR <- 2409; VR <- 221^2; # the restricted estimator and its squared standard error
YU <- 2217; VU <- 257^2; # the unrestricted estimator and its squared standard error
VUR <- 211^2; # the covariance between the restricted and robust esitmators
```
### 2. Calculate the minimal set of ingredients needed for adaptation: over-ID test statistic and the correlation coefficient
We calculate the over-identification (over-ID) test statistic to be `tO=1.2`:
```r
YO <- YR - YU; VO <- VR - 2*VUR + VU;
tO <- YO / sqrt(VO);
```
We compute the correlation coefficient to be `corr=-0.52`:
```r
VUO <- (VUR - VU); corr <- VUO/sqrt(VO)/sqrt(VU);
```

### 3. Calculate the Adaptive Estimate based on Interpolation
Based on the correlation coefficient, we interpolate the adaptive estimator based on pre-tabulated results in the `/lookup_tables/` directory. To do so, we just need to feed our estimates to the wrapper function `calculate_adaptive_estimates()`, which returns the adaptive estimate. This function also approximates the adaptive estimator based on soft-thresholding, which can be thought of as analogous to LASSO shrinkage. The other wrapper function `calculate_max_regret()` returns the adaptation regret of the estimators, which we abbreviate as "max regret."   

```r
adaptive_estimate_results <- calculate_adaptive_estimates(YR, VR, YU, VU, VUR)
max_regret_results <- calculate_max_regret(corr)

```
The results returned by these two functions are summarized in the table below. In this case, the pre-test estimator, which chooses between the restricted estimator $Y_{R}$ and robust estimator $Y_{U}$ based on the over-identification statistic, exhibits a large max regret of 68%. Intuitively, while the pre-test may perform well if the bias is very large -- in which case $Y_{U}$ will be selected -- or very small -- in which case $Y_{R}$ will be selected -- there exist intermediate values of bias at which the pre-test estimator becomes very noisy because it has low power.

In contrast, the adaptive estimator provides exhibits a max regret of only 15%, indicating near oracle performance. That is, the adaptive estimator exposes the researcher to worst case MSE only 15% above what they would face if the magnitude of any confounding trend were known ex-ante.


| Hospitalization Year=0 | $Y_{U}$    | $Y_{R}$ | $Y_O$  |   GMM   | Adaptive | Soft-threshold | Pre-test  |
|-----------|------------|---------|--------|---------|----------|-----------|-------|
| Estimate   | 2,217   | 2,409  | 192     | 2,379    | 2,302     | 2,287 |  2,409  |
|Std Error  | (257)   | (221)  | (160)   | (219)    |           |       |       |
|Max Regret | 38%     | ∞      |       |  ∞       | 15%        | 15%       | 68%   |
| Threshold  |        |         |        |         |          | 0.52      | 1.96  |

Here the GMM estimate is the efficient estimate when `YR` is unbiased (i.e., when no trend is present).  It is calculated as:
```
GMM <- YU - VUO /VO * YO;
```

## Matlab scripts

`/lookup_tables/` contains pre-tabulated solutions to the adaptation problem over a grid of correlation coefficients.

`adaptive_estimate.m` : a function that takes the input of the pre-tabulated solutions and outputs the adaptive estimates via a lookup table.

`const_estimate.m` : a function that takes the input of the pre-tabulated solutions and outputs the adaptive estimates via a lookup table where the worst-case MSE is constrained to be no more than {5%,15%,20%} of $Y_U$.  The default is set to 20%.

`minimax_locus_plot.m` : a function that takes the input of the pre-tabulated solutions and outputs the locus of $B$-minimax estimates via a lookup table.  The oracle MSE function (risk of each $B$-minimax estimator) as well as the risk function of the adaptive estimator are also plotted.

`example.m` : an example script that calls `adaptive_estimate.m` and `const_estimate.m` and outputs the adaptive results in a formatted csv file; it also calls `minimax_locus_plot.m` and outputs the $B$-minimax locus plot along with risk functions.

## R scripts

`example.R` : an example script that reads the pre-tabulated solutions and outputs the adaptive estimates via an instantaneous lookup table.
 
