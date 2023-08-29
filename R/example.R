# DESCRIPTION =========================================================
# Author: Liyang Sun
#
#  This script contains functions to calculate the adaptive estimators
#  described in Armstrong, Kline and Sun (2023) for scalar parameters.
#
#  First, this script contains function that are used to construct the
#  adaptive estimator given the input of robust and restricted estimators.
#
#  Second, this script contains function that calculates the 
#  worst-case adaptation regret.

# PRELIMINARIES =======================================================
library(R.matlab)
library(akima) # spline interpolation package

# Read in the lookup table -----------------------------------
policy <- readMat("../Matlab/lookup_tables/policy.mat")
const.policy <- readMat("../Matlab/lookup_tables/const_policy.mat")
thresholds <- readMat('../Matlab/lookup_tables/thresholds.mat')
const.thresholds <- readMat('../Matlab/lookup_tables/const_thresholds.mat')

# Read in the robust and restrictive estimates and their variance-covariance matrix
YR <- 2408.8413; VR <- 220.6057^2; # the restricted estimator and its squared standard error
YU <- 2217.2312; VU <- 257.3961^2; # the robust estimator and its squared standard error
VUR <- 211.4772^2; # the covariance between the restricted and robust esitmators

# Load the function for computing 
source("calculate_adaptive_estimates.R")

# Call the function with the provided inputs
adaptive_estimate_results <- calculate_adaptive_estimates(YR, VR, YU, VU, VUR)

# Print the results
print('The over-ID test statistic is based on')
adaptive_estimate_results$YO
print('The efficient estimate (under the restriction of no confounding bias) is')
adaptive_estimate_results$GMM
print('The adaptive estimate is')
print(adaptive_estimate_results$adaptive_nonlinear)
print('The adaptive estimate, constrained to increase the worst-case risk by no more than 20%, is')
print(adaptive_estimate_results$adaptive_nonlinear_const)
print('The adaptive soft-thresholded estimate is')
print(adaptive_estimate_results$adaptive_st)

tO <- adaptive_estimate_results$YO[1]/adaptive_estimate_results$YO[2]

# Define the values
row1 <- c(YU, YR, adaptive_estimate_results$YO[1], adaptive_estimate_results$GMM[1], 
          adaptive_estimate_results$adaptive_nonlinear, adaptive_estimate_results$adaptive_st, 
          YU*(abs(tO)>1.96)+YR*(abs(tO)<1.96))
row2 <- c(sqrt(VU), sqrt(VR), adaptive_estimate_results$YO[2], adaptive_estimate_results$GMM[2], NA, NA, NA)
row4 <- c( NA, NA, NA, NA, NA, NA,adaptive_estimate_results$st,1.96)


# Calculate the correlation coefficient, which determines the worst-case bias-variance tradeoff
VO <- VR - 2 * VUR + VU
VUO <- (VUR - VU)
corr <- VUO / sqrt(VO) / sqrt(VU)
print('The correlation coefficient is')
corr

source("calculate_max_regret.R")
max_regret_results <- calculate_max_regret(corr)

# Define the values
row3 <- c(max_regret_results$max_regret_YU-1, Inf, NA,Inf, max_regret_results$max_regret_nonlinear-1, 
          max_regret_results$max_regret_st-1, max_regret_results$max_regret_ttest-1)


# Create a matrix
table_data <- rbind(row1, row2, row3, row4)

# Create a LaTeX table using xtable package
library(xtable)

colnames(table_data) <- c("YU", "YR", "YO", "GMM", "Adaptive", "Soft-threshold", "Pre-test")

# Convert to LaTeX table format
latex_table <- xtable(table_data, caption = "Summary for Robustness Checks")

# Print the LaTeX table
print(latex_table, include.rownames = FALSE, hline.after = c(-1, 0, nrow(table_data)))


