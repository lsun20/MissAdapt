# DESCRIPTION =========================================================
#
# This script contains functions to calculate adaptive estimators as 
# described in Armstrong, Kline, and Sun (2023) for scalar parameters. 
# The Vignette offers a detailed step-by-step guide for this script. 
# We encourage you to explore it. If you encounter any issues, please 
# don't hesitate to reach out
#
# First, this script calls the function "calculate_adaptive_estimates.R" 
# to construct the adaptive estimator based on robust and restricted 
# estimators.
#
# Second, this script calls the function "calculate_max_regret.R" to 
# compute the worst-case adaptation regret.
#
# Third, this script invokes the function 
# "plot_adaptive_and_minimax_estimates.R" to generate plots illustrating 
# the locus of B-minimax estimates and associated risk functions.
#
# Ensure that the required lookup tables are stored in the correct 
# directory "../Matlab/lookup_tables/XXX," where XXX refers to various 
# Matlab objects provided by the authors. These objects contain 
# pre-tabulated solutions for the adaptive problem across a fine grid. 
# The functions in this script rely on these inputs to interpolate the 
# pre-tabulated solutions for calculating the adaptive estimator
# for specific applications.
#
# PRELIMINARIES =======================================================


rm(list=ls())
library(R.matlab)
library(akima) # spline interpolation package
library(pracma) # meshgrid function
library(xtable)
 
# LOAD DATA FROM THE APPLICATION AS IN VIGNETTE PART 1 ================
# Read in the robust and restrictive estimates and their variance-covariance matrix
YR <- 2408.8413; VR <- 220.6057^2; # the restricted estimator and its squared standard error
YU <- 2217.2312; VU <- 257.3961^2; # the robust estimator and its squared standard error
VUR <- 211.4772^2; # the covariance between the restricted and robust estimators

# COMPUTE THE ADAPTIVE ESTIMATOR AND OUTPUT THE TABLE SHOWN IN THE 
# VIGNETTE PART 3 =====================================================
# Load the function for computing 
source("calculate_adaptive_estimates.R")

# Call the function with the provided inputs
adaptive_estimate_results <- calculate_adaptive_estimates(YR, VR, YU, VU, VUR)

# Print the results
print('The over-ID test statistic is based on Y_O and its std err')
adaptive_estimate_results$YO
tO <- adaptive_estimate_results$YO[1]/adaptive_estimate_results$YO[2]
print('The efficient estimate (under the restriction of no confounding bias) is and its std err')
adaptive_estimate_results$GMM
print('The adaptive estimate is')
print(adaptive_estimate_results$adaptive_nonlinear)
print('The adaptive estimate, constrained to increase the worst-case risk by no more than 20%, is')
print(adaptive_estimate_results$adaptive_nonlinear_const)
print('The adaptive soft-thresholded estimate is')
print(adaptive_estimate_results$adaptive_st)

# Define the values
row1 <- c(YU, YR, adaptive_estimate_results$YO[1],
          adaptive_estimate_results$GMM[1], 
          adaptive_estimate_results$adaptive_nonlinear,
          adaptive_estimate_results$adaptive_st, 
          YU*(abs(tO)>1.96)+YR*(abs(tO)<1.96))
row2 <- c(sqrt(VU), sqrt(VR), adaptive_estimate_results$YO[2],
          adaptive_estimate_results$GMM[2], NA, NA, NA)
row4 <- c( NA, NA, NA, NA, NA,adaptive_estimate_results$st,1.96)


# Calculate the correlation coefficient, which determines the worst-case bias-variance tradeoff
VO <- VR - 2 * VUR + VU
VUO <- (VUR - VU)
corr <- VUO / sqrt(VO) / sqrt(VU)
print('The correlation coefficient is')
corr

source("calculate_max_regret.R")
max_regret_results <- calculate_max_regret(VR, VU, VUR) # needs the correlation coefficient only

# Define the values
row3 <- c(max_regret_results$max_regret_YU-1, Inf, NA,Inf,
          max_regret_results$max_regret_nonlinear-1, 
          max_regret_results$max_regret_st-1,
          max_regret_results$max_regret_ttest-1)

# Create a matrix
table_data <- rbind(row1, row2, row3, row4)

# Create a LaTeX table using xtable package
colnames(table_data) <- c("YU", "YR", "YO", "GMM", "Adaptive",
                          "Soft-threshold", "Pre-test")

# Convert to LaTeX table format
latex_table <- xtable(table_data, caption = "Summary for Robustness Checks")

# Print the LaTeX table
print(latex_table, include.rownames = FALSE, hline.after = c(-1, 0,
                                                             nrow(table_data)))

# PLOT THE FIGURE OF MINIMAX ESTIMATES AS SHOWN IN THE VIGNETTE PART 4
# =====================================================================
source("plot_adaptive_and_minimax_estimates.R")
plot_adaptive_and_minimax_estimates(YR, YU, VR, VU, VUR)

