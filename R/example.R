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
YR <- 33.52740; VR <- 1.811236^2;
YU <- 52.94609; VU <- 2.539776^2;
VUR <- 1.811236^2; # BLP example

# Calculate the over-id test statistic
YO <- YR - YU;
VO <- VR - 2*VUR + VU;
VUO <- (VUR - VU);
print('The over-id test statistic is')
tO <- YO / sqrt(VO);
print(tO)
print('The efficient estimator is')
CUE <- YU - VUO /VO * YO;
print(CUE)
print('The correlation coefficient is')
corr <- VUO/sqrt(VO)/sqrt(VU);
print(corr)


# Interpolate the adaptive estimate based on the estimated correlation coeff
Sigma_UO_grid <- abs(tanh(seq(-3,-0.05,0.05)));
Ky <- length(policy$y.grid); psi.grid <- rep(Ky,0);
for (i in seq(1,Ky,1)) {
  psi.function <- splinefun(Sigma_UO_grid,policy$psi.mat[i,], method = "fmm", ties = mean)
  psi.grid[i] <- psi.function(abs(corr)) # interpolate across corr. coeff.
}

psi.grid.extrap <- splinefun(policy$y.grid,psi.grid, method = "natural")
t_tilde <- psi.grid.extrap(tO) # extrapolate for the adaptive estimate

adaptive_nonlinear <- VUO/sqrt(VO) * t_tilde + CUE

print('The adaptive estimate is')
print(adaptive_nonlinear)

# Interpolate the adaptive estimate with a 20% constraint on the worst-case risk 
Ky <- length(policy$y.grid); psi.const.grid <- rep(Ky,0);
for (i in seq(1,Ky,1)) {
  psi.const.function <- splinefun(Sigma_UO_grid,const.policy$psi1.mat[i,], method = "fmm", ties = mean)
  psi.const.grid[i] <- psi.const.function(abs(corr)) # interpolate across corr. coeff.
}

psi.const.grid.extrap <- splinefun(policy$y.grid,psi.const.grid, method = "natural")
t_tilde_const <- psi.const.grid.extrap(tO) # extrapolate for the adaptive estimate

adaptive_nonlinear_const <- VUO/sqrt(VO) * t_tilde_const + CUE

print('The adaptive estimate, constrained to increase the worst-case risk by no more than 20%,  is')
print(adaptive_nonlinear_const)

# Interpolate the soft-threshold estimate based on the estimated correlation coeff
print('The adaptive soft threshold is')
st.function <- splinefun(Sigma_UO_grid,thresholds$st.mat, method = "fmm", ties = mean)
st <- st.function(abs(corr))
print(st)
print('The adaptive soft-thresholded estimate is')
adaptive_st <- VUO/sqrt(VO) * ((tO > st)*(tO - st) + (tO < -st)*(tO + st)) + CUE
print(adaptive_st)

