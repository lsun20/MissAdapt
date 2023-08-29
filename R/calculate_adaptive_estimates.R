calculate_adaptive_estimates <- function(YR, VR, YU, VU, VUR) {
  
  # Read in the lookup table
  policy <- readMat("../Matlab/lookup_tables/policy.mat")
  const.policy <- readMat("../Matlab/lookup_tables/const_policy.mat")
  thresholds <- readMat('../Matlab/lookup_tables/thresholds.mat')
  const.thresholds <- readMat('../Matlab/lookup_tables/const_thresholds.mat')
  
  # Calculate the over-id test statistic
  YO <- YR - YU
  VO <- VR - 2 * VUR + VU
  VUO <- (VUR - VU)
  tO <- YO / sqrt(VO)
  GMM <- YU - VUO / VO * YO; V_GMM <- VU - VUO/VO*VUO;
  corr <- VUO / sqrt(VO) / sqrt(VU)
  
  # Interpolate the adaptive estimate based on the estimated correlation coeff
  Sigma_UO_grid <- abs(tanh(seq(-3, -0.05, 0.05)))
  Ky <- length(policy$y.grid)
  psi.grid <- rep(Ky, 0)
  for (i in seq(1, Ky, 1)) {
    psi.function <- splinefun(Sigma_UO_grid, policy$psi.mat[i,], method = "fmm", ties = mean)
    psi.grid[i] <- psi.function(abs(corr))
  }
  
  psi.grid.extrap <- splinefun(policy$y.grid, psi.grid, method = "natural")
  t_tilde <- psi.grid.extrap(tO)
  adaptive_nonlinear <- VUO / sqrt(VO) * t_tilde + GMM
  
  # Interpolate the adaptive estimate with a 20% constraint on the worst-case risk
  psi.const.grid <- rep(Ky, 0)
  for (i in seq(1, Ky, 1)) {
    psi.const.function <- splinefun(Sigma_UO_grid, const.policy$psi1.mat[i,], method = "fmm", ties = mean)
    psi.const.grid[i] <- psi.const.function(abs(corr))
  }
  
  psi.const.grid.extrap <- splinefun(policy$y.grid, psi.const.grid, method = "natural")
  t_tilde_const <- psi.const.grid.extrap(tO)
  adaptive_nonlinear_const <- VUO / sqrt(VO) * t_tilde_const + GMM
  
  # Interpolate the soft-threshold estimate based on the estimated correlation coeff
  st.function <- splinefun(Sigma_UO_grid, thresholds$st.mat, method = "fmm", ties = mean)
  st <- st.function(abs(corr))
  adaptive_st <- VUO / sqrt(VO) * ((tO > st) * (tO - st) + (tO < -st) * (tO + st)) + GMM
  
  results <- list(
    st = st,
    GMM = c(GMM, sqrt(V_GMM)),
    YO = c(YO, sqrt(VO)),
    adaptive_nonlinear = adaptive_nonlinear,
    adaptive_nonlinear_const = adaptive_nonlinear_const,
    adaptive_st = adaptive_st
  )
  
  return(results)
}
