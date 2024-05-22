calculate_adaptive_estimates <- function(YR, VR, YU, VU, VUR) {
  
  # Read in the lookup table
  if(!file.exists("../Matlab/lookup_tables/policy.mat")){stop("Check if ../Matlab/lookup_tables/policy.mat exists")}
  if(!file.exists("../Matlab/lookup_tables/const_policy.mat")){stop("Check if ../Matlab/lookup_tables/const_policy.mat exists")}
  if(!file.exists("../Matlab/lookup_tables/thresholds.mat")){stop("Check if ../Matlab/lookup_tables/thresholds.mat exists")}
  if(!file.exists("../Matlab/lookup_tables/const_thresholds.mat")){stop("Check if ../Matlab/lookup_tables/const_thresholds.mat exists")}
  
  policy <- readMat("../Matlab/lookup_tables/policy.mat")
  const.policy <- readMat("../Matlab/lookup_tables/const_policy.mat")
  thresholds <- readMat('../Matlab/lookup_tables/thresholds.mat')
  const.thresholds <- readMat('../Matlab/lookup_tables/const_thresholds.mat')
  mse <- readMat('../Matlab/lookup_tables/emse_corr.mat')
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

  # Interpolate the hard-threshold estimate based on the estimated correlation coeff
  ht.function <- splinefun(Sigma_UO_grid, thresholds$ht.mat, method = "fmm", ties = mean)
  ht <- ht.function(abs(corr))
  adaptive_ht <- VUO / sqrt(VO) * ((tO > ht) * (tO ) + (tO < -ht) * (tO)) + GMM
  
  # Calculate the ERM estimate as well as the adaptive
  erm <- VUO / sqrt(VO) * tO * (tO^2/(tO^2+1)) + GMM
  lambda.function <-  splinefun(Sigma_UO_grid, mse$MSE.lambda.mat, method = "fmm", ties = mean)
  lambda <- lambda.function(abs(corr))
  adaptive_erm <- VUO / sqrt(VO) * tO * (tO^2/(tO^2+lambda)) + GMM
  results <- list(
    st = st, ht = ht,
    GMM = c(GMM, sqrt(V_GMM)),
    YO = c(YO, sqrt(VO)),
    adaptive_nonlinear = adaptive_nonlinear,
    adaptive_nonlinear_const = adaptive_nonlinear_const,
    adaptive_st = adaptive_st, adaptive_ht = adaptive_ht,
    erm = erm,
    erm_lambda = lambda,
    adaptive_erm = adaptive_erm
  )
  
  return(results)
}
