calculate_coverage <- function(YR, VR, YU, VU, VUR) {
  
  # Read in the lookup table
  if(!file.exists("../Matlab/lookup_tables/policy.mat")){stop("Check if ../Matlab/lookup_tables/policy.mat exists")}
  if(!file.exists("../Matlab/lookup_tables/const_policy.mat")){stop("Check if ../Matlab/lookup_tables/const_policy.mat exists")}
  if(!file.exists("../Matlab/lookup_tables/thresholds.mat")){stop("Check if ../Matlab/lookup_tables/thresholds.mat exists")}
  if(!file.exists("../Matlab/lookup_tables/const_thresholds.mat")){stop("Check if ../Matlab/lookup_tables/const_thresholds.mat exists")}
  
  policy <- readMat("../Matlab/lookup_tables/policy.mat")
  
  # Calculate the over-id test statistic
  YO <- YR - YU
  VO <- VR - 2 * VUR + VU
  VUO <- (VUR - VU)
  tO <- YO / sqrt(VO)
  GMM <- YU - VUO / VO * YO; V_GMM <- VU - VUO/VO*VUO;
  corr <- VUO / sqrt(VO) / sqrt(VU); corr2 <- corr^2
  # Interpolate the adaptive estimate based on the estimated correlation coeff
  Sigma_UO_grid <- abs(tanh(seq(-3, -0.05, 0.05)))
  Ky <- length(policy$y.grid)
  psi.grid <- rep(Ky, 0)
  for (i in seq(1, Ky, 1)) {
    psi.function <- splinefun(Sigma_UO_grid, policy$psi.mat[i,], method = "fmm", ties = mean)
    psi.grid[i] <- psi.function(abs(corr))
  }
  
  psi.grid.extrap <- splinefun(policy$y.grid, psi.grid, method = "natural")
  # Interpolate the critical values for the adaptive estimator
  flci_adaptive_cv <-  readMat("../Matlab/lookup_tables/flci_adaptive_cv.mat")
  Sigma_UO_grid <- tanh(seq(-3, -0.05, 0.05))
  flci_c_function <-  splinefun(Sigma_UO_grid, flci_adaptive_cv$min.c.vec[1,], method = "fmm", ties = mean)#under B=0
  flci_c <- flci_c_function(corr)
  
  risk <-  readMat("../Matlab/lookup_tables/risk.mat")
  
  b_grid <- risk$b.grid;  
  Kb <- length(b_grid)
   # Use simulation to calculate the coverage over b_grid
  set.seed(1)
  sims <- 100000
  t <- rnorm(sims);  x <- rnorm(sims)
   
  stat <- matrix(0, nrow = sims, ncol = Kb)
  YRtest <- matrix(0, nrow = sims, ncol = Kb)
  
  for (j in 1:Kb) {
    t_b <- t + b_grid[j]
    t_tilde <- psi.grid.extrap(t_b)  # Lookup table to fill in the estimate
    # t_tilde <- t_b  # If just using Y_U
    stat[, j] <- abs(sqrt(1 - corr2) * x + corr * (t_tilde - b_grid[j]))
    YRtest[, j] <- abs(sqrt(1 - corr2) * x + corr * ((1+VO/VUO)*t - b_grid[j]))
  }
  
  rejection <- colMeans(stat > 1.96)
  # max_cov_sigmaU <- 1 - min(rejection)
  # min_cov_sigmaU <- 1 - max(rejection)
  max_cov_sigmaU <- 1 - max(rejection[abs(b_grid)<=1])
  min_cov_sigmaU <- 1 - max(rejection[abs(b_grid)<=2])
  
  rejection <- colMeans(stat > flci_c)
  # max_cov_flci <- 1 - min(rejection)
  # min_cov_flci <- 1 - max(rejection)
  max_cov_flci <- 1 - min(rejection[abs(b_grid)<=1])
  min_cov_flci <- 1 - max(rejection[abs(b_grid)<=2])
  
  rejection <- colMeans(YRtest > (1.96*sqrt(VR/VU)))
  # max_cov_YR  <- 1 - min(rejection)
  # min_cov_YR  <- 1 - max(rejection)
  max_cov_YR  <- 1 - min(rejection[abs(b_grid)<=1])
  min_cov_YR  <- 1 - max(rejection[abs(b_grid)<=2])
  
  results <- list(
    adaptive_sigmaU = c(min_cov_sigmaU, max_cov_sigmaU),
    adaptive_flci = c(min_cov_flci, max_cov_flci),
    flci_c = flci_c,
    YRtest =  c(min_cov_YR, max_cov_YR)
  )
  
  return(results)
}
  
  