calculate_coverage <- function(YR, VR, YU, VU, VUR, B) {
  
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
  if (B == 0) {
    B_index <- 1
  } else {
    B_index <- which(flci_adaptive_cv$B.grid == B)
  }
  
  Sigma_UO_grid <- tanh(seq(-3, -0.05, 0.05))
  flci_c_function <-  splinefun(Sigma_UO_grid, flci_adaptive_cv$min.c.vec[B_index,], method = "fmm", ties = mean)#under B=0
  flci_c <- flci_c_function(corr)
  
  thresholds <- readMat('../Matlab/lookup_tables/thresholds.mat')
  # Interpolate the soft-threshold estimate based on the estimated correlation coeff
  st.function <- splinefun(Sigma_UO_grid, thresholds$st.mat, method = "fmm", ties = mean)
  st <- st.function(abs(corr))
  # Interpolate the critical values for the adaptive soft-threshold estimator
  flci_adaptive_st_cv <-  readMat("../Matlab/lookup_tables/flci_adaptive_st_cv.mat")
  flci_st_c_function <-  splinefun(Sigma_UO_grid, flci_adaptive_st_cv$min.st.c.vec[B_index,], method = "fmm", ties = mean) 
  flci_st_c <- flci_st_c_function(corr)
  
  risk <-  readMat("../Matlab/lookup_tables/risk.mat")
  
  b_grid <- risk$b.grid;  
  Kb <- length(b_grid)
   # Use simulation to calculate the coverage over b_grid
  set.seed(1)
  sims <- 100000
  t <- rnorm(sims);  x <- rnorm(sims)
   
  stat <- matrix(0, nrow = sims, ncol = Kb)
  st_stat <- matrix(0, nrow = sims, ncol = Kb)
  YRtest <- matrix(0, nrow = sims, ncol = Kb)
  pretest <- matrix(0, nrow = sims, ncol = Kb)
  
  for (j in 1:Kb) {
    t_b <- t + b_grid[j]
    t_tilde <- psi.grid.extrap(t_b)  # Lookup table to fill in the estimate
    # t_tilde <- t_b  # If just using Y_U
    stat[, j] <- abs(sqrt(1 - corr2) * x + corr * (t_tilde - b_grid[j]))
    t_st <- (t_b > st) * (t_b - st) + (t_b < -st) * (t_b + st)
    st_stat[, j] <- abs(sqrt(1 - corr2) * x + corr * (t_st - b_grid[j]))
    YRtest[, j] <- abs(sqrt(1 - corr2) * x + corr * ((1+VO/VUO)*t_b - b_grid[j]))
    pretest[, j] = (abs(t_b)>1.96) * (abs(sqrt(1 - corr2) * x + corr * (t_b - b_grid[j])) > 1.96 ) + (abs(t_b)<1.96) * (abs(sqrt(1 - corr2) * x + corr * ((1+VO/VUO)*t_b - b_grid[j])) > (1.96*sqrt(VR/VU)));
  }
 
  
  rejection <- colMeans(stat > 1.96)
  max_cov_sigmaU <- 1 - min(rejection)
  min_cov_sigmaU <- 1 - max(rejection)
  # max_cov_sigmaU <- 1 - max(rejection[abs(b_grid)<=1])
  # min_cov_sigmaU <- 1 - max(rejection[abs(b_grid)<=2])
  
  rejection <- colMeans(stat > flci_c)
  max_cov_flci <- 1 - min(rejection)
  min_cov_flci <- 1 - max(rejection)
  # max_cov_flci <- 1 - min(rejection[abs(b_grid)<=1])
  # min_cov_flci <- 1 - max(rejection[abs(b_grid)<=2])
  
  rejection <- colMeans(st_stat > 1.96)
  max_cov_st_sigmaU <- 1 - min(rejection)
  min_cov_st_sigmaU <- 1 - max(rejection)
  
  rejection <- colMeans(st_stat > flci_st_c)
  max_cov_st_flci <- 1 - min(rejection)
  min_cov_st_flci <- 1 - max(rejection)
  
  rejection <- colMeans(YRtest > (1.96*sqrt(VR/VU)))
  max_cov_YR  <- 1 - min(rejection)
  min_cov_YR  <- 1 - max(rejection)
  # max_cov_YR  <- 1 - min(rejection[abs(b_grid)<=1])
  # min_cov_YR  <- 1 - max(rejection[abs(b_grid)<=2])
  
  pretest_rejection <- colMeans(pretest)
  max_cov_pretest  <- 1 - min(pretest_rejection)
  min_cov_pretest  <- 1 - max(pretest_rejection)
  
  results <- list(
    adaptive_sigmaU = c(min_cov_sigmaU, max_cov_sigmaU),
    adaptive_flci = c(min_cov_flci, max_cov_flci),
    adaptive_st_sigmaU = c(min_cov_st_sigmaU, max_cov_st_sigmaU),
    adaptive_st_flci = c(min_cov_st_flci, max_cov_st_flci),
    flci_c = flci_c, flci_st_c = flci_st_c,
    YRtest =  c(min_cov_YR, max_cov_YR),
    pretest = c(min_cov_pretest, max_cov_pretest)
  )
  
  return(results)
}
  
  