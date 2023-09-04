calculate_max_regret <- function(VR, VU, VUR) {
  VO <- VR - 2 * VUR + VU
  VUO <- (VUR - VU)
  corr <- VUO / sqrt(VO) / sqrt(VU)
  Sigma_UO_grid <- abs(tanh(seq(-3, -0.05, 0.05))) # grid of the correlation coefficient
  
  # Interpolate the risk function for the adaptive estimator
  if(!file.exists("../Matlab/lookup_tables/risk.mat")){stop("Check if ../Matlab/lookup_tables/risk.mat exists")}
  risk <-  readMat("../Matlab/lookup_tables/risk.mat")
  b_grid <- risk$b.grid; risk_mat<- risk$risk.mat;
  Kb <- length(b_grid)
  risk_function_adaptive <- numeric(Kb)
  
  for (i in 1:Kb) {
    psi.const.function <- splinefun(Sigma_UO_grid, risk_mat[i,], method = "fmm", ties = mean)
    risk_function_adaptive[i] <- psi.const.function(abs(corr))
  }
  
  # Calculate the oracle risk function
  if(!file.exists("../Matlab/lookup_tables/minimax_rho_B9.csv")){stop("Check if ../Matlab/lookup_tables/minimax_rho_B9.csv exists")}
  rho_tbl <- read.csv('../Matlab/lookup_tables/minimax_rho_B9.csv', header = FALSE)
  rho_b_over_sigma_function<- splinefun(rho_tbl[,1], rho_tbl[,2], method = "fmm", ties = mean)
  rho_b_over_sigma <- rho_b_over_sigma_function(abs(b_grid))
  risk_oracle <- rho_b_over_sigma + 1 / corr^2 - 1
  
  max_regret_YU = max((1/corr^2)/risk_oracle);
  message('The robust estimator Y_U has worst-case adaptation regret \n', round(max_regret_YU,2), '\n')
  
  max_regret_nonlinear <- max(risk_function_adaptive / risk_oracle)
  message('The adaptive estimator minimizes the worst-case adaptation to be \n', round(max_regret_nonlinear,2), '\n')
  
  # Add soft threshold risk functions
  Eb <- function(l) 1 + l^2 + (b_grid^2 - 1 - l^2) * (pnorm(l - b_grid) - pnorm(-l - b_grid)) + (-b_grid - l) * dnorm(l - b_grid) - (l - b_grid) * dnorm(-l - b_grid)
  
  if(!file.exists("../Matlab/lookup_tables/thresholds.mat")){stop("Check if ../Matlab/lookup_tables/thresholds.mat exists")}
  thresholds <- readMat('../Matlab/lookup_tables/thresholds.mat')
  # Interpolate the soft-threshold estimate based on the estimated correlation coeff
  st.function <- splinefun(Sigma_UO_grid, thresholds$st.mat, method = "fmm", ties = mean)
  st <- st.function(abs(corr))
  #message('The adaptive soft threshold is \n', round(st round(, '\n')
  
  risk_function_st_adaptive <- Eb(st) + 1 / corr^2 - 1
  
  max_regret_st <- max(risk_function_st_adaptive / risk_oracle)
  message('The adaptive soft-threshold estimator has worst-case adaptation regret\n', round(max_regret_st,2), '\n')
  
  # Use simulation to calculate the risk function for the pre-test estimator that switches between Y_U and Y_R
  set.seed(1)
  sims <- 100000
  x <- rnorm(sims)
  x_b <- matrix(x, sims, Kb) + matrix(b_grid, sims, Kb, byrow = TRUE)
  
  Ebsims_ht <- function(l) colSums(((x_b > l) * x_b + (x_b < l & x_b > -l) * (1 + VO / VUO) * x_b
                                    + (x_b < -l) * x_b - matrix(b_grid, sims, Kb, byrow = TRUE))^2) / sims
  
  risk_function_ht_ttest <- (Ebsims_ht(1.96) + 1 / corr^2 - 1)
  
  # Calculate max regret for various estimators
  max_regret_ttest <- max(risk_function_ht_ttest / risk_oracle)
  message('The pre-test estimator has worst-case adaptation regret\n', round(max_regret_ttest,2), '\n')
  
  results <- list(
    max_regret_YU = max_regret_YU,
    max_regret_nonlinear = max_regret_nonlinear,
    max_regret_st = max_regret_st,
    max_regret_ttest = max_regret_ttest
  )
  
  return(results)
}
