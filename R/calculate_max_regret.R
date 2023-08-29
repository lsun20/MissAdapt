calculate_max_regret <- function(corr) {
  Sigma_UO_grid <- abs(tanh(seq(-3, -0.05, 0.05)))
  rho_tbl <- read.csv('../Matlab/lookup_tables/minimax_rho_B9.csv', header = FALSE)
  # Calculate the oracle risk function
  rho_b_over_sigma_function<- splinefun(rho_tbl[,1], rho_tbl[,2], method = "fmm", ties = mean)
  rho_b_over_sigma <- rho_b_over_sigma_function(abs(b_grid))
  risk_oracle <- rho_b_over_sigma + 1 / corr^2 - 1
  
  max_regret_YU = max((1/corr^2)/risk_oracle);
  
  # Risk function
  risk <-  readMat("../Matlab/lookup_tables/risk.mat")
  b_grid <- risk$b.grid; risk_mat<- risk$risk.mat;
  Kb <- length(b_grid)
  risk_function_adaptive <- numeric(Kb)
  
  for (i in 1:Kb) {
    psi.const.function <- splinefun(Sigma_UO_grid, risk_mat[i,], method = "fmm", ties = mean)
    risk_function_adaptive[i] <- psi.const.function(abs(corr))
  }
  
  max_regret_nonlinear <- max(risk_function_adaptive / risk_oracle)
  cat('The worst-case adaptation regret is\n', max_regret_nonlinear, '\n')
  
  # Add soft threshold risk functions
  Eb <- function(l) 1 + l^2 + (b_grid^2 - 1 - l^2) * (pnorm(l - b_grid) - pnorm(-l - b_grid)) + (-b_grid - l) * dnorm(l - b_grid) - (l - b_grid) * dnorm(-l - b_grid)
  
  
  thresholds <- readMat('../Matlab/lookup_tables/thresholds.mat')
  cat('The adaptive soft threshold is\n')
  # Interpolate the soft-threshold estimate based on the estimated correlation coeff
  st.function <- splinefun(Sigma_UO_grid, thresholds$st.mat, method = "fmm", ties = mean)
  st <- st.function(abs(corr))
  cat(st, '\n')
  
  risk_function_st_adaptive <- Eb(st) + 1 / corr^2 - 1
  
  max_regret_st <- max(risk_function_st_adaptive / risk_oracle)
  cat('The adaptive soft-threshold has max regret\n', max_regret_st, '\n')
  
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
  cat('The pre-test has worst-case adaptation regret\n', max_regret_ttest, '\n')
  
  results <- list(
    max_regret_YU = max_regret_YU,
    max_regret_nonlinear = max_regret_nonlinear,
    max_regret_st = max_regret_st,
    max_regret_ttest = max_regret_ttest
  )
  
  return(results)
}
