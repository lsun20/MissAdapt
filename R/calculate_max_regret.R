calculate_max_regret <- function(VR, VU, VUR) {
  VO <- VR - 2 * VUR + VU
  VUO <- (VUR - VU)
  corr <- VUO / sqrt(VO) / sqrt(VU)
  Sigma_UO_grid <- abs(tanh(seq(-3, -0.05, 0.05))) # grid of the correlation coefficient
  
  # Interpolate the risk function for the adaptive estimator
  if(!file.exists("../Matlab/lookup_tables/risk.mat")){stop("Check if ../Matlab/lookup_tables/risk.mat exists")}
  risk <-  readMat("../Matlab/lookup_tables/risk.mat")
  mse <- readMat('../Matlab/lookup_tables/emse_corr.mat')
  
  b_grid <- risk$b.grid; risk_mat<- risk$risk.mat;
  Kb <- length(b_grid)
  risk_function_adaptive <- numeric(Kb)
  
  for (i in 1:Kb) {
    risk.spline.function <- splinefun(Sigma_UO_grid, risk_mat[i,], method = "fmm", ties = mean)
    risk_function_adaptive[i] <- corr^2 * risk.spline.function(abs(corr))
  }
  
  # Calculate the oracle risk function
  if(!file.exists("../Matlab/lookup_tables/minimax_rho_B9.csv")){stop("Check if ../Matlab/lookup_tables/minimax_rho_B9.csv exists")}
  rho_tbl <- read.csv('../Matlab/lookup_tables/minimax_rho_B9.csv', header = FALSE)
  rho_b_over_sigma_function<- splinefun(rho_tbl[,1], rho_tbl[,2], method = "fmm", ties = mean)
  rho_b_over_sigma <- rho_b_over_sigma_function(abs(b_grid))
  risk_oracle <- corr^2 * rho_b_over_sigma + 1 - corr^2
  
  max_regret_YU = max(1/risk_oracle);
  message('The unrestricted estimator Y_U has max regret \n', round(max_regret_YU,2), '\n')
  
  max_regret_nonlinear <- max(risk_function_adaptive / risk_oracle)
  message('The adaptive estimator minimizes the max regret to be \n', round(max_regret_nonlinear,2), '\n')
  
  # Add soft threshold risk functions
  Eb <- function(l) 1 + l^2 + (b_grid^2 - 1 - l^2) * (pnorm(l - b_grid) - pnorm(-l - b_grid)) + (-b_grid - l) * dnorm(l - b_grid) - (l - b_grid) * dnorm(-l - b_grid)
  
  if(!file.exists("../Matlab/lookup_tables/thresholds.mat")){stop("Check if ../Matlab/lookup_tables/thresholds.mat exists")}
  thresholds <- readMat('../Matlab/lookup_tables/thresholds.mat')
  # Interpolate the soft-threshold estimate based on the estimated correlation coeff
  st.function <- splinefun(Sigma_UO_grid, thresholds$st.mat, method = "fmm", ties = mean)
  st <- st.function(abs(corr))
  #message('The adaptive soft threshold is \n', round(st round(, '\n')
  
  risk_function_st_adaptive <- corr^2 * Eb(st) + 1 - corr^2 
  
  max_regret_st <- max(risk_function_st_adaptive / risk_oracle)
  message('The adaptive soft-threshold estimator has max regret\n', round(max_regret_st,2), '\n')
  
 
  # Use simulation to calculate the risk function for the pre-test estimator that switches between Y_U and Y_R
  set.seed(1)
  sims <- 100000
  x <- rnorm(sims)
  x_b <- matrix(x, sims, Kb) + matrix(b_grid, sims, Kb, byrow = TRUE)
  
  Ebsims_ht <- function(l) colSums(((x_b > l) * x_b + (x_b < l & x_b > -l) * (1 + VO / VUO) * x_b
                                    + (x_b < -l) * x_b - matrix(b_grid, sims, Kb, byrow = TRUE))^2) / sims
  
  risk_function_ht_ttest <- corr^2 * Ebsims_ht(1.96) + 1 - corr^2
  # Use analytical formula to calculate the risk function for the pre-test estimator that switches between Y_U and GMM
  Eb_ht <- function(l) {
    result <- 1+(b_grid^2-1) *(pnorm(l-b_grid)-pnorm(-l-b_grid))+ 
      (l-b_grid) *dnorm(l-b_grid) - (-l-b_grid) *dnorm(-l-b_grid)
    return(result)
  }

  risk_function_ht_ttest <-   corr^2 *Eb_ht(1.96) + 1 - corr^2
  
                         
  # Calculate max regret for various estimators
  max_regret_ttest <- max(risk_function_ht_ttest / risk_oracle)
  message('The pre-test estimator has max regret\n', round(max_regret_ttest,2), '\n')
  
  # Interpolate the hard-threshold estimate based on the estimated correlation coeff
  ht.function <- splinefun(Sigma_UO_grid, thresholds$ht.mat, method = "fmm", ties = mean)
  ht <- ht.function(abs(corr))
   
  risk_function_ht_adaptive <- corr^2 * Eb_ht(ht) + 1 - corr^2 
  
  max_regret_ht <- max(risk_function_ht_adaptive / risk_oracle)
  message('The adaptive hard-threshold estimator has max regret\n', round(max_regret_ht,2), '\n')
  # Calculate the adaptive ERM lambda
  lambda.function <-  splinefun(Sigma_UO_grid, mse$MSE.lambda.mat, method = "fmm", ties = mean)
  lambda <- lambda.function(abs(corr))
  
  Ebsims_MSE = function(l)  colSums((x_b^3 / (x_b^2 + l) - matrix(b_grid, sims, Kb, byrow = TRUE))^2) / sims
  # Calculate max regret for various estimators
  risk_function_erm <- corr^2 * Ebsims_MSE(1) + 1 - corr^2 
  risk_function_adaptive_erm <- corr^2 * Ebsims_MSE(lambda) + 1 - corr^2 
  max_regret_erm <- max(risk_function_erm / risk_oracle)
  max_regret_adaptive_erm <- max(risk_function_adaptive_erm / risk_oracle)
  
  results <- list(
    max_regret_YU = max_regret_YU,
    max_regret_nonlinear = max_regret_nonlinear,
    max_regret_st = max_regret_st,
    st = st,
    max_regret_ht = max_regret_ht,
    ht = ht,
    max_regret_ttest = max_regret_ttest,
    max_regret_erm = max_regret_erm,
    max_regret_adaptive_erm = max_regret_adaptive_erm,

    lambda = lambda,
    max_risk_adaptive = max(risk_function_adaptive),
    max_risk_st_adaptive = max(risk_function_st_adaptive),
    max_risk_ht_adaptive = max(risk_function_ht_adaptive),
    max_risk_ht_ttest = max(risk_function_ht_ttest),
    max_risk_erm = max(risk_function_erm),
    max_risk_adaptive_erm = max(risk_function_adaptive_erm)
  )
  
  return(results)
}
