plot_adaptive_and_minimax_estimates <- function(YR, YU, VR, VU, VUR) {
 
  # Calculate the over-id test statistic
  YO <- YR - YU
  VO <- VR - 2 * VUR + VU
  VUO <- (VUR - VU)
  tO <- YO / sqrt(VO)
  GMM <- YU - VUO / VO * YO; V_GMM <- VU - VUO/VO*VUO;
  corr <- VUO / sqrt(VO) / sqrt(VU)
  
  # Load lookup tables
  policy <- readMat("../Matlab/lookup_tables/policy.mat")
  risk <- readMat("../Matlab/lookup_tables/risk.mat")
  minimax <- readMat('../Matlab/lookup_tables/priors_bounded_normal_mean.mat')   # Load priors
  
  B = 9;
  # Define variables
  Sigma_UO_grid <- abs(tanh(seq(-3, -0.05, 0.05)))
  b_grid <- seq(-B, B, 0.05)
  y_grid <- c(tO)
  Ky <- length(y_grid)
  Kb <- length(b_grid)
  y_grid_ub <- y_grid + 0.05
  y_grid_lb <- y_grid - 0.05
  Pi <- matrix(0, nrow = Ky, ncol = Kb)
  
  # Calculate Pi
  for (i in 1:Kb) {
    Pi[,i] <- pnorm((y_grid_ub  + y_grid ) / 2, mean = b_grid[i], sd = 1) -
      pnorm((y_grid_lb  + y_grid ) / 2, mean = b_grid[i], sd = 1)
  }
  
  # Define Psi function
  Psi <- function(mu_grid) {
    (Pi %*% (mu_grid * (b_grid))) / (Pi %*% mu_grid)
  }
  

  t_minimax <- numeric(length(minimax$B.grid))
  
  # Calculate t_minimax
  for (i in 1:length(minimax$B.grid)) {
    t_minimax[i] <- Psi(minimax$prior.bounded.normal.mean.mat[, i])  # Minimax estimator
  }
  
  # Calculate Y_minimax
  Y_minimax <- (VUO / sqrt(VO)) * t_minimax + GMM
  B_grid <- c(0, minimax$B.grid)
  Y_minimax <- c(GMM, Y_minimax)
  
  # Calculate adaptive estimates
  Ky <- length(policy$y.grid)
  psi.grid <- rep(Ky, 0)
  for (i in seq(1, Ky, 1)) {
    psi.function <- splinefun(Sigma_UO_grid, policy$psi.mat[i,], method = "fmm", ties = mean)
    psi.grid[i] <- psi.function(abs(corr))
  }
  
  psi.grid.extrap <- splinefun(policy$y.grid, psi.grid, method = "natural")
  t_tilde <- psi.grid.extrap(tO)
  adaptive_nonlinear <- VUO / sqrt(VO) * t_tilde + GMM
  
  # Calculate adaptive risk
  Kb <- length(B_grid)
  risk_function_adaptive <- numeric(Kb)

  for (i in 1:Kb) {
    risk.spline.function <- splinefun(Sigma_UO_grid, risk$risk.mat[abs(risk$b.grid - B_grid[i])<0.001,], method = "fmm", ties = mean)
    risk_function_adaptive[i] <- corr^2 * risk.spline.function(abs(corr))
  }
   
  # Calculate oracle risk
  # Calculate the oracle risk function
  if(!file.exists("../Matlab/lookup_tables/minimax_rho_B9.csv")){stop("Check if ../Matlab/lookup_tables/minimax_rho_B9.csv exists")}
  rho_tbl <- read.csv('../Matlab/lookup_tables/minimax_rho_B9.csv', header = FALSE)
  rho_b_over_sigma_function<- splinefun(rho_tbl[,1], rho_tbl[,2], method = "fmm", ties = mean)
  rho_b_over_sigma <- rho_b_over_sigma_function(abs(B_grid))
  #rho_b_over_sigma <- interp1(rho_tbl[, 1], rho_tbl[, 2], abs(B_grid), method = "spline")
  risk_oracle <- corr^2 * rho_b_over_sigma + 1 - corr^2
  
  # Plot the figure
  B_grid_scale <- sqrt(B_grid)
 
  # Save the figure
  figurename <- paste("minimax_locus_sigmatb_", round(abs(corr) * 100) / 100, "_B", B, ".png", sep = "")
  png(figurename, width = 640, height = 420, units = "px")
  # Create the plot area
  par(mar = c(5, 5, 4, 5), pty = 'm',cex=1.2)  # Adjust margins as needed
  plot(B_grid_scale, Y_minimax, type = "n", las = 1,  
       xlab = "B/√Σ_O", xaxt = "n", ylab = "Estimate", ylim = c(min(YU, YR) * 0.99, max(YU, YR) * 1.01))
   
  axis(1, at = c(0, 1, 2, 3), labels = c("0", "1", "4", "9"),las=1)
  # Add triangles at specified coordinates
  points(B_grid_scale, Y_minimax, pch = 2, col = "black", cex = 1)
  
  abline(h=YR, col = "black", pch = 2); text(x=0, y=YR*1.0025, 'YR', adj = 0)
  abline(h=GMM, col = "black", pch = 4); text(x=0, y=GMM*1.0025, 'GMM', adj = 0)
  abline(h=adaptive_nonlinear, col ="black", pch = 3); text(x=0, y=adaptive_nonlinear*1.0025, 'Adaptive', adj = 0)
  abline(h=YU, col = "black", pch = 1); text(x=0, y=YU*1.0025, 'YU', adj = 0)
  
  par(new = TRUE) 
  plot(B_grid_scale, risk_oracle, type = "n", xlab = "", ylab = "", 
       axes = FALSE,ylim = c(min(risk_function_adaptive) - 0.2, max(risk_function_adaptive) + 0.2))
  axis(4, at = pretty(range(risk_oracle,risk_function_adaptive)), col.axis = "blue", las = 1)
  mtext('Mean squared error relative to YU', side=4, line=3,cex=1.2)
  
  points(B_grid_scale, risk_oracle, type = "l", col = "blue", lwd = 2)
  points(B_grid_scale, risk_function_adaptive, type = "l", col = "blue", lty = 3, lwd = 2 )
  abline(h=1,lty=2 )
  
  par(new = FALSE)
  legend("topright", legend = c("B-minimax estimates", "Oracle risk", "Adaptive risk"), 
         col = c("black", "blue", "blue"), pch = c(2, NA, NA), lty = c(NA,  1, 3), lwd = c(1,2,2),bty = "n")
 
  dev.off()
  message("Figure saved successfully as ",figurename)
}
