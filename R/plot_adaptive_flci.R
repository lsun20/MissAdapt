plot_adaptive_flci <- function(YR, YU, VR, VU, VUR) {
 
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
  #minimax <- readMat('../Matlab/lookup_tables/priors_bounded_normal_mean.mat')   # Load priors
  flci <- readMat("../Matlab/lookup_tables/flci_adaptive_cv.mat")
  flci_minimax <- readMat("../Matlab/lookup_tables/flci_minimax_cv.mat")
  
  Sigma_UO_grid <- abs(tanh(seq(-3, -0.05, 0.05)))
  
  
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
  
  # Calculate adaptive flci critical values
  KB <- length(flci$B.grid); B <- 9;
  cv_function_adaptive <- numeric(KB)

  for (i in 1:KB) {
    cv.spline.function <- splinefun(Sigma_UO_grid, flci$min.c.vec[i,], method = "fmm", ties = mean)
    cv_function_adaptive[i] <- cv.spline.function(abs(corr))
  }
  ubound.adaptive<-  adaptive_nonlinear + cv_function_adaptive*sqrt(VU)  #upper conf bound
  lbound.adaptive<- adaptive_nonlinear - cv_function_adaptive*sqrt(VU)  #lower conf bound
  
  # Calculate B-minimax flci critical values
  KB <- length(flci_minimax$B.minimax.grid); B <- 9;
  cv_minimax_function_adaptive <- numeric(KB)
  
  for (i in 1:KB) {
    cv.minimax.spline.function <- splinefun(Sigma_UO_grid, flci_minimax$min.c.vec[i,], method = "fmm", ties = mean)
    cv_minimax_function_adaptive[i] <- cv.minimax.spline.function(abs(corr))
  }
  
  # Plot the critical values
  # Create a data frame with the required values
  B_grid_scale <- flci$B.grid
  B_grid_scale[1] <- 0 #0.01 is approx 0 
  B_grid_scale <- as.vector(B_grid_scale)
  cv_minimax_function_adaptive <- c(1.96*sqrt(1-corr^2),t(cv_minimax_function_adaptive))
  df <- data.frame(B_grid_scale,  cv_function_adaptive,cv_minimax_function_adaptive)
  
  plot <- ggplot(df, aes(x = B_grid_scale )) +
    geom_line(aes( y = cv_function_adaptive, linetype = 'dotted'),
                color = 'black', size = 1) +
    geom_line(aes(y =  cv_minimax_function_adaptive, linetype = 'solid'),
                color = 'black', size = 1) +
    geom_hline(yintercept = 1.96, linetype = "solid", color = "red") +
    labs(x = "B", y = "FLCI critical values") +
    theme_minimal()  
  
  plot <- plot +  theme(legend.position = "bottom") +
    scale_linetype_manual(name = "CI type", values = c('dotted','solid'),
                          labels = c(TeX('Adaptive FLCI'), "Minimax FLCI"))   
  print(plot)
  figurename <-  paste("adaptive_minimax_flci_sigmatb_", round(abs(corr) * 100) / 100, "_B", B, ".png", sep = "")
  ggsave(filename = figurename,
         plot, width = 10, height = 10, units = "cm")
  
  # Plot the figure
  B_grid_scale <- sqrt(flci$B.grid)
  B_grid_scale[1] <- 0 #0.01 is approx 0 
  
  figurename <- paste("adaptive_flci_sigmatb_", round(abs(corr) * 100) / 100, "_B", B, ".png", sep = "")
  par(mar = c(5, 5, 4, 5), pty = 'm',cex=1.2)  # Adjust margins as needed
  
  # Create a data frame with the required values
  df <- data.frame(B_grid_scale, adaptive_nonlinear,ubound.adaptive,lbound.adaptive, GMM, YU, V_GMM, VU)
  
  # Plot using ggplot2
  plot <-  ggplot(df, aes(x = B_grid_scale)) +
  geom_hline(aes(yintercept = adaptive_nonlinear, linetype = 'adaptive'), color = 'black', size = 1) +
  geom_hline(aes(yintercept = GMM, linetype = 'GMM'), color = 'black', size = 1) +
  geom_hline(aes(yintercept = YU, linetype = 'YU'), color = 'black', size = 1) +
  geom_ribbon(aes(ymin = lbound.adaptive,
                    ymax = ubound.adaptive),  alpha=.15,
                color = 'black', linetype = 'solid', size = 1) +
    geom_hline(yintercept = GMM + 1.96 * sqrt(V_GMM), color = 'black', linetype = 'dashed', size = 1) +
    geom_hline(yintercept = GMM - 1.96 * sqrt(V_GMM), color = 'black', linetype = 'dashed', size = 1) +
    geom_hline(yintercept = YU + 1.96 * sqrt(VU), color = 'black', linetype = 'dotted', size = 1) +
    geom_hline(yintercept = YU - 1.96 * sqrt(VU), color = 'black', linetype = 'dotted', size = 1) +
    labs(x = 'B', y = '95% confidence intervals') +
    scale_x_continuous(expand = c(0, 0),breaks = c(0, 1, 2, 3), labels = c("0", "1", "4", "9")) +
    ylim((min(YU, YR)-1.96*sqrt(VU)) * 0.99, (max(YU, YR)+1.96*sqrt(VU)) * 1.01)+
    theme_minimal()  
    
  plot <- plot +  theme(legend.position = "right") +
    scale_linetype_manual(name = "Estimator", values = c("solid", "dashed", "dotted"),
                        labels = c("Adaptive", "GMM", "YU"))  
 
  ggsave(filename = figurename,
         plot, width = 15, height = 10, units = "cm")
  
  dev.off()
  message("Figure saved successfully as ",figurename)
}
