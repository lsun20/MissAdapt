plot_shrinkage_estimates <- function(YR, YU, VR, VU, VUR) {
 
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
  thresholds <- readMat('../Matlab/lookup_tables/thresholds.mat')
  
  B = 9;
  # Define variables
  Sigma_UO_grid <- abs(tanh(seq(-3, -0.05, 0.05)))
  b_grid <- seq(-B, B, 0.05)
  y_grid <- c(tO)
  Ky <- length(y_grid)
  Kb <- length(b_grid)
  
  
  # Calculate adaptive estimates
  Ky <- length(policy$y.grid)
  psi.grid <- rep(Ky, 0)
  for (i in seq(1, Ky, 1)) {
    psi.function <- splinefun(Sigma_UO_grid, policy$psi.mat[i,], method = "fmm", ties = mean)
    psi.grid[i] <- psi.function(abs(corr))
  }
  
  # Calculate adaptive soft-threshold
  # Interpolate the soft-threshold estimate based on the estimated correlation coeff
  st.function <- splinefun(Sigma_UO_grid, thresholds$st.mat, method = "fmm", ties = mean)
  st <- st.function(abs(corr))
  psi.st = ((policy$y.grid > st)*(policy$y.grid - st) +
              (policy$y.grid < -st)*(policy$y.grid + st)) + 0;
  st_str <- round(st,2);
  
  # Calculate pre-test
  psi.pretest <- ((policy$y.grid> 1.96)*(policy$y.grid) + 
                    (policy$y.grid < - 1.96)*(policy$y.grid)) + 0;
  # Calculate empirical risk minimizer
  psi.dCDH <-(policy$y.grid)^2/(policy$y.grid^2+1)*policy$y.grid
  
  figurename <- paste("policy_plot_sigmatb_", round(abs(corr) * 100) / 100, "_B", B, ".png", sep = "")
  
  # Create a data frame
  y_grid_scale <- policy$y.grid
  a <- psi.grid
  c <- psi.st
  b <- psi.pretest
  d <- psi.dCDH
  
  df <- data.frame(
    y_grid_scale,a,b,c,d
  )
  
  df_long <- tidyr::gather(df, key = "series", value = "y", -y_grid_scale)
  
  shrinkage_linetypes <- c('dashed', 'blank', 'solid' ,'dotdash')
  plot <- ggplot(filter(df_long, series != "b"), aes(x = y_grid_scale, y = y, linetype = series, color = series)) +
    geom_line(size = 1) +
    geom_point(data = filter(df_long, series == "b"), size = 1, shape = 1) +
    geom_vline(xintercept = c(st, -st, 1.96, -1.96), linetype = "solid", color = "grey") +
    geom_segment(aes(x = -3, xend = 3, y = -3, yend = 3), linetype = "dashed", color = "grey")+
    labs(
      x = TeX("$T_O$"),
      y = TeX("Estimate for scaled bias b/$\\sigma_O$"),
      title = NULL
    ) +
    
    theme_minimal() +
    
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.text = element_text(size = 14),
      legend.title = element_blank()
    ) +
    
    scale_x_continuous(
      breaks = c(-3, -2, -1, 0, st, 1.96, 3),
      limits = c(-3, 3),
      labels = c('-3', '-2', '-1', '0', st_str, '1.96', '3')
    ) +
    
    scale_y_continuous(
      breaks = c(-3, -2, -1, 0, 1, 2, 3),
      limits = c(-3, 3),
      labels = c('-3', '-2', '-1', '0', '1', '2', '3')
    ) +
    scale_color_manual(values = cbp1, name = "",
                       labels = c(TeX('Adaptive'),'Pre-test',
                                  TeX('Soft-threshold'),'ERM')) +
    scale_linetype_manual(values = shrinkage_linetypes, name = "",
                          labels = c(TeX('Adaptive'),'Pre-test',
                                     TeX('Soft-threshold'),'ERM')) +
    guides(color = guide_legend(ncol = 3))  # Adjust the number of columns in the legend
  print(plot)
  ggsave(filename = figurename,
         plot, width = 10, height = 10, units = "cm")
  dev.off()
  message("Figure saved successfully as ",figurename)
  
  figurename <- paste("weight_plot_sigmatb_", round(abs(corr) * 100) / 100, "_B", B, ".png", sep = "")
  
  df <- data.frame(
    y_grid_scale,a_w=a/y_grid_scale,b_w=b/y_grid_scale,c_w=c/y_grid_scale,d_w=d/y_grid_scale
  )
  df  <- filter(df, y_grid_scale != 0)
  df_long <- tidyr::gather(df, key = "series", value = "y", -y_grid_scale)
  
  shrinkage_linetypes <- c('dashed', 'blank', 'solid' ,'dotdash')
  plot <- ggplot(filter(df_long, series != "b_w"), aes(x = y_grid_scale, y = y, linetype = series, color = series)) +
    geom_line(size = 1) +
    geom_point(data = filter(df_long, series == "b_w"), size = 1, shape = 1) +
    geom_vline(xintercept = c(st, -st, 1.96, -1.96), linetype = "solid", color = "grey") +
    labs(
      x = TeX("$T_O$"),
      y = TeX('Weight placed on $Y_U$'),
      title = NULL
    ) +
    
    theme_minimal() +
    
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.text = element_text(size = 14),
      legend.title = element_blank()
    ) +
    
    scale_x_continuous(
      breaks = c(-3, -2, -1, 0, st, 1.96, 3),
      limits = c(-3, 3),
      labels = c('-3', '-2', '-1', '0', st_str, '1.96', '3')
    ) +
    scale_color_manual(values = cbp1, name = "",
                       labels = c(TeX('Adaptive'),'Pre-test',
                                  TeX('Soft-threshold'),'ERM')) +
    scale_linetype_manual(values = shrinkage_linetypes, name = "",
                          labels = c(TeX('Adaptive'),'Pre-test',
                                     TeX('Soft-threshold'),'ERM')) +
    guides(color = guide_legend(ncol = 3))  # Adjust the number of columns in the legend
  print(plot)
  ggsave(filename = figurename,
         plot, width = 10, height = 10, units = "cm")
  dev.off()
  message("Figure saved successfully as ",figurename)
  
}
