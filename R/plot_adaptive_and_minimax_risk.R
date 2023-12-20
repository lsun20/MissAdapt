plot_adaptive_and_minimax_risk <- function(YR, YU, VR, VU, VUR) {
 
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
  minimax <- readMat('../Matlab/lookup_tables/risk_bounded_normal_mean.mat')   # Load minimax risk
  B_minimax_grid = seq(0.1,9,0.1);
  b_grid <- minimax$b.max.grid;

  B = 9;
  # Define variables
  Sigma_UO_grid <- abs(tanh(seq(-3, -0.05, 0.05)))
  
  # Calculate risk of YR
  risk_function_YR <- (VO* b_grid^2 + VR)/VU;
  
  # Calculate adaptive risk
  Kb <- length(b_grid)
  risk_function_adaptive <- numeric(Kb)

  for (i in 1:Kb) {
    risk.spline.function <- splinefun(Sigma_UO_grid, risk$risk.mat[i,], method = "fmm", ties = mean)
    risk_function_adaptive[i] <- corr^2 * risk.spline.function(abs(corr))
  }
   
  # Calculate oracle risk
  # Calculate the oracle risk function
  if(!file.exists("../Matlab/lookup_tables/minimax_rho_B9.csv")){stop("Check if ../Matlab/lookup_tables/minimax_rho_B9.csv exists")}
  rho_tbl <- read.csv('../Matlab/lookup_tables/minimax_rho_B9.csv', header = FALSE)
  rho_b_over_sigma_function<- splinefun(rho_tbl[,1], rho_tbl[,2], method = "fmm", ties = mean)
  rho_b_over_sigma <- rho_b_over_sigma_function(abs(b_grid))
  #rho_b_over_sigma <- interp1(rho_tbl[, 1], rho_tbl[, 2], abs(B_grid), method = "spline")
  risk_oracle <- corr^2 * rho_b_over_sigma + 1 - corr^2
  
  # Plot the figure
  b_grid_scale <- sign(b_grid)*sqrt(abs(b_grid));
  

  # Assuming b_grid_scale, risk_function_YR, VU, b_max_grid_scale, corr_str, risk_function_minimax,
  # B_minimax_grid, risk_function_oracle are defined appropriately in your R environment.
  library(ggplot2)
  
  # Assuming b_grid_scale, risk_function_YR, VU, b_max_grid_scale, corr_str, risk_function_minimax,
  # B_minimax_grid, risk_function_oracle are defined appropriately in your R environment.
  # Define color and linetype preferences
  cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442")
  my_linetypes <- c('solid', 'dotted', 'dashed', 'dotdash','longdash')
  library(latex2exp)
  
  # Create a data frame
  df <- data.frame(
    b_grid_scale = b_grid_scale,
    a = rep(1, Kb),
    b = risk_function_YR,
    c =  corr^2 * minimax$risk.function.minimax[, B_minimax_grid == 1] + 1 - corr^2,
    d =  corr^2 * minimax$risk.function.minimax[, B_minimax_grid == 4] + 1 - corr^2,
    e = risk_oracle
  )
  
  # Melt the data frame to long format
  df_long <- tidyr::gather(df, key = "series", value = "y", -b_grid_scale)
  
  figurename <- paste("risk_plot_sigmatb_", round(abs(corr) * 100) / 100, "_B", B, ".png", sep = "")
  par(mar = c(5, 5, 4, 5), pty = 'm',cex=1.2)  # Adjust margins as needed
  
  # Plotting using ggplot2
  plot <- ggplot(df_long, aes(x = b_grid_scale, y = y, linetype = series, color = series)) +
    geom_line(size=1.2) +
    labs(x = TeX("b/SD($Y_R$-$Y_U$)"), y = TeX("MSE Relative to $\\sigma^2_U$"), title = NULL) +
    scale_x_continuous(breaks = c(-3, -2, -1, 0, 1, 2, 3), labels = c('-9', '-4', '-1', '0', '1', '4', '9')) +
    ylim(c(min(risk_oracle-0.1), 1.5)) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.text = element_text(size = 12),  # Adjust legend text size
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 14)
    ) +
    scale_color_manual(values = cbp1, name = "",
                       labels = c(TeX('$Y_U$'),TeX('$Y_R$'),
                       TeX('Minimax under |b/SD($Y_R$-$Y_U$)|$\\leq$1'),
                       TeX('Minimax under |b/SD($Y_R$-$Y_U$)|$\\leq$4'),'Oracle')) +
    scale_linetype_manual(values = my_linetypes, name = "",
                          labels = c(TeX('$Y_U$'),TeX('$Y_R$'),
                        TeX('Minimax under |b/SD($Y_R$-$Y_U$)|$\\leq$1'),
                        TeX('Minimax under |b/SD($Y_R$-$Y_U$)|$\\leq$4'),'Oracle')) +
    guides(color = guide_legend(ncol = 3))
  print(plot)
  ggsave(filename = figurename,
         plot, width = 15, height = 10, units = "cm")
  
  
  dev.off()
  message("Figure saved successfully as ",figurename)
}
