plot_adaptive_and_minimax_priors <- function(YR,YU,VR,VU,VUR) {
  
  
  # YR is the restricted estimate
  # YU is the unbiased estimate
  # B_grid is the upper bound on the scaled bias e.g. [0.5 1 2];
  #global corr_sq %if passing in corr_str as a string
  # Calculate the over-id test statistic
  YO <- YR - YU
  VO <- VR - 2 * VUR + VU
  VUO <- (VUR - VU)
  tO <- YO / sqrt(VO)
  GMM <- YU - VUO / VO * YO; V_GMM <- VU - VUO/VO*VUO;
  corr <- VUO / sqrt(VO) / sqrt(VU)

  
  # Define variables
  Sigma_UO_grid <- abs(tanh(seq(-3, -0.05, 0.05)))
   
  B = 9;
  b_grid <- seq(-B, B, 0.05)# minimax?
  
  # Form worst-case prior
  prior_wtd <- readMat('../Matlab/lookup_tables/priors_wtd.mat');      
  prior <- readMat('../Matlab/lookup_tables/priors.mat');        
  
  b_grid <-prior_wtd$b.grid
  # spline interpolate priors for each b_grid based on observed
  Kb <- length(b_grid)
  prior_function_adaptive <- numeric(Kb)
  prior_adaptive <- numeric(Kb) # no weight
  
  for (i in 1:Kb) {
    prior.spline.function <- splinefun(Sigma_UO_grid, 
                                       prior_wtd$adaptive.prior.w[i,], method = "fmm", ties = mean)
    prior_function_adaptive[i] <- prior.spline.function(abs(corr)) # reweighted by minimax risk
    prior.spline.function <- splinefun(Sigma_UO_grid, 
                                       prior$x.mat[i,], method = "fmm", ties = mean)
    prior_adaptive[i] <- prior.spline.function(abs(corr))
  }
  #cbind(b_grid,prior_function_adaptive,prior_adaptive)[350:450,]
  
  # Create a data frame with x and y values
  data <- data.frame(x = b_grid, y = prior_function_adaptive)
  data$x <- round(data$x*5)/5
  data_sum <- aggregate(y ~ x, data, sum)
  
  minimax <- readMat('../Matlab/lookup_tables/priors_bounded_normal_mean.mat')   # Load priors but harmonize the b_grid
  index_low <- minimax$B.grid==1
  minimax_low <-  data.frame(x = minimax$b.max.grid, z = minimax$prior.bounded.normal.mean.mat[,index_low])
  minimax_low$x <- round(minimax_low$x*5)/5
  minimax_low_sum <- aggregate(z ~ x, minimax_low, sum)
  
  index_high <- minimax$B.grid == 5
  minimax_high <- data.frame(x = minimax$b.max.grid, z = minimax$prior.bounded.normal.mean.mat[, index_high])
  minimax_high$x <- round(minimax_high$x * 5) / 5
  minimax_high_sum <- aggregate(z ~ x, minimax_high, sum)
  
  data_sum$z <- minimax_low_sum$z
  data_sum$w <- minimax_high_sum$z
  
  data_sum$w[data_sum$w < 0.001] <- NA
  data_sum$y[data_sum$y < 0.001] <- NA
  data_sum$z[data_sum$z < 0.001] <- NA
  
  
  
  data_sum_long <- data_sum %>%
    pivot_longer(!x, names_to = "prior", values_to = "probability") %>%
    drop_na()
  
  data_sum_long$prior[data_sum_long$prior == "y"] <- "Adaptive"
  data_sum_long$prior[data_sum_long$prior == "z"] <- "Least favorable given B=1"
  data_sum_long$prior[data_sum_long$prior == "w"] <- "Least favorable given B=5"
  
  figurename <- paste("minimax_prior_sigmatb_", round(abs(corr) * 100) / 100, "_B", B, ".png", sep = "")
  par(mar = c(5, 5, 4, 5), pty = 'm',cex=1.2)  # Adjust margins as needed
  
  plot <- ggplot(data = data_sum_long,aes(x = x,y=probability, shape=prior, fill=prior )) +
    geom_bar(stat = "identity", position = position_identity(), alpha = 0.3) +
    geom_point(size = 3)  +
    labs(x = "b", y = "Prior") +
    theme_minimal() +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
    scale_x_continuous(breaks = seq(min(b_grid), max(b_grid), by = 2))  +
    theme(legend.position="bottom", legend.title=element_blank()) 
  
  print(plot)
  
   
  ggsave(filename = figurename,
         plot, width = 15, height = 10, units = "cm")
  
  dev.off()

}