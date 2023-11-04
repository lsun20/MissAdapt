function minimax_prior_plot(YR,YU,VR,VU,VUR)
% YR is the restricted estimate
% YU is the unbiased estimate
% B_grid is the upper bound on the scaled bias e.g. [0.5 1 2];
global corr_sq %if passing in corr_str as a string
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



# Calculate adaptive risk
Kb <- length(b_grid)
prior_function_adaptive <- numeric(Kb)
prior_adaptive <- numeric(Kb) # no weight

for (i in 1:Kb) {
  prior.spline.function <- splinefun(Sigma_UO_grid, 
                                     prior_wtd$adaptive.prior.w[i,], method = "fmm", ties = mean)
  prior_function_adaptive[i] <- prior.spline.function(abs(corr))
  prior.spline.function <- splinefun(Sigma_UO_grid, 
                                     prior$x.mat[i,], method = "fmm", ties = mean)
  prior_adaptive[i] <- prior.spline.function(abs(corr))
}
cbind(b_grid,prior_function_adaptive,prior_adaptive)[350:450,]

# Create a data frame with x and y values
data <- data.frame(x = b_grid, y = prior_function_adaptive)
data$x <- round(data$x*5)/5
data_sum <- aggregate(y ~ x, data, sum)

minimax <- readMat('../Matlab/lookup_tables/priors_bounded_normal_mean.mat')   # Load priors
index1 <- minimax$B.grid==1
minimax1 <-  data.frame(x = minimax$b.max.grid, z = minimax$prior.bounded.normal.mean.mat[,index1])
minimax1$x <- floor(minimax1$x*5)/5
minimax1_sum <- aggregate(z ~ x, minimax1, sum)

index5 <- minimax$B.grid == 5
minimax5 <- data.frame(x = minimax$b.max.grid, z = minimax$prior.bounded.normal.mean.mat[, index5])
minimax5$x <- floor(minimax5$x * 5) / 5
minimax5_sum <- aggregate(z ~ x, minimax5, sum)

data_sum$z <- minimax1_sum$z
data_sum$w <- minimax5_sum$z


plot <- ggplot(data_sum, aes(x = x,y=y )) +
  labs(x = "b", y = "Probability") +
  theme_minimal() 
# Add points for y
plot <- plot + geom_point(data = subset(data_sum, y > 0.001),aes(y = y), shape = 1) +
  geom_bar(aes(y = y),stat = "identity", fill = "blue", alpha = 0.3)

# Add bars for z

plot <- plot + geom_point(data = subset(minimax1_sum, z > 0.001), 
                          aes(y = z), shape = 2) +
  geom_bar(data = minimax1_sum, aes(y = z), stat = "identity", fill = "green", alpha = 0.3)
# Add bars for w

plot <- plot + geom_point(data = subset(minimax5_sum, z > 0.001), 
                          aes(y = z), shape = 3) +
  geom_bar(data = minimax5_sum, aes(y = z), stat = "identity", fill = "red", alpha = 0.3)


# Remove vertical grid lines
plot <- plot + theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
  scale_x_continuous(breaks = seq(min(b_grid), max(b_grid), by = 2)) +
  scale_colour_manual(name = "Treatment & State",
                      labels = c("Control, Non-F", "Control, Flwr", "Exclosure, Non-F"),
                      values = c("blue", "red", "green")) +   
  scale_shape_manual(name = "Treatment & State",
                     labels = c("Control, Non-F", "Control, Flwr", "Exclosure, Non-F"),
                     values = c(1,2,3))

print(plot)

###



load('sim_results/minimax_bounded_normal_mean.mat');


%% plot adaptive prior
clf(figure(1))
fig = figure(1)
b_grid_scale = b_grid;
b_max_grid_scale =b_max_grid;
yyaxis left
p1=plot(b_grid_scale,prior_function_adaptive,'-',...
        'DisplayName','Adaptive','LineJoin','miter','LineWidth',2 ); ylim([0 0.6])


xlabel('b/\Sigma_O^{1/2}');
%             xticks([-9 -4  -1  0 1 4 9]);
xlim([-9 9]);
ax=gca; ax.FontSize = 14;
%             ylabel('Adaptive prior','FontSize',14 );  
hold on

yyaxis right
p2=plot(b_max_grid_scale,prior_bounded_normal_mean_mat(:,B_grid == 1),':',...
        'DisplayName','B/\Sigma_O^{1/2}=1','LineJoin','miter','LineWidth',2 );
hold on

p3=plot(b_max_grid_scale,prior_bounded_normal_mean_mat(:,B_grid == 5),'--',...
        'DisplayName','B/\Sigma_O^{1/2}=5','LineJoin','miter','LineWidth',2 );
ylim([0 0.6])
%             ylabel('Worst-case prior for B-minimax','FontSize',14,'rotation',-90,'VerticalAlignment','bottom'); ylim([0 0.6])
yticks([])
hold off

legend([p1 p2 p3]...
       ,'Location','southoutside','Orientation','horizontal',...
       'FontSize',14 ,'NumColumns', 3);
legend('boxoff');
%% Save figures
figurename = strcat( 'sim_results/minimax_prior_sigmatb_',...
                     string(round(abs(corr(1))*100)/100),'_B',string(B),'.png');
set(fig,'Units','Inches');
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fig,figurename,'-dpng','-r0')

end