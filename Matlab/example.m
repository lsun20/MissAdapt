%% Provide the restricted and the unrestricted estimates
% and the correlation between the restricted and unrestricted estimators. 
 

YR = [33.52740]; VR = [1.811236^2];
YU = [52.94609]; VU = [ 2.539776^2];
VUR = [1.811236^2];
YO = YR - YU;
VO = VR - 2*VUR + VU;
VUO = (VUR - VU);
disp('The over-id test statistic is')
tO = YO./sqrt(VO);
disp(tO)
disp('The efficient estimator is')
GMM = YU - VUO./VO .* YO;
disp(GMM)
disp('The correlation coefficient is')
corr = VUO./sqrt(VO)./sqrt(VU);
disp(corr)

%% Read in normalized estimates and scale them up for the application
% will write to sim_results/ a formatted table for each parameter
%%
[results] = adaptive_estimate(YR(1),YU(1),VR(1),VU(1),VUR(1),corr(1))
T = array2table(results)
T.Properties.VariableNames(1:7) = {'Y_U','Y_R',...
    'Y_O','GMM','Fully nonlinear','Adaptive soft-threshold',...
    'Pre-test'}
T.Properties.RowNames(1:4) = {'Estimate','Std Error','Max Regret','Threshold'}
writetable(T,'results.csv','WriteRowNames',true)
% plot risk functions
%adaptive_plot(YR(1),YU(1),VR(1),VU(1),VUR(1),9,corr(1))

%% Form locus of B-minimax estimates
%minimax_locus_plot(YR(1),YU(1),VR(1),VU(1),VUR(1),9,corr(1))
%minimax_prior_plot(YR(1),YU(1),VR(1),VU(1),VUR(1),9,corr(1))

 






 

