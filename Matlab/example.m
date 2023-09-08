%% Provide the restricted and the unrestricted estimates
% and the correlation between the restricted and unrestricted estimators. 
 

YR = [2408.8413]; VR = [220.6057^2];
YU = [2217.2312]; VU = [257.3961^2];
VUR = [211.4772^2];

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
% will write the results to a formatted table for each parameter
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
minimax_locus_plot(YR(1),YU(1),VR(1),VU(1),VUR(1),corr(1))
%minimax_prior_plot(YR(1),YU(1),VR(1),VU(1),VUR(1),9,corr(1))

%% Constrain the worst-case risk to be no more than 20% of Y_U
[results] = const_estimate(YR(1),YU(1),VR(1),VU(1),VUR(1),corr(1))

T = array2table(results)
T.Properties.VariableNames(1:2) = {'Fully nonlinear','Adaptive soft-threshold'}
T.Properties.RowNames(1:4) = {'Estimate','Max Regret','Max Risk','Threshold'}
writetable(T,'const_results.csv','WriteRowNames',true)


 






 

