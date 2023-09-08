% DESCRIPTION =========================================================
%
% This script calculates adaptive estimators as 
% described in Armstrong, Kline, and Sun (2023) for scalar parameters. 
% The Vignette offers a detailed step-by-step guide for this script. 
% We encourage you to explore it. If you encounter any issues, please 
% don't hesitate to reach out.
%
% First, this script calls the function "adaptive_estimates.m" 
% to construct the adaptive estimator based on robust and restricted 
% estimators. The same function also computes the worst-case adaptation 
% regret.
%
% Second, this script invokes the function 
% "minimax_locus_plot.m" to generate plots illustrating 
% the locus of B-minimax estimates and associated risk functions.
%
% Ensure that the required lookup tables are stored in the correct 
% directory "./lookup_tables/XXX," where XXX refers to various 
% Matlab objects provided by the authors. These objects contain 
% pre-tabulated solutions for the adaptive problem across a fine grid. 
% The functions in this script rely on these inputs to interpolate the 
% pre-tabulated solutions for calculating the adaptive estimator
% for specific applications.
%

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
 
%% Constrain the worst-case risk to be no more than 20% of Y_U
[results] = const_estimate(YR(1),YU(1),VR(1),VU(1),VUR(1),corr(1))

T = array2table(results)
T.Properties.VariableNames(1:2) = {'Fully nonlinear','Adaptive soft-threshold'}
T.Properties.RowNames(1:4) = {'Estimate','Max Regret','Max Risk','Threshold'}
writetable(T,'const_results.csv','WriteRowNames',true)

%% Form locus of B-minimax estimates
minimax_locus_plot(YR(1),YU(1),VR(1),VU(1),VUR(1),corr(1))
 



 






 

