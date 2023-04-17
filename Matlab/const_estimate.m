function [results_scale_minimax] = ...
    const_estimate(YR,YU,VR,VU,VUR,corr_str,B_grid)
% YR is the restricted estimate
% YU is the unbiased estimate
% corr_str can be a string  e.g. 0.52 = '052' that is used to look up the
% csv
% corr_str can also be the numerical corr coef that is used to interpolate
% the estimates
% B_grid is the (optional) upper bound on the scaled bias e.g. [0.5 1 2];
if nargin < 7 || isempty(B_grid)
     B_grid = 9;
end
global corr_sq %if passing in corr_str as a string
rho_tbl = readmatrix('lookup_tables/minimax_rho_B9.csv');

    %% Over-id test
    YO = YR - YU;
    VO = VR - 2*VUR + VU;
    VUO = (VUR - VU);
    disp('The over-id test statistic is')
    tO = YO/sqrt(VO);
    disp(tO)
    disp('The efficient estimator is')
    CUE = YU - VUO/VO * YO;
    disp(CUE)
    disp('The efficient estimator has variance')
    V_CUE = VU - VUO^2/VO;
    disp(V_CUE)
    disp('The correlation coefficient is')
    corr = VUO/sqrt(VO)/sqrt(VU);
    disp(corr)
    disp('The worst-case regret of Y_U is')
    disp(VU / V_CUE)
    
    %% Loop over upper bounds - if interpolate the corr_str
    if isnumeric(corr_str)
        Sigma_UO_grid = abs(tanh((-3:0.05:-0.05)));  % take out zero correlation coeff
        if abs(corr_str) > max(Sigma_UO_grid) || abs(corr_str) < min(Sigma_UO_grid)
            error("Extrapolating beyond the grid! " + ...
                "Pre-tabulations are available for absolute values of " + ...
                "correlation coefficient between " + ...
            min(Sigma_UO_grid) + " and " + max(Sigma_UO_grid) + ". " + ...
            "Optimize numerically for the given rho value.")
        end

        for i = 1:length(B_grid)
            B = B_grid(i);
            B
       

            %% Form nonlinear adaptive estimates
            % Looks for the nonlinear estimates stored in the /lookup_tables/
            % subdirectory
            
            load('lookup_tables/const_policy.mat');      %1=20%, 2=15%, 3= 5%    
            % spline interpolate nonlinear estimate and apply scaling
            [X_grid,Y_grid] = meshgrid(Sigma_UO_grid,y_grid);
            disp('The adaptive estimate is')   
            t_tilde = interp2(X_grid,Y_grid,psi1_mat,abs(corr_str),tO,'spline');
            adaptive_nonlinear = VUO/sqrt(VO) * t_tilde + CUE;
            % get the gradient and calculate the SURE
            Ky = length(y_grid); psi = zeros(Ky,1);
            for i = 1:Ky
               psi(i) = interp1(Sigma_UO_grid,psi1_mat(i,:),abs(corr_str),'spline');
            end
            SURE = VUO^2/VO*((t_tilde - tO)^2 + 2*interp1(y_grid,gradient(psi,y_grid),tO,'spline') - 1) + ...
                VU - VUO^2/VO;
            disp(interp1(y_grid,gradient(psi,y_grid),tO,'spline'))
            
         
            %% Risk function
            load('lookup_tables/const_risk.mat');
            Kb = length(b_grid); risk_function_adaptive = zeros(Kb,1); 
            for i = 1:Kb
               risk_function_adaptive(i) = interp1(Sigma_UO_grid,risk1_mat(i,:),abs(corr_str),'spline');
            end
            % calculate the oracle risk function
            rho_b_over_sigma = interp1(rho_tbl(:,1),rho_tbl(:,2),abs(b_grid),'spline');
            risk_oracle = rho_b_over_sigma + 1/corr_str^2 -1;              
            
            penalty_nonlinear = max(risk_function_adaptive./risk_oracle);
             % report adapatation penalty
            disp('The adaptive penalty estimate is')
            disp(penalty_nonlinear)

            %% Form thresholded estimates (minimax loss)
            load('lookup_tables/const_thresholds.mat');
            disp('The adaptive soft threshold is')
            st = interp1(Sigma_UO_grid,st1_mat,abs(corr_str),'spline');
             disp(st)
            disp('The adaptive soft-thresholded estimate is')
            adaptive_st = VUO/sqrt(VO) * ((tO > st)*(tO - st) + (tO < -st)*(tO + st)) + CUE;
            disp(adaptive_st)
            SURE_st = VUO^2/VO*( (abs(tO) > st)*(st^2+2) + (abs(tO) < st)*(tO^2) - 1) + ...
                VU - VUO^2/VO; 
            Kb = length(b_grid); 
            Eb = @(l) 1+l^2+...
             (b_grid.^2-1-l^2).*(normcdf(l-b_grid)-normcdf(-l-b_grid))+...
             (-b_grid-l).*normpdf(l-b_grid) - (l-b_grid).*normpdf(-l-b_grid);
            risk_function_st_adaptive = Eb(st) + 1/corr_str^2 - 1;

            %% Calculate penalty for various estimators
            disp('The adaptive soft-threshold has penalty')
            disp(max(risk_function_st_adaptive./risk_oracle))
            
            %% Calculate the wrost-case risk increase for various estimators
            disp('The adaptive estimator has worst case risk')
            max_risk_adaptive = corr_str^2*max(risk_function_adaptive )  
            disp('The adaptive soft-threshold has worst case risk')
            max_risk_st_adaptive = corr_str^2*max(risk_function_st_adaptive)
            %% put everything in a matrix
            results = zeros(4,2); 
            results(1,1)=adaptive_nonlinear; results(2,1)=penalty_nonlinear;
            results(3,1)=max_risk_adaptive ;
            results(1,2)=adaptive_st; results(2,2)=max(risk_function_st_adaptive./risk_oracle); 
            results(3,2)=max_risk_st_adaptive; results(4,2)= st;
            results_scale_minimax = results;

        end
    else
    %% Loop over upper bounds - if look up the corr_str file
        for i = 1:length(B_grid)
            B = B_grid(i);
            B
            %% Form nonlinear adaptive estimates
            % Looks for the nonlinear estimates stored in the /lookup_tables/
            % subdirectory


            Tbl = readtable(strcat('lookup_tables/const_adaptive_psi_sigmatb_',corr_str,'_B',string(B),'.csv'));
            t_grid = Tbl.y_grid; psi = Tbl.x1_grid;
%             bayes = Tbl.bayes; psi_bimodal = Tbl.psi_bimodal_grid;
            % spline interpolate nonlinear estimate and apply scaling
            disp('The adaptive estimate is')
            t_tilde = interp1(t_grid,psi,tO,'spline');
            adaptive_nonlinear = VUO/sqrt(VO) * t_tilde + CUE;
            disp(adaptive_nonlinear)
            % get the gradient and calculate the SURE
            SURE = VUO^2/VO*((t_tilde - tO)^2 + 2*interp1(t_grid,gradient(psi,t_grid),tO,'spline') - 1) + ...
                VU - VUO^2/VO;
    
            Tbl = readtable(strcat('lookup_tables/const_risk_sigmatb_',corr_str,'_B',string(B),'.csv'));
            b_grid=Tbl.b_grid; Kb = length(b_grid);
            rho_b_over_sigma = interp1(rho_tbl(:,1),rho_tbl(:,2),abs(b_grid),'spline');
            risk_oracle = rho_b_over_sigma + 1/corr_sq  -1;
            risk_function_adaptive = Tbl.risk_function_x1;
%             risk_function_bimodal = Tbl.risk_function_bimodal;
             
%             risk_bimodal = Tbl.risk_bimodal;
            penalty_nonlinear = max(risk_function_adaptive./risk_oracle);
%             penalty_bimodal = max(risk_function_bimodal./risk_bimodal);
            % report adapatation penalty
            disp('The adaptive penalty estimate is')
            disp(penalty_nonlinear)
    %% Calculate the wrost-case risk increase for various estimators
            disp('The adaptive estimator has worst case risk')
            disp(corr_sq*max(risk_function_adaptive ) )
%             disp('The adaptive soft-threshold has worst case risk')
%             disp(corr_str^2*max(risk_function_st_adaptive))
            %% Form thresholded estimates (minimax loss)
            cons = 1/corr_sq - 1; 
            %% Define risk function for soft threshold at l
            Eb = @(l) 1+l^2+...
            (b_grid.^2-1-l^2).*(normcdf(l-b_grid)-normcdf(-l-b_grid))+...
            (-b_grid-l).*normpdf(l-b_grid) - (l-b_grid).*normpdf(-l-b_grid);

            Tbl = readmatrix(strcat('lookup_tables/const_mu_st_sigmatb_',corr_str,'_B',string(B),'.csv'));
            disp('The adaptive soft threshold is')
            st = Tbl(1,1);
            disp(st)
            disp('The adaptive soft-thresholded estimate is')
            adaptive_st = VUO/sqrt(VO) * ((tO > st)*(tO - st) + (tO < -st)*(tO + st)) + CUE;
            disp(adaptive_st)
            SURE_st = VUO^2/VO*( (abs(tO) > st)*(st^2+2) + (abs(tO) < st)*(tO^2) - 1) + ...
                VU - VUO^2/VO;
            disp('The adaptive soft-threshold has penalty')
            penalty_st = max((corr_sq*Eb(st)+1-corr_sq)./(corr_sq*rho_b_over_sigma + 1-corr_sq));
            disp(penalty_st)
            disp('The adaptive soft-thresholded has worst case risk')
            disp(max(corr_sq*Eb(st)+1-corr_sq))
 
            %% put everything in a matrix
            results = zeros(5,2); 
            results(1,1)=adaptive_nonlinear; results(2,1)=penalty_nonlinear;
            results(4,1)=SURE; results(5,1)=corr_sq*max(risk_function_adaptive);
            results(1,2)=adaptive_st; results(2,2)=penalty_st; results(3,2)=st;
            results(4,2)=SURE_st; results(5,2)=max(corr_sq*Eb(st)+1-corr_sq);

            writematrix(results,'const_results.csv')
            results_scale_minimax = results(:)';

        end
    end
end