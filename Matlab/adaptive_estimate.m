function [results] = ...
    adaptive_estimate(YR,YU,VR,VU,VUR,B_grid,corr_str)
% YR is the restricted estimate
% YU is the unbiased estimate
% B_grid is the upper bound on the scaled bias e.g. [0.5 1 2];
% corr_str can be a string  e.g. 0.52 = '052' that is used to look up the
% csv
% corr_str can also be the numerical corr coef that is used to interpolate
% the estimates
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
    disp('The correlation coefficient is')
    corr = VUO/sqrt(VO)/sqrt(VU);
    disp(corr)
    
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
        rho_tbl = readmatrix('lookup_tables/minimax_rho_B9.csv');

        for i = 1:length(B_grid)
            B = B_grid(i);
            B

            %% Form nonlinear adaptive estimates
            % Looks for the nonlinear estimates stored in the /lookup_tables/
            % subdirectory
            
            load('lookup_tables/policy.mat');          
            % spline interpolate nonlinear estimate and apply scaling
            [X_grid,Y_grid] = meshgrid(Sigma_UO_grid,y_grid);
            disp('The adaptive estimate is')   
            t_tilde = interp2(X_grid,Y_grid,psi_mat,abs(corr_str),tO,'spline');
            adaptive_nonlinear = VUO/sqrt(VO) * t_tilde + CUE;
            % get the gradient and calculate the SURE
            Ky = length(y_grid); psi = zeros(Ky,1);
            for i = 1:Ky
               psi(i) = interp1(Sigma_UO_grid,psi_mat(i,:),abs(corr_str),'spline');
            end
            SURE = VUO^2/VO*((t_tilde - tO)^2 + 2*interp1(y_grid,gradient(psi,y_grid),tO,'spline') - 1) + ...
                VU - VUO^2/VO;
            disp(interp1(y_grid,gradient(psi,y_grid),tO,'spline'))
            
         
            %% Risk function
            load('lookup_tables/risk.mat');
            Kb = length(b_grid); risk_function_adaptive = zeros(Kb,1);
            for i = 1:Kb
               risk_function_adaptive(i) = interp1(Sigma_UO_grid,risk_mat(i,:),abs(corr_str),'spline');
            end
            % calculate the oracle risk function
            rho_b_over_sigma = interp1(rho_tbl(:,1),rho_tbl(:,2),abs(b_grid),'spline');
            risk_oracle = rho_b_over_sigma + 1/corr_str^2 -1;

            
            penalty_nonlinear = max(risk_function_adaptive./risk_oracle);
            % report adapatation penalty
            disp('The worst-case adaptation regret is')
            disp(penalty_nonlinear)

            %% Form thresholded estimates (minimax loss)
            load('lookup_tables/thresholds.mat');
            disp('The adaptive soft threshold is')
            st = interp1(Sigma_UO_grid,st_mat,abs(corr_str),'spline');
            disp(st)
            disp('The adaptive soft-thresholded estimate is')
            adaptive_st = VUO/sqrt(VO) * ((tO > st)*(tO - st) + (tO < -st)*(tO + st)) + CUE;
            disp(adaptive_st)
            SURE_st = VUO^2/VO*( (abs(tO) > st)*(st^2+2) + (abs(tO) < st)*(tO^2) - 1) + ...
                VU - VUO^2/VO;
        
            disp('The adaptive hard threshold is')
            ht = interp1(Sigma_UO_grid,ht_mat,abs(corr_str),'spline');
            disp(ht)
            disp('The adaptive hard-thresholded estimate is')
            adaptive_ht = VUO/sqrt(VO) * ((tO > ht)*(tO) + (tO < -ht)*(tO)) + CUE;
            disp(adaptive_ht)

            disp('The pre-test (1.96) estimate is')
            %pretest_ht = VUO/sqrt(VO) * ((tO > 1.96)*(tO) + (tO < -1.96)*(tO)) + CUE;
            pretest_ht = YR - sqrt(VO) * ((tO > 1.96)*(tO) + (tO < -1.96)*(tO));
            disp(pretest_ht)

            load('lookup_tables/risk_thresholds.mat');
            Kb = length(b_grid); 
            risk_function_st_adaptive = zeros(Kb,1);
            risk_function_ht_adaptive = zeros(Kb,1);
            risk_function_ht_ttest = zeros(Kb,1);
            for i = 1:Kb
               risk_function_st_adaptive(i) = interp1(Sigma_UO_grid,risk_st_mat(i,:),abs(corr_str),'spline');
               risk_function_ht_adaptive(i) = interp1(Sigma_UO_grid,risk_ht_mat(i,:),abs(corr_str),'spline');
               risk_function_ht_ttest(i) = interp1(Sigma_UO_grid,risk_ht_ttest_mat(i,:),abs(corr_str),'spline');
            end
            %% Use similation to calculate the risk function for the pre-test estimator that switches btw Y_U and Y_R
            sims = 100000;
            x = normrnd(0,1,[sims,1]);
            x_b = x*ones(1,Kb) + ones(sims,1)*b_grid';
            Ebsims_ht = @(l) sum(((x_b > l).*x_b + (x_b < l & x_b > -l)*(1+VO/VUO).*x_b + (x_b < -l).*x_b...
                -ones(sims,1)*b_grid').^2,1)/sims;
            risk_function_ht_ttest = (Ebsims_ht(1.96) + 1/corr_str^2 - 1)';
            %% Calculate penalty for various estimators
            disp('The adaptive soft-threshold has penalty')
            disp(max(risk_function_st_adaptive./risk_oracle))
            disp('The adaptive hard-threshold has penalty')
            disp(max(risk_function_ht_adaptive./risk_oracle))
            disp('The pre-test has worst-case adaptation regret')
            disp(max(risk_function_ht_ttest./risk_oracle))
            %% put everything in a matrix
            results = zeros(4,3); 
            results(1,1)=adaptive_nonlinear; results(2,1)=penalty_nonlinear;
            results(4,1)=SURE;
            results(1,2)=adaptive_st; results(2,2)=max(risk_function_st_adaptive./risk_oracle); 
            results(3,2)=st; results(4,2)=SURE_st;
            %results(1,3)=adaptive_ht; results(2,3)=max(risk_function_ht_adaptive./risk_oracle); results(3,3)=ht;
            results(1,3)=pretest_ht; results(2,3)=max(risk_function_ht_ttest./risk_oracle); results(3,3)=1.96;
 

        end
    else
    %% Loop over upper bounds - if look up the corr_str file
        for i = 1:length(B_grid)
            B = B_grid(i);
            B
            %% Form nonlinear adaptive estimates
            % Looks for the nonlinear estimates stored in the /lookup_tables/
            % subdirectory
            
            Tbl = readtable(strcat('lookup_tables/minimax_adaptive_psi_sigmatb_',corr_str,'_B',string(B),'.csv'));
            t_grid = Tbl.y_grid; psi = Tbl.psi;
            bayes = Tbl.bayes;  
            % spline interpolate nonlinear estimate and apply scaling
            disp('The adaptive estimate is')
            t_tilde = interp1(t_grid,psi,tO,'spline');
            adaptive_nonlinear = VUO/sqrt(VO) * t_tilde + CUE;
            disp(adaptive_nonlinear)
            % get the gradient and calculate the SURE
            SURE = VUO^2/VO*((t_tilde - tO)^2 + 2*interp1(t_grid,gradient(psi,t_grid),tO,'spline') - 1) + ...
                VU - VUO^2/VO;
    
            Tbl = readtable(strcat('lookup_tables/risk_and_oracle_risk_sigmatb_',corr_str,'_B',string(B),'.csv'));
            b_grid=Tbl.b_grid; Kb = length(b_grid);
            risk_function_adaptive = Tbl.risk_function_adaptive; 
            risk_oracle = Tbl.risk_oracle; 
            penalty_nonlinear = max(risk_function_adaptive./risk_oracle); 
            % report adapatation penalty
            disp('The worst-case adaptation regret is')
            disp(penalty_nonlinear)
    
            %% Form thresholded estimates (minimax loss)
            Tbl = readmatrix(strcat('lookup_tables/minimax_mu_st_ht_sigmatb_',corr_str,'_B',string(B),'.csv'));
            disp('The adaptive soft threshold is')
            st = Tbl(1,1);
            disp(st)
            disp('The adaptive soft-thresholded estimate is')
            adaptive_st = VUO/sqrt(VO) * ((tO > st)*(tO - st) + (tO < -st)*(tO + st)) + CUE;
            disp(adaptive_st)
            SURE_st = VUO^2/VO*( (abs(tO) > st)*(st^2+2) + (abs(tO) < st)*(tO^2) - 1) + ...
                VU - VUO^2/VO;
            disp('The adaptive soft-threshold has penalty')
            disp(Tbl(2,1))
            disp('The adaptive hard threshold is')
            ht = Tbl(4,1);
            disp(ht)
            disp('The adaptive hard-thresholded estimate is')
            adaptive_ht = VUO/sqrt(VO) * ((tO > ht)*(tO) + (tO < -ht)*(tO)) + CUE;
            disp(adaptive_ht)
            disp('The adaptive hard-threshold has penalty')
            disp(Tbl(5,1))
            disp('The pre-test (1.96) estimate is')
            %pretest_ht = VUO/sqrt(VO) * ((tO > 1.96)*(tO) + (tO < -1.96)*(tO)) + CUE;
            pretest_ht = YR - sqrt(VO) * ((tO > 1.96)*(tO) + (tO < -1.96)*(tO));
            disp('The pre-test has worst-case adaptation regret')
            %disp(Tbl(8,1)) %switches between Y_U and GMM
            %% Use similation to calculate the risk function for the pre-test estimator that switches btw Y_U and Y_R
            sims = 100000;
            x = normrnd(0,1,[sims,1]);
            x_b = x*ones(1,Kb) + ones(sims,1)*b_grid';
            Ebsims_ht = @(l) sum(((x_b > l).*x_b + (x_b < l & x_b > -l)*(1+VO/VUO).*x_b + (x_b < -l).*x_b...
                -ones(sims,1)*b_grid').^2,1)/sims;
            risk_function_ht_ttest = (Ebsims_ht(1.96) + 1/corr^2 - 1)';
            % overwrite the pre-tabulated penalty term
            Tbl(8,1) = max(risk_function_ht_ttest./risk_oracle);
            %% put everything in a matrix
            results = zeros(4,3); 
            results(1,1)=adaptive_nonlinear; results(2,1)=penalty_nonlinear;
            results(4,1)=SURE;
            results(1,2)=adaptive_st; results(2,2)=Tbl(2,1); results(3,2)=st;
            results(4,2)=SURE_st;
            %results(1,3)=adaptive_ht; results(2,3)=Tbl(5,1); results(3,3)=ht;
            results(1,3)=pretest_ht; results(2,3)=Tbl(8,1); results(3,3)=1.96;
 
          
        end
    end
end