function [penalty_nonlinear, penalty_st, penalty_ttest] = calculate_max_regret(corr)
     Sigma_UO_grid = abs(tanh((-3:0.05:-0.05)));  % take out zero correlation coeff
        if abs(corr) > max(Sigma_UO_grid) || abs(corr) < min(Sigma_UO_grid)
            error("Extrapolating beyond the grid! " + ...
                "Pre-tabulations are available for absolute values of " + ...
                "correlation coefficient between " + ...
            min(Sigma_UO_grid) + " and " + max(Sigma_UO_grid) + ". " + ...
            "Optimize numerically for the given rho value.")
        end
        rho_tbl = readmatrix('lookup_tables/minimax_rho_B9.csv');

	%% Risk function
    load('lookup_tables/risk.mat');
    Kb = length(b_grid);
    risk_function_adaptive = zeros(Kb,1);
    
    for i = 1:Kb
        risk_function_adaptive(i) = interp1(Sigma_UO_grid, risk_mat(i,:), abs(corr), 'spline');
    end
    
    % calculate the oracle risk function
    rho_b_over_sigma = interp1(rho_tbl(:,1), rho_tbl(:,2), abs(b_grid), 'spline');
    risk_oracle = rho_b_over_sigma + 1 / corr^2 - 1;

    max_regret_nonlinear = max(risk_function_adaptive ./ risk_oracle);
    disp('The worst-case adaptation regret is');
    disp(max_regret_nonlinear);
    
    %% Add soft threshold risk functions
    Eb = @(l) 1 + l^2 + ...
        (b_grid.^2 - 1 - l^2) .* (normcdf(l - b_grid) - normcdf(-l - b_grid)) + ...
        (-b_grid - l) .* normpdf(l - b_grid) - (l - b_grid) .* normpdf(-l - b_grid);
    risk_function_st_adaptive = Eb(st) + 1 / corr^2 - 1;
    
    load('lookup_tables/thresholds.mat');
    disp('The adaptive soft threshold is');
    st = interp1(Sigma_UO_grid, st_mat, abs(corr), 'spline');
    disp(st);
    
    max_regret_st = max(risk_function_st_adaptive ./ risk_oracle);
    disp('The adaptive soft-threshold has max regret');
    disp(max(risk_function_st_adaptive ./ risk_oracle));
    
    %% Use simulation to calculate the risk function for the pre-test estimator that switches between Y_U and Y_R
    rng(1);
    sims = 100000;
    x = normrnd(0,1,[sims,1]);
    x_b = x * ones(1, Kb) + ones(sims, 1) * b_grid';
    Ebsims_ht = @(l) sum(((x_b > l) .* x_b + (x_b < l & x_b > -l) * (1 + VO / VUO) .* x_b ...
        + (x_b < -l) .* x_b - ones(sims, 1) * b_grid').^2, 1) / sims;
    risk_function_ht_ttest = (Ebsims_ht(1.96) + 1 / corr^2 - 1)';
    
    %% Calculate max regret for various estimators
    max_regret_ttest = max(risk_function_ht_ttest ./ risk_oracle);
    disp('The pre-test has worst-case adaptation regret');
    disp(max_regret_ttest);
end
