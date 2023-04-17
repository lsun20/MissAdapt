function minimax_locus_plot(YR,YU,VR,VU,VUR,corr_str,B_adaptive_grid)
% YR is the restricted estimate
% YU is the unbiased estimate
% B_adaptive_grid is the (optional) upper bound on the scaled bias e.g. [0.5 1 2];
if nargin < 7 || isempty(B_adaptive_grid)
     B_adaptive_grid = 9;
end
global corr_sq %if passing in corr_str as a string
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

        rho_tbl = readmatrix('lookup_tables/minimax_rho_B9.csv');

        for i = 1:length(B_adaptive_grid)
            B = B_adaptive_grid(i);
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
            %% Form locus of B-minimax estimates based on pre-tabulated values
%             load('lookup_tables/minimax_bounded_normal_mean.mat');
%             Y_minimax = zeros(length(B_grid),1);
%             for i = 1:length(B_grid)
%                 t_minimax = interp1(y_grid,psi_function_minimax(:,i),tO,'spline');
%                 Y_minimax(i) = VUO/sqrt(VO) * t_minimax + CUE;
%             end
%             B_grid = [0 B_grid];
%             Y_minimax = [CUE Y_minimax']';
            %% Form locus of B-minimax estimates based on posterior means
            B=9;
            global b_grid y_grid Ky Kb
            b_grid = (-B:0.05:B)';
            y_grid = [tO];
            Ky = length(y_grid);
            Kb = length(b_grid);
            Sigma = [1 1;1 1];
            %Pi = pmf(Sigma);
            y_grid_ub = y_grid+0.05; y_grid_lb = y_grid-0.05;
            Pi = normcdf((y_grid_ub+y_grid)/2*ones(1,Kb),ones(Ky,1)*b_grid',1) - ...
             normcdf((y_grid_lb+y_grid)/2*ones(1,Kb),ones(Ky,1)*b_grid',1); 
            Psi = @(mu_grid) (Pi*(mu_grid.*(b_grid))) ./ (Pi*mu_grid);
            load('lookup_tables/priors_bounded_normal_mean.mat');
            t_minimax = zeros(length(B_grid),1);
            for i = 1:length(B_grid)
                t_minimax(i) = Psi(prior_bounded_normal_mean_mat(:,i)); % minimax estimator
                
            end
            Y_minimax = VUO/sqrt(VO) * t_minimax + CUE;
            B_grid = [0 B_grid];
            Y_minimax = [CUE Y_minimax']';
            %% Form adaptive risk
            Kb = length(B_grid); risk_function_adaptive = zeros(Kb,1);
            load('lookup_tables/risk.mat');
            for i = 1:Kb
               risk_function_adaptive(i) = interp1(Sigma_UO_grid,risk_mat(abs(b_grid-B_grid(i))<0.001,:),abs(corr(1)),'spline');
            end
            risk_function_adaptive = corr_str^2*risk_function_adaptive;
            %% Form oracle risk
            rho_b_over_sigma = interp1(rho_tbl(:,1),rho_tbl(:,2),abs(B_grid),'spline');
            risk_oracle = corr_str^2*rho_b_over_sigma + 1 - corr_str^2;
            %% plot the locus and the riskfig = figure(1)
            clear fig
            fig = figure(1)
            B_grid_scale = sqrt(B_grid);
            yyaxis left
            yline([YR CUE adaptive_nonlinear YU],'-',{'Y_R','GMM','Adaptive','Y_U'},...
                'LabelHorizontalAlignment','left','FontSize',14 )
%             yline([YR CUE YU],'-',{'Y_R','GMM','Y_U'},...
%                 'LabelHorizontalAlignment','left','FontSize',14 )
            hold on
%             p1=plot(B_grid_scale,Y_minimax,'--',...
%                 'DisplayName','B-minimax estimates','LineWidth',2 );
            p1=scatter(B_grid_scale,Y_minimax,50,'k','^',...
                'DisplayName','B-minimax estimates');             
            
            xlabel('B/\Sigma_O^{1/2}');  
            xticks([0 1 2 3])
            xticklabels({'0','1','4','9'})
            ylabel('Estimate','FontSize',14 );  
            
            yyaxis right
                     
            p2=plot(B_grid_scale,risk_function_adaptive,':',...
                'DisplayName','Adaptive risk','LineJoin','miter','LineWidth',2 );
            
            p3=plot(B_grid_scale,risk_oracle,'--',...
                'DisplayName','Oracle risk','LineJoin','miter','LineWidth',2 );
            hold on   
            yline(1,'--','Color',[0.5 0.5 0.5])
            
            ylabel('Mean squared error relative to Y_U','FontSize',14,'rotation',-90,'VerticalAlignment','bottom'); 
            ylim([min(risk_oracle)-0.1 max(risk_function_adaptive)]+0.1)
            hold off
            
            legend([p1 p2 p3]...
                ,'Location','southoutside','Orientation','horizontal',...
                'FontSize',14 ,'NumColumns', 3);
            legend('boxoff');
            %% Save figures
            figurename = strcat( 'minimax_locus_sigmatb_',...
                string(round(abs(corr(1))*100)/100),'_B',string(B),'.png');
            set(fig,'Units','Inches');
            pos = get(fig,'Position');
            set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
            print(fig,figurename,'-dpng','-r0')
            
        end
    else
%% Loop over upper bounds
        for i = 1:length(B_grid)
            B = B_grid(i);
            B
            %% Form nonlinear adaptive estimates
            % Looks for the nonlinear estimates stored in the /lookup_tables/
            % subdirectory
    
            Tbl = readtable(strcat('lookup_tables/risk_and_oracle_risk_sigmatb_',corr_str,'_B',string(B),'.csv'));
            risk_function_adaptive =  VUO^2/VO* Tbl.risk_function_adaptive;
            risk_function_oracle = VUO^2/VO* Tbl.risk_oracle;
            
            Tbl = readtable(strcat('lookup_tables/risk_st_and_oracle_risk_sigmatb_',corr_str,'_B',string(B),'.csv'));
            bias_grid = Tbl.bias_grid; Kb = length(bias_grid);
    
            risk_function_st_adaptive = VUO^2/VO*  Tbl.risk_function_st_adaptive;
            risk_function_ht_ttest = VUO^2/VO*  Tbl.risk_function_ht_ttest;
            risk_function_ht_adaptive = VUO^2/VO*  Tbl.risk_function_ht_adaptive;
    
            load('lookup_tables/additive.mat');
            risk_function_add =  VUO^2/VO* (risk_function_add+ 1/corr_sq -1);
            %% report relative risk results
            results = zeros(3,1); results(2,1) = risk_function_adaptive(bias_grid == 0)/VU;
            results(1,1) = max(risk_function_adaptive/VU);
            breakeven = interp1(risk_function_adaptive(bias_grid>0)/VU,bias_grid(bias_grid>0),1,'spline');
            results(3,1) = breakeven;
            results_scale_minimax = results(:)'
            %% plot adaptive and oracle risk and adaptation penalty (ratio of adaptive to oracle risk)
            fig = figure(1)
    %         plot(bias_grid,risk_function_adaptive,':',...
    %             bias_grid,risk_function_st_adaptive,'-.',...
    %             bias_grid,risk_function_oracle,'-',...
    %             bias_grid,(risk_function_adaptive)./(risk_function_oracle),':',... 
    %             bias_grid,(risk_function_st_adaptive)./(risk_function_oracle),'-.',... 
    %         'LineJoin','miter','LineWidth',2 );
    %         legend({'Adaptive risk','Soft-threshold risk','Oracle risk',...
    %             'Adaptation penalty','Soft-threshold adaptation penalty'},'Location','southoutside','Orientation','horizontal',...
    %             'FontSize',14 );
    %         legend('boxoff');
    %         xlabel('b/\Sigma_O^{1/2}'); xlim([min(bias_grid) max(bias_grid)])
    %         ylabel('Mean squared error','FontSize',14 );
            b_grid_scale = sign(bias_grid).*sqrt(abs(bias_grid));
            plot(b_grid_scale,risk_function_adaptive/VU,':',...
                 b_grid_scale,risk_function_add/VU,'-',...
                b_grid_scale,risk_function_st_adaptive/VU,'-.',...
                b_grid_scale,risk_function_ht_ttest/VU,'-.',...
                b_grid_scale,risk_function_oracle/VU,'-',...
                b_grid_scale,ones(Kb,1),'--',...
            'LineJoin','miter','LineWidth',2 );
            legend({'Adaptive multiplicative','Adaptive additive',...
                'Soft-threshold','Pre-test','Oracle','Y_U'}...
                ,'Location','southoutside','Orientation','horizontal',...
                'FontSize',14 ,'NumColumns', 3);
            legend('boxoff');
            
            xlabel('b/\Sigma_O^{1/2}'); xlim([min(b_grid_scale) max(b_grid_scale)])
            ylabel('Mean squared error','FontSize',14 ); ylim([0.9 1.1])
            xticks([-3 -2 -1 0 1 2 3])
            xticklabels({'-9','-4','-1','0','1','4','9'})
            ylabel('Mean squared error relative to \Sigma_U','FontSize',14 );
            %% Save figures
            figurename = strcat( 'risk_plot_sigmatb_',corr_str,'_B',string(B),'.png');
            set(fig,'Units','Inches');
            pos = get(fig,'Position');
            set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
            print(fig,figurename,'-dpng','-r0')
        end
    end

end

