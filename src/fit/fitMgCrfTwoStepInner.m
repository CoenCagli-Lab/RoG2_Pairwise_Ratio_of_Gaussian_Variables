function [x_pair,x_ind] = fitMgCrfTwoStepInner(mu_hat,sigma_hat,cont,mu_eta,var_eta,rho_eta,tuning,opts)
%% fitMgCrfTwoStepInner


%%
        fun_temp = @(x)(mgNegLogLikeParam(x,mu_hat,sigma_hat,cont,mu_eta,tuning));

        [pars_0,lb_pars,ub_pars] = mgOptimBounds(max(mu_hat,[],2),var_eta);
%% ---------------------------------------------------------------------- % 
% % these lines fix rho_eta, var_eta to be the empirically recorded values,
% % and not parameters to be optimized. This makes a quantitative 
% % difference in the quality of fits but not a qualitative difference in 
% % the results.
%% ---------------------------------------------------------------------- % 
%       pars_0(10:12) = [reshape(var_eta,1,2),rho_eta];
%       lb_pars(10:12) = [reshape(var_eta,1,2),rho_eta];
%       ub_pars(10:12) = [reshape(var_eta,1,2),rho_eta];
%%
%% If negative log-likelihood returns NaN, computation will error out
% % Continue to next loop
        if isnan(fun_temp(pars_0))
            x_pair = NaN(size(pars_0));
            x_ind = NaN(size(pars_0));
            return
        end
%% Optimize the parameters that are independent of correlation structrue (all besides rho parameters)
        x_ind = fmincon(fun_temp,pars_0,[],[],[],[],lb_pars,ub_pars,[],opts);
        % % Fix these parameters, change the range for the rho parameters
        lb_pars = x_ind;
        ub_pars = x_ind;
        x_ind(end) = rho_eta;
        lb_pars([7:8 11]) = [-1 -1 rho_eta];
        ub_pars([7:8 11]) = [1 1 rho_eta];
%% ---------------------------------------------------------------------- % 
% %  IMPORTANT - Currently set so that rho_eta == 0, as the data analyzed 
% %  for this paper did not have this. You should change out these lines if 
% %  you have recorded correlations in spontaneous activity.
% ----------------------------------------------------------------------- %
%%
%         lb_pars(end) = -1;
%         ub_pars(end) = 1; 
%%        
%% Optimize full model
        x_pair = fmincon(fun_temp,x_ind,[],[],[],[],lb_pars,ub_pars,[],opts);

end