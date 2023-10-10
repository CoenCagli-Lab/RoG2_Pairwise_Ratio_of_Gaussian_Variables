function [x_pair,x_ind] = fitRogCrfOneStepInner(mu_hat,sigma_hat,cont,mu_eta,var_eta,rho_eta,tuning,opts)
% FITROGONESTEPINNER Inner loop for fitting the RoG in one step

        fun_temp = @(x)(rogNegLogLikeParam(x,mu_hat,sigma_hat,cont,mu_eta,tuning));

        [pars_0,lb_pars,ub_pars] = rogOptimBounds(max(mu_hat,[],2),var_eta);
%% ---------------------------------------------------------------------- % 
% % these lines fix rho_eta, var_eta to be the empirically recorded values,
% % and not parameters to be optimized. This makes a quantitative 
% % difference in the quality of fits but not a qualitative difference in 
% % the results.
%% ---------------------------------------------------------------------- % 
%       pars_0(17:19) = [reshape(var_eta,1,2),rho_eta];
%       lb_pars(17:19) = [reshape(var_eta,1,2),rho_eta];
%       ub_pars(17:19) = [reshape(var_eta,1,2),rho_eta];
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
        % lb_pars([1:12 17:18]) = x_ind([1:12 17:18]);
        % ub_pars([1:12 17:18]) = x_ind([1:12 17:18]);
%% ---------------------------------------------------------------------- % 
% %  We only considered the case where rho_n1d2=rho_n2d1=0
% %  Change these lines to look at the more general case
% % ----------------------------------------------------------------------- %
%        lb_pars(13:16) = [-1 -1 -1 -1];
%        ub_pars(13:16) = [1 1 1 1]; 
        pars_0(end) = rho_eta;
        lb_pars([13:16 19]) = [-1 -1 0 0 rho_eta]; 
        ub_pars([13:16,19]) = [1 1 0 0 rho_eta]; 
        x_ind(end) = rho_eta;
%% ---------------------------------------------------------------------- % 
% %  IMPORTANT - Currently set so that rho_eta == 0, as the data analyzed 
% %  for this paper did not have this. You should change out these lines if 
% %  you have recorded correlations in spontaneous activity.
% ----------------------------------------------------------------------- %
%%
%         lb_pars(end) = -1;
%         ub_pars(end) = 1; 
%%        
%% Optimize the pairwise rho parameters

        x_pair = fmincon(fun_temp,pars_0,[],[],[],[],lb_pars,ub_pars,[],opts);
end