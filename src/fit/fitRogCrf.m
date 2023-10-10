function [x_pair] = fitRogCrf(spikes,cont,mu_eta,cov_eta,tuning,opts)
    arguments
        spikes (2,:,:);
        cont (1,:);
        mu_eta (2,1);
        cov_eta (2,2);
        tuning function_handle = @contrastResponse;
        opts optim.options.Fmincon = optimoptions(@fmincon,...
                        'Algorithm','sqp',...
                        'MaxIter',1000,...
                        'Display','off',...
                        'TolX',1e-6,...
                        'TolFun',1e-6,...
                        'TolCon',1e-3);
    end

    if size(spikes,2) ~= size(cont(:))
        error('Number of Contrasts does not match!')
    end

%     if isequal(size(cov_eta),[2,2])
        var_eta = diag(cov_eta);
        if any(var_eta == 0)
            rho_eta = 0;
        else
            rho_eta = cov_eta(1,2)/((prod(var_eta))^(0.5));
        end
%     else
%         error("Must input covariance matrix for noise")
%     end


        %% Empirical (hat) and Null mean and covariance
        [mu_hat,sigma_hat,~] = mncovEmpirical(spikes);
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
        if isnan(fun_temp(pars_0))
            x_pair = NaN(size(pars_0));
            return 
        end
%% Optimize the parameters that are independent of correlation structrue (all besides rho parameters)
        x_ind = fmincon(fun_temp,pars_0,[],[],[],[],lb_pars,ub_pars,[],opts);
%         [x_ind,~,~,~,~,~,hessInd] = fmincon(fun_temp,pars_0,[],[],[],[],lb_pars,ub_pars,[],opts);

        % % Fix these parameters, change the range for the rho parameters
        lb_pars([1:12 17:18]) = x_ind([1:12 17:18]);
        ub_pars([1:12 17:18]) = x_ind([1:12 17:18]);
%% ---------------------------------------------------------------------- % 
% %  We only considered the case where rho_n1d2=rho_n2d1=0
% %  Change these lines to look at the more general case
% % ----------------------------------------------------------------------- %
%        lb_pars(13:16) = [-1 -1 -1 -1];
%        ub_pars(13:16) = [1 1 1 1]; 

        lb_pars([13:16 19]) = [-1 -1 0 0 rho_eta]; % change rho_eta to -1
        ub_pars([13:16,19]) = [1 1 0 0 rho_eta]; % change rho_eta to 1
%% ---------------------------------------------------------------------- % 
% %  IMPORTANT - Currently set so that rho_eta == 0, as the data analyzed 
% %  for this paper did not have this. You should change out these lines if 
% %  you have recorded correlations in spontaneous activity.
% ----------------------------------------------------------------------- %
%%
%         x_ind(end) = rho_eta;
%         lb_pars(end) = -1;
%         ub_pars(end) = 1; 
%%        
%% Optimize the pairwise rho parameters
    x_pair = fmincon(fun_temp,x_ind,[],[],[],[],lb_pars,ub_pars,[],opts);
%% Alternative: Minimize least square between parametric noise correlation function and true
%%
        % fun_lsq = @(pars)(lsqNoiseCorr(pars,cont,mu_eta,corr_hat,tuning));
        % x_pair = fmincon(fun_lsq,x_ind,[],[],[],[],lb_pars,ub_pars,[],opts);
%         prob = createOptimProblem('fmincon','x0',x_ind,'objective',fun_temp,'lb',lb_pars,'ub',ub_pars,'options',opts);
%         gs = GlobalSearch('Display','off');
%         x_pair = run(gs,prob);
%         [x_pair,~,~,~,~,~,hessPair] = fmincon(fun_temp,x_ind,[],[],[],[],lb_pars,ub_pars,[],opts);
%         hess = struct('pair',hessPair,'ind',hessInd);
%         hess1D = sqrt(diag(inv(hess1D)));
%         hess = sqrt(diag(inv(hess)));
%         hess([1:12 17:19]) = hess1D([1:12 17:19]);
end
% function out = lsqNoiseCorr(params,contrasts,mu_eta,corrEmp,tuning)
%     [~,~,noiseCorrFit] = mncovRogCrf(params,contrasts,mu_eta,tuning);
%     out = sum((noiseCorrFit-corrEmp).^2,"all");
% end