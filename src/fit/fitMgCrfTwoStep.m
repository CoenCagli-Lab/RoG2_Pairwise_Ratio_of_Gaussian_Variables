function [x_pair,x_ind] = fitMgCrfTwoStep(spikes,cont,mu_eta,cov_eta,optional)
arguments
    spikes (2,:,:);
    cont (1,:);
    mu_eta (2,1);
    cov_eta (2,2);
    optional.tuning function_handle = @contrastResponse;
    optional.optimopts optim.options.Fmincon = optimoptions(@fmincon,...
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

var_eta = diag(cov_eta);
if any(var_eta == 0)
    rho_eta = 0;
else
    rho_eta = cov_eta(1,2)/((prod(var_eta))^(0.5));
end

%% Empirical (hat) and Null mean and covariance
[mu_hat,sigma_hat,~] = mncovEmpirical(spikes);
[x_pair,x_ind] = fitMgCrfTwoStepInner(mu_hat,sigma_hat,cont,mu_eta,var_eta,rho_eta,optional.tuning,optional.optimopts);
end