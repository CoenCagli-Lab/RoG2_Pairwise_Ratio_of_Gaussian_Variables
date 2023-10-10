function [x_pair,x_ind] = fitMgCrfOneStep(spikes,cont,mu_eta,cov_eta,tuning,opts)
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
[x_pair,x_ind] = fitMgCrfOneStepInner(mu_hat,sigma_hat,cont,mu_eta,var_eta,rho_eta,tuning,opts);
end