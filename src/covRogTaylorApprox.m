function [covResponse,noiseCorrResponse] = covRogTaylorApprox(alpha_n,beta_n,mu_n,alpha_d,beta_d,mu_d,rho,cov_eta)
%% covRogTaylorApprox Creates Covariance Matrix for a Ratio of Gaussians Model using the Taylor Series Approximation
% Assumes variances for individual numerator/denominator variables are parametrized by a power law relationship to the mean
% Covariance Between R_1 = N_1/D_1 and R_2 = N_2/D_2
%   INPUT
%       alpha_n (2,1) - alpha parameter for numerators in power law of variance
%       beta_n (2,1) - beta parameter for numerators in power law of variance
%       mu_n (2,number of stimuli) - mean for numerators
%       alpha_d (2,1) - alpha parameter for denominators in power law of variance
%       beta_d (2,1) - beta parameter for denominators in power law of variance
%       mu_d (2,number of stimuli) - mean for denominators
%       rho (4,1) - rho parameters [rho_n,rho_d,rho_n1d2,rho_n2d1]
%       cov_eta (2,2) - covariance of additive noise
%   OUTPUT
%       covResponse (number of stimuli,2,2) - covariance of the Ratio of Gaussians model
%       noiseCorrResponse (number of stimuli,1) - noise correlation for Ratio of Gaussians model
%%

    var_d = varPowerLaw(alpha_d,beta_d,mu_d);
    var_n = varPowerLaw(alpha_n,beta_n,mu_n);

    rvar = varRogTaylorApprox(mu_n,var_n,mu_d,var_d)+diag(cov_eta);

    t1 = mu_n(1,:).*mu_n(2,:);
    t2 = mu_d(1,:).*mu_d(2,:);
    t3 = rho(1)*sqrt(var_n(1,:).*var_n(2,:));
    t4 = rho(2)*sqrt(var_d(1,:).*var_d(2,:));
    t5 = rho(3)*sqrt(var_n(1,:).*var_d(2,:));
    t6 = mu_n(1,:).*mu_d(2,:);
    t7 = rho(4)*sqrt(var_n(2,:).*var_d(1,:));
    t8 = mu_n(2,:).*mu_d(1,:);

   rcov = t1./t2.*(t3./t1+t4./t2-t5./t6-t7./t8)+cov_eta(1,2);
%    if size(rvar,1) ~= 2
%        rvar = rvar.';
%    end

%    if size(rcov,1) ~= 1
%        rcov = rcov.';
%    end

    stim_num = size(mu_n,2);
    covResponse = NaN(stim_num,2,2);
    noiseCorrResponse = NaN(stim_num,1);
    for ind=1:stim_num
        covResponse(ind,:,:) = [rvar(1,ind) rcov(ind); rcov(ind) rvar(2,ind)];
%         covResponse(ind,:,:) = diag(rvar(:,ind))+diag(rcov(ind),-1)+diag(rcov(ind),1);
        noiseCorrResponse(ind) = rcov(ind)/((rvar(1,ind).*rvar(2,ind))^(0.5));
    end

end
