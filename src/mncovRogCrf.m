function [mu_r,sigma_r,corr_r] = mncovRogCrf(pars,contrasts,mu_eta,tuning)
%% mncovRogCrf Returns the Mean/Covariance Approximations for the pairwise Ratio of Gaussians Model using the contrast response parametrization
%   INPUT
%       pars (1,19) = [Rmax(2) eps(2) alphan1 betan1 alphad1 betad1 alphan2 betan2 alphasd2 betad2 rhon rhod rhon1d2 rhon2d1 var_eta(2) rho_eta]
%       contrasts = levels of contrast [0-100]
%       mu_eta (2,1) = mean spont. activity
%   OUTPUT
%       mu_r = mean of the Ratio of Gaussians (2,number of contrasts)
%       sigma_r = covariance of the Ratio of Gaussians (number of contrasts,2,2)
%       corr_r = correlation of the Ratio of Gaussins (number of contrasts,1)
%%
%     if ~exist('mu_eta','var')
%         mu_eta = zeros(2,1);
%     end
    arguments
        pars (1,19);
        contrasts (1,:);
        mu_eta (2,1);
        tuning function_handle = @contrastResponse;
    end

    [mu_n,mu_d] = tuning(pars(1:2),pars(3:4),contrasts);
%     [mu_n,mu_d] = contrastResponse(pars(1:2),pars(3:4),contrasts);

    c_eta = pars(19)*sqrt(pars(17)*pars(18));
%     sigma_eta = diag(pars([17 18]))+diag(c_eta,1)+diag(c_eta,-1);
    sigma_eta = [pars(17) c_eta; c_eta pars(18)];

    mu_r = mu_n./mu_d+reshape(mu_eta,2,1);
    mu_r(mu_r<0) = NaN;
    [sigma_r,corr_r] = covRogTaylorApprox(...
                                    pars([5 9]).',...
                                    pars([6 10]).',...
                                    mu_n,...
                                    pars([7 11]).',...
                                    pars([8 12]).',...
                                    mu_d,...
                                    pars(13:16),...
                                    sigma_eta);
end
