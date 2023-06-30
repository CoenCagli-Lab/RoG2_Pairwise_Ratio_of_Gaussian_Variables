function [MU_R,SIG_R,CORR_R] = mncovMgCrf(params,cont,mu_eta,tuning)
%% mncovMgCrf Returns the Mean/Covariance Approximations for the modulated Gaussians Model using the contrast response parametrization
%   INPUT
%       pars (1,11) = [Rmax(2) eps(2) sigma_g(2) rho_P rho_G var_eta(2) rho_eta]
%       contrasts = levels of contrast [0-100]
%       mu_eta (2,1) = mean spont. activity
%   OUTPUT
%       mu_r = mean of the modulated Gaussians (2,number of contrasts)
%       sigma_r = covariance of the modulated Gaussians (number of contrasts,2,2)
%       corr_r = correlation of the modulated Gaussians (number of contrasts,1)
%     if ~exist('mu_eta','var')
%         mu_eta = zeros(2,1);
%     end
arguments
        params (1,11);
        cont (1,:);
        mu_eta (2,1);
        tuning function_handle = @contrastResponse;
    end
    [mu_n,mu_d] = tuning(params(1:2),params(3:4),cont);

    MU_R = mu_n./mu_d+reshape(mu_eta,2,1);
    MU_R(MU_R<0) = NaN;
    g = reshape(params(5:6),2,[]);
    v_eta = reshape(params(9:10),2,[]);
    sigma_R2 = MU_R+g.^2.*MU_R.^2+v_eta;
%     sigma_R2(sigma_R2<0) = NaN;
    pMU_R = prod(MU_R,1).^(0.5);
%     pMU_R(imag(pMU_R)~=0) = NaN;
%     pMU_R = real(pMU_R);
    cv = params(7)*pMU_R+params(8)*prod(g.*MU_R,1)+params(11)*prod(v_eta.^(0.5),1);

    nStims = numel(cont);
    SIG_R = NaN(nStims,2,2);
    for iStim=1:nStims
        SIG_R(iStim,:,:) = [   sigma_R2(1,iStim) cv(iStim);...
                    cv(iStim)   sigma_R2(2,iStim);];
    end
    CORR_R = cv./sqrt(prod(sigma_R2,1));
end
