function nll = mgLogLikeParam(pars,MU_HAT,SIGMA_HAT,contrasts,mu_eta,tuning)
%% mgLogLikeParam Computes the parametrized negative Log-likelihood of a bivarate Gaussian model with Poisson-Gamma inspired moments.
% This model is an adaptation of Goris (2014) to allow for full
% specification of the joint distribution. Parametrized by the contrast
% response function with an adaptive gain for overdispersion.
% Moments are identical to the Goris (2014) moments:
%   mu_R_i = R_i^max*c^2/(epsilon_i^2+c^2)
%   sigma_R_i^2 = mu_R_i+g_i^2*mu_r_i^2
%   cov(R_1,R_2) = rho_p*(mu_R_1*mu_R_2)^(1/2)+rho_g*g_1*g_2*mu_R_1*mu_R_2
% 
%%
    arguments
        pars (1,11);
        MU_HAT (2,:,:);
        SIGMA_HAT (:,2,2);
        contrasts (1,:);
        mu_eta (2,1);
        tuning function_handle = @contrastResponse;
    end
    [MU,SIGMA] = mncovMgCrf(pars,contrasts,mu_eta,tuning);

    v1 = SIGMA(:,1,1);
    v2 = SIGMA(:,2,2);
    cv12 = SIGMA(:,1,2);
    v1h = SIGMA_HAT(:,1,1);
    v2h = SIGMA_HAT(:,2,2);
    cv12h = SIGMA_HAT(:,1,2);

    delmu = reshape(MU_HAT - MU,numel(v1),[]);
    t1 = v1.*v2-cv12.^2;

    T1 = log(abs(t1));
    T2 =v2.*(v1h+delmu(:,1).^2)./t1;
    T3= v1.*(v2h+delmu(:,2).^2)./t1;
    T4 = -2*cv12.*(cv12h+prod(delmu,2))./t1;

    nll = 0.5*(T1+T2+T3+T4);
    nll = mean(nll(:),'default','omitnan');


end