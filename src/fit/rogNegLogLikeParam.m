function nll = rogNegLogLikeParam(pars,MU_HAT,SIGMA_HAT,contrasts,mu_eta,tuning) 
%% NLL_GAUSS_PARAM Computes the negative log likelihood for the Ratio of Gaussians model
% parametrized by contrast response function
%%
    arguments
        pars (1,19);
        MU_HAT (2,:,:);
        SIGMA_HAT (:,2,2);
        contrasts (1,:);
        mu_eta (2,1);
        tuning function_handle = @contrastResponse;
    end
    [MU,SIGMA] = mncovRogCrf(pars,contrasts,mu_eta,tuning);


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