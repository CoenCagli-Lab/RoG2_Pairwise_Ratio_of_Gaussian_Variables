function [nll_mggauss,pfit,gof_mggauss] = fitMgCrfNFoldCv(nFold,spikes,cont,mu_eta,cov_eta,optional)
%% fitMgCrfNFoldCv N-fold Cross Validated Fit for the Modulated Gaussian model using the contrast response function parametrization. Fits both the pairwise and independent models.
%   INPUT
%       nFold - number of folds to cross-validate over (set = number of trials to get LOOCV)
%       spikes (2,contrast levels, number of trials) - pairwise activity
%       cont - contrast vector (0-100%)
%       mu_eta (2,1) - spontaneous mean firing rate
%       cov_eta (2,2) - covariance of additive noise
%       optional.optimopts - fmincon options (defualts sqp, 1000 function iterations)
%       optional.kind {"two-step","one-step"} - flag for whether to fit the MG model in two-steps (fitting the independent model, then freezing those parameters and only fitting the correlation parameters (i.e., the rho parameters) for the pairwise model) or in one-step (fitting the independent model and pairwise model seperately). Defaults to "two-step".
%       --------
%       optional.tuning (unused, testing only) - function handle for different tunings
%       optional.seed (unused, testing only)  - seed for training-test cross-validation splitting
%   OUTPUT
%       nll_mggauss - median cross validated negative log-likelihood for the modulated Gaussian model
%       pfit_pair (1,11) -  best fit parameters for the pairwise model NOTE: by setting pfit_pair([7:8 11]) = 0, this is the best fit parameters for the independent model
%       gof_mggauss - median cross validated goodness of fit normalized between an oracle model (negative log-likelihood using the training set mean and variance) and a null model (negative log likelihood using the training set mean and variance with no contrast dependence)
%%
%%
arguments
    nFold (1,1);
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
    optional.seed uint32 = [];
    optional.kind string {mustBeMember(optional.kind,["one-step","two-step"])} = "two-step";
end

nTrials = size(spikes,3);
if size(spikes,2) ~= size(cont(:))
    error('Number of Contrasts does not match!')
end

var_eta = diag(cov_eta);
if any(var_eta == 0)
    rho_eta = 0;
else
    rho_eta = cov_eta(1,2)/((prod(var_eta))^(0.5));
end

%% Cross-Validation Splits
% nFold = nTrials;
if isempty(optional.seed)
    rng('default');
else
    rng(optional.seed);
end
cvid = 1 + mod((1:nTrials)',nFold);
cvid = cvid(randperm(nTrials));

nllmg_pair = NaN(nFold,1);
nllmg_norm = NaN(nFold,1);
% nllind_norm = NaN(nFold,1);
%% Fitting...
% pfit = NaN(nFold,11);
for iFold=1:nFold

    %% CV Setting
    cv_set = (cvid~=iFold);
    r_train = spikes(:,:,cv_set);
    r_test = spikes(:,:,~cv_set);



    %% Empirical (hat) and Null mean and covariance
    [MU_HAT_TRAIN,SIGMA_HAT_TRAIN] = mncovEmpirical(r_train);
    switch optional.kind
        case "two-step"
            [x_pair,~] = fitMgCrfTwoStepInner( ...
                MU_HAT_TRAIN, ...
                SIGMA_HAT_TRAIN, ...
                cont, ...
                mu_eta, ...
                var_eta, ...
                rho_eta, ...
                optional.tuning, ...
                optional.optimopts);
        case "one-step"
            [x_pair,~] = fitMgCrfOneStepInner( ...
                MU_HAT_TRAIN, ...
                SIGMA_HAT_TRAIN, ...
                cont, ...
                mu_eta, ...
                var_eta, ...
                rho_eta, ...
                optional.tuning, ...
                optional.optimopts);
    end
    [MU_NB,SIGMA_NB] = mncovMgCrf(x_pair,cont,mu_eta,optional.tuning);
    %% Compute the goodness of fit scores
    tmp = -log_bvnpdf(reshape(r_test,2,[]).',mean(r_train,[2 3],'omitnan').',diag(var(reshape(r_train,2,[]).','omitnan')));
    tmp(isinf(tmp)) = NaN;
    nllnull = mean(tmp,'omitnan');
    %% ---------------------------------------------------------------------- %
    % These are set this way to reproduce the independent RoG oracle/null
    % structure and to allow for a fair comparison between independent and pairwise models
    % ----------------------------------------------------------------------- %
    SIGMA_HAT_TRAIN(:,1,2) = 0;
    SIGMA_HAT_TRAIN(:,2,1) = 0;

    nlloracle = rogNegLogLikeFull(r_test,MU_HAT_TRAIN,SIGMA_HAT_TRAIN);
    nllmg = rogNegLogLikeFull(r_test,MU_NB,SIGMA_NB);

    nllmg_norm(iFold) = (nllmg-nllnull)./(nlloracle-nllnull); %% currently unused
    nllmg_pair(iFold) = nllmg;
    % pfit(iFold,:) = x;


end
%% Fit the model on all trials to get the parameter estimates
%%% Can also take parameter estimates as mean across CV trials
switch optional.kind
    case "two-step"
        pfit = fitMgCrfTwoStep(spikes,cont,mu_eta,cov_eta,tuning=optional.tuning,optimopts=optional.optimopts);
    case "one-step"
        pfit = fitMgCrfOneStep(spikes,cont,mu_eta,cov_eta,tuning=optional.tuning,optimopts=optional.optimopts);
end
nll_mggauss = median(nllmg_pair,1,'omitnan');
gof_mggauss = median(nllmg_norm,1,'omitnan');
% pfit_pair = median(pfit,1,'omitnan');

end