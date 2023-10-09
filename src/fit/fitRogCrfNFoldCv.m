function [nll_pair,nll_ind,pfit_pair,gof_pair,gof_ind,gof_diff] = fitRogCrfNFoldCv(nFold,spikes,cont,mu_eta,cov_eta,tuning,opts)
%% fitRogCrfNFoldCv N-fold Cross Validated Fit for the Ratio of Gaussians model using the contrast response function parametrization. Fits both the pairwise and independent models.
%   INPUT
%       nFold - number of folds to cross-validate over (set = number of trials to get LOOCV)
%       spikes (2,contrast levels, number of trials) - pairwise activity
%       cont - contrast vector (0-100%)
%       mu_eta (2,1) - spontaneous mean firing rate
%       cov_eta (2,2) - covariance of additive noise
%       tuning - function handle, unused
%       opts - fmincon options
%   OUTPUT
%       nll_pair, nll_ind - median cross validated negative log-likelihood for the pairwise/independent model
%       pfit_pair (1,19) -  best fit parameters for the pairwise model NOTE: by setting pfit_pair([13:16 19]) = 0, this is the best fit parameters for the independent model
%       gof_pair, gof_ind - median cross validated goodness of fit normalized between an oracle model (negative log-likelihood using the training set mean and variance) and a null model (negative log likelihood using the training set mean and variance with no contrast dependence)
%       gof_diff - median difference in goodness of fit. Largely unused as median(gofPair - gofInd) is approximately median(gofPair)-median(gofInd)
%
%%
    arguments
        nFold (1,1);
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


    %% Set optimization parameters

    % nFold = nTrials;
    cvid = 1 + mod((1:nTrials)',nFold);
    cvid = cvid(randperm(nTrials));

    nllrog_pair = NaN(nFold,1);
    nllrog_ind = NaN(nFold,1);
    nllpair_norm = NaN(nFold,1);
    nllind_norm = NaN(nFold,1);
%% Can take the median across CV folds to get parameter estimates
  pfit = NaN(nFold,19);


    for iFold=1:nFold

        %% CV Setting
        cv_set = (cvid~=iFold);
        r_train = spikes(:,:,cv_set);
        r_test = spikes(:,:,~cv_set);


        %% Empirical (hat) and Null mean and covariance
        [MU_HAT_TRAIN,SIGMA_HAT_TRAIN] = mncovEmpirical(r_train);
        fun_temp = @(x)(rogNegLogLikeParam(x,MU_HAT_TRAIN,SIGMA_HAT_TRAIN,cont,mu_eta,tuning));

        [pars_0,lb_pars,ub_pars] = rogOptimBounds(max(MU_HAT_TRAIN,[],2),var_eta);
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
            continue
        end
%% Optimize the parameters that are independent of correlation structrue (all besides rho parameters)
        x_ind = fmincon(fun_temp,pars_0,[],[],[],[],lb_pars,ub_pars,[],opts);
        % % Fix these parameters, change the range for the rho parameters
        lb_pars([1:12 17:18]) = x_ind([1:12 17:18]);
        ub_pars([1:12 17:18]) = x_ind([1:12 17:18]);
%% ---------------------------------------------------------------------- % 
% %  We only considered the case where rho_n1d2=rho_n2d1=0
% %  Change these lines to look at the more general case
% % ----------------------------------------------------------------------- %
%        lb_pars(13:16) = [-1 -1 -1 -1];
%        ub_pars(13:16) = [1 1 1 1]; 

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

        x_pair = fmincon(fun_temp,x_ind,[],[],[],[],lb_pars,ub_pars,[],opts);
%% Compute the goodness of fit scores
        [MU_PAIR,SIGMA_PAIR] = mncovRogCrf(x_pair,cont,mu_eta,tuning);
        [MU_IND,SIGMA_IND] = mncovRogCrf(x_ind,cont,mu_eta,tuning);

        tmp = -log_bvnpdf(reshape(r_test,2,[]).',mean(r_train,[2 3],'omitnan').',diag(var(reshape(r_train,2,[]),0,2,'omitnan')));
        tmp(isinf(tmp)) = NaN;
        nllnull = mean(tmp,'omitnan');
%% ---------------------------------------------------------------------- % 
% These are set this way to reproduce the independent RoG oracle/null 
% structure and to allow for a fair comparison between independent and pairwise models
% ----------------------------------------------------------------------- %
        SIGMA_HAT_TRAIN(:,1,2) = 0;
        SIGMA_HAT_TRAIN(:,2,1) = 0;
%%
        nlloracle = rogNegLogLikeFull(r_test,MU_HAT_TRAIN,SIGMA_HAT_TRAIN);
        nllpair = rogNegLogLikeFull(r_test,MU_PAIR,SIGMA_PAIR);
        nllind = rogNegLogLikeFull(r_test,MU_IND,SIGMA_IND);


        nllpair_norm(iFold) = (nllpair-nllnull)./(nlloracle-nllnull);
        nllrog_pair(iFold) = nllpair;
        nllrog_ind(iFold) = nllind;
        nllind_norm(iFold) = (nllind-nllnull)./(nlloracle-nllnull);
      % pfit(iFold,:) = x_pair;

    end
%% Fit the model on all trials to get the parameter estimates
% % Can also take parameter estimates as mean across CV trials
    pfit_pair = fitRogCrf(spikes,cont,mu_eta,cov_eta,tuning);
%   pfit_pair = median(pfit_pair,1,'omitnan');
%% Compute the median scores (goodness of fit, negative log likelihood)
    nll_pair = median(nllrog_pair,1,'omitnan');
    nll_ind = median(nllrog_ind,1,'omitnan');
    gof_pair = median(nllpair_norm,1,'omitnan');
    gof_ind = median(nllind_norm,1,'omitnan');
    gof_diff = median(nllpair_norm-nllind_norm,1,'omitnan');

end
