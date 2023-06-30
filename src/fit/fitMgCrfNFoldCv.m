function [nll_mggauss,pfit,gof_mggauss] = fitMgCrfNFoldCv(nFold,spikes,cont,mu_eta,cov_eta,tuning,opts)
%% fitMgCrfNFoldCv N-fold Cross Validated Fit for the Modulated Gaussian model using the contrast response function parametrization. Fits both the pairwise and independent models.
%   INPUT
%       nFold - number of folds to cross-validate over (set = number of trials to get LOOCV) %TODO 8/31/22 - Deal with cases where some trials are NaN (mismatched trials) - Possibly preprocess? 
%       spikes (2,contrast levels, number of trials) - pairwise activity
%       cont - contrast vector (0-100%)
%       mu_eta (2,1) - spontaneous mean firing rate
%       cov_eta (2,2) - covariance of additive noise
%       tuning - function handle, unused
%       opts - fmincon options
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
        tuning function_handle = @contrastResponse;
        opts struct = optimoptions(@fmincon,...
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



    % nFold = nTrials;
    cvid = 1 + mod((1:nTrials)',nFold);
    cvid = cvid(randperm(nTrials));

    nllmg_pair = NaN(nFold,1);
    nllmg_norm = NaN(nFold,1);
    nllind_norm = NaN(nFold,1);
%% Can take the median across CV folds to get parameter estimates
   pfit = NaN(nFold,11);



    for iFold=1:nFold

        %% CV Setting
        cv_set = (cvid~=iFold);
        r_train = spikes(:,:,cv_set);
        r_test = spikes(:,:,~cv_set);
       
       

        %% Empirical (hat) and Null mean and covariance
        [MU_HAT_TRAIN,SIGMA_HAT_TRAIN] = mncovEmpirical(r_train);
        fun_temp = @(x)(mgLogLikeParam(x,MU_HAT_TRAIN,SIGMA_HAT_TRAIN,cont,mu_eta));

        [pars_0,lb_pars,ub_pars] = mgOptimBounds(max(MU_HAT_TRAIN,2),var_eta);
%% ---------------------------------------------------------------------- % 
% % these lines fix rho_eta, var_eta to be the empirically recorded values,
% % and not parameters to be optimized. This makes a quantitative 
% % difference in the quality of fits but not a qualitative difference in 
% % the results.
%% ---------------------------------------------------------------------- % 
%       pars_0(10:12) = [reshape(var_eta,1,2),rho_eta];
%       lb_pars(10:12) = [reshape(var_eta,1,2),rho_eta];
%       ub_pars(10:12) = [reshape(var_eta,1,2),rho_eta];
%%
%% If negative log-likelihood returns NaN, computation will error out
% % Continue to next loop
        if isnan(fun_temp(pars_0))
            continue
        end
%% Optimize the parameters that are independent of correlation structrue (all besides rho parameters)
        x_ind = fmincon(fun_temp,pars_0,[],[],[],[],lb_pars,ub_pars,[],opts);
        % % Fix these parameters, change the range for the rho parameters
        lb_pars = x_ind;
        ub_pars = x_ind;
        lb_pars(7:8) = [-1 -1];
        ub_pars(7:8) = [1 1];
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
%% Optimize full model
        x = fmincon(fun_temp,x_ind,[],[],[],[],lb_pars,ub_pars,[],opts);
        [MU_NB,SIGMA_NB] = mncovMgCrf(x,cont,mu_eta);
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
        pfit(iFold,:) = x;


    end

    nll_mggauss = median(nllmg_pair,1,'omitnan');
    gof_mggauss = median(nllmg_norm,1,'omitnan');
    pfit = median(pfit,1,'omitnan');

end
% function [pars_0,lb,ub] = mgOptimBounds(spikes,var_eta,rho_eta)
%         rm_0 = reshape(max(mean(spikes,3,'default','omitnan'),[],2),1,[]);
%         % % % par 17-18: var_et
%         var_et_0 = var_eta(:).';
% 
%         pars_0 = [...
%             rm_0,...
%             25,...
%             25,...
%             1,...
%             1,...
%             0,...
%             0,...
%             var_et_0,...
%             rho_eta];
%         lb = [...
%             0.5*rm_0,...
%             1,...
%             1,...
%             0.01,...
%             0.01,...
%             -1,...
%             -1,...
%             0.1*var_et_0,...
%             -1];
%         ub = [...
%             2*rm_0,...
%             100,...
%             100,...
%             20,...
%             20,...
%             1,...
%             1,...
%             10*var_et_0,...
%             1];
% end
% 
% function [MU_R,SIG_R] = mncovMgCrf(params,cont,mu_eta)
%     [mu_n,mu_d] = contrastResponse(params(1:2),params(3:4),cont);
% 
%     MU_R = mu_n./mu_d+reshape(mu_eta,2,[]);
% 
%     g = reshape(params(5:6),2,[]);
%     v_eta = reshape(params(9:10),2,[]);
%     sigma_R2 = MU_R+g.^2.*MU_R.^2+v_eta;
% 
%     cv = params(7)*prod(MU_R,1).^(0.5)+params(8)*prod(g.*MU_R,1)+params(11)*prod(v_eta.^(0.5),1);
% 
%     nStims = numel(cont);
%     SIG_R = NaN(nStims,2,2);
%     for iStim=1:nStims
%         SIG_R(iStim,:,:) = [   sigma_R2(1,iStim) cv(iStim);...
%                     cv(iStim)   sigma_R2(2,iStim);];
%     end
% end