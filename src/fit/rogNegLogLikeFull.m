function nll = rogNegLogLikeFull(r,MU,SIGMA)
%% rogNegLogLikeFull Negative Log Likelihood for the pairwise ratio of Gaussians model
% NOTE: This is the full log likelihood, and not the approximation used for
% function fitting. See rogNegLogLikeParam for that function. This is
% computed directly using a bivariate normal pdf
%%
    [~,nS,nT] = size(r);

    nll = NaN(nS,nT);
    for iS = 1:nS
        nll(iS,:) = -log_bvnpdf(squeeze(r(:,iS,:)).',MU(:,iS).',squeeze(SIGMA(iS,:,:)));
    end
    nll(isinf(nll)) = NaN;
    nll = mean(nll(:),'default','omitnan');
end