function [nll] = modulatedPoissonLogLike(params,response,contrasts,mu_eta)
    [~,nS,nT] = size(response);
    response(response < -1) = NaN;
    nC = numel(contrasts);
    [mu_n,mu_d] = contrastResponse(params(1:2),params(3:4),contrasts);
    
    mu_r = mu_n./mu_d + reshape(mu_eta,2,1);
    sigg2 = reshape(params(5:6).^2,2,1);
%     sigg2 = repmat(sigg2(:),1,nC);
%     r = ones(size(mu_r)) ./ sigg2;
    rparam = 1 ./ sigg2;
    sparam = sigg2 .* mu_r;
    gammaR = gammaln(rparam);
    nll = NaN(2,nS,nT);
%     s = sigg2 .* mu_r;
    for iS = 1:nS
        sLoop = sparam(:,iS);
        Lplus = log1p(sLoop);
        R = squeeze(response(:,iS,:));
        nll(:,iS,:) = rparam .* Lplus + R .* Lplus  - R .* log(sLoop) + gammaR + gammaln(R + 1) - gammaln(R + rparam);
    end

    nll = mean(nll(:),'omitnan');
end