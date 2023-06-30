function NLL = mvnneglogpdf(X,MU,SIGMA,NormConst)
    arguments
        X (:,:) double;
        MU (:,:) double;
        SIGMA (:,:,:) double;
        NormConst logical = 0;
    end

    [nDists,nVars] = size(X);
    
    X0 = X-MU;
    
    mahaSqrt = zeros(nDists,nVars);
    logdetcov = zeros(nDists,1);

    if ismatrix(SIGMA)
        if nDists == 1
            SIGMA = repmat(SIGMA,1,1,2);
        else
            SIGMA = repmat(SIGMA,1,1,nDists);
        end
    end

    for iDist = 1:nDists
        U = cholcov(SIGMA(:,:,iDist),0);
        mahaSqrt(iDist,:) = X0(iDist,:) / U;
        logdetcov(iDist) = sum(log(diag(U)));
    end
    
    maha = sum(mahaSqrt.^2,2);

    NLL = 0.5*maha + logdetcov;

    if NormConst
        NLL = NLL + nVars*log(2*pi)/2;
    end
end

