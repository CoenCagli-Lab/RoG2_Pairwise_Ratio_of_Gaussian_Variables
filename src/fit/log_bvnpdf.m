function res = log_bvnpdf(X,MU,SIGMA)
% log_bvnpdf Log probability density for a bivariate Gaussian
% NOTE: log(2*pi) is removed since we always look at this as a relative value, so this constant factor gets subtracted off anyway
%%
    sig1 = SIGMA(1,1)^0.5;
    sig2 = SIGMA(2,2)^0.5;
    rho = SIGMA(1,2)/(sig1*sig2);
    CENT = X-MU;
    z = CENT(:,1).^2./(sig1^2)-2*rho*prod(CENT,2)./(sig1*sig2)+CENT(:,2).^2./(sig2^2);

%     res = -log(2*pi)-log(sig1)-log(sig2)-1/2*log(1-rho^2)-z/(2*(1-rho^2));
    res = -log(sig1)-log(sig2)-1/2*log(1-rho^2)-z/(2*(1-rho^2));
end