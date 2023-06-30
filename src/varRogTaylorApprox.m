function vd = varRogTaylorApprox(mu_n,var_n,mu_d,var_d,rho)
%% varRogTaylorApprox Variance for a Ratio of Univariate Gaussian Random Variables

    if nargin <5
        rho = 0;
    end
    t1 = mu_n.*mu_n;
    t2 = mu_d.*mu_d;
    vd = t1./t2.*(var_n./t1+var_d./t2-rho.*sqrt(var_n.*var_d)./(mu_n.*mu_d));
    
end