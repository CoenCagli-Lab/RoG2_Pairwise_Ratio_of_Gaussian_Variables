function var_out = varPowerLaw(alpha,beta,mu)
%% varPowerLaw Generates variance according to a power relation

    n = size(mu,2);
%         var_out = exp(log(alpha)+beta.*log(mu));
    alpha = repmat(alpha,1,n);
    beta = repmat(beta,1,n);
%     var_out = alpha.*(mu.^beta);
    var_out = exp(log(alpha)+beta.*log(mu));

    
end