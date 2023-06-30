function [rcov_mat,var_n,var_d] = covRogFull(alpha_n,beta_n,mu_n,alpha_d,beta_d,mu_d,rho)
%% covRogFull Creates a Covariance Matrix for a Ratio of Gaussians Model
%   Assumes variances for individual numerator/denominator variables are parametrized by a power law relationship to the mean
%   Covariance Between Vectors [N_1 D_1] and [N_2 D_2]
%   Format is (stimulus_number,4,4). If stimulus_number == 1, then the output is just (4,4) Covariance Matrix
%%
    arguments  
        alpha_n (2,1) {mustBeNumeric,mustBePositive};
        beta_n (2,1) {mustBeNumeric,mustBePositive} ;
        mu_n (2,:) {mustBeNumeric,mustBePositive};
        alpha_d (2,1) {mustBeNumeric,mustBePositive};
        beta_d (2,1) {mustBeNumeric,mustBePositive};
        mu_d (2,:) {mustBeNumeric,mustBePositive};
        rho (4,1) {mustBeNumeric};
    end

    var_d = varPowerLaw(alpha_d,beta_d,mu_d);
    var_n = varPowerLaw(alpha_n,beta_n,mu_n);

    stim_num = size(mu_n,2);
    rvar_n = rho(1).*(var_n(1,:).*var_n(2,:)).^(0.5);
    rvar_d = rho(2).*(var_d(1,:).*var_d(2,:)).^(0.5);
    rvar_n1d2 = rho(3).*(var_n(1,:).*var_d(2,:)).^(0.5);
    rvar_n2d1 = rho(4).*(var_n(2,:).*var_d(1,:)).^(0.5);
    rcov_mat = zeros(stim_num,4,4);
    for ind=1:stim_num
%       off_diag = [0    0    rvar_n(ind)    rvar_n1d2(ind);
%                   0    0    rvar_n2d1(ind) rvar_d(ind);
%                   0    0    0    0;
%                   0    0    0    0;];
%       rcov_mat(ind,:,:) = off_diag+off_diag.'+diag([var_n(1,ind) var_d(1,ind) var_n(2,ind) var_d(2,ind)]);

      rcov_mat(ind,:,:) = [var_n(1,ind)    0          rvar_n(ind)    rvar_n1d2(ind);
                           0           var_d(1,ind)   rvar_n2d1(ind) rvar_d(ind);
                           rvar_n(ind)    rvar_n2d1(ind)    var_n(2,ind)    0;
                            rvar_n1d2(ind)    rvar_d(ind)        0        var_d(2,ind)];
    end
end
