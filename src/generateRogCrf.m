function [spikes,D,parameters] = generateRogCrf(num_trials,num_stims,parameters)
%% generateRogCrf produces spike counts/rates for the pairwise Ratio of Gaussians model assuming a contrast response parametrization
%   INPUTS
%       num_trials - number of trials (default = 1000)
%       num_stims - number of contrast levels, 0-100% (default = 10)
%       parameters (name-value pairs) -
%           Rmax (2,1) - contrast response max rate (default U[5,50])
%           eps (2,1) - contrast response semi-saturation (default U[10,100])
%           alpha_n (2,1) - power-law variance proportion for numerator (default U[0.1,1])
%           beta_n (2,1) - power-law variance exponent for numerator (default U[1,2])
%           alpha_d (2,1) - power-law variance proportion for denominator (default U[0.1,1])
%           beta_d (2,1) - power-law variance exponent for denominator (default U[1,2])
%           rho_n - correlation between numerators (default U[-0.5,0.5])
%           rho_d - correlation between denominators (default U[-0.5,0.5])
%           rho_n1d2 - correlation between N1 and D2 (default 0)
%           rho_n2d1 - correlation between N2 and D1 (default 0)
%           rho (4,1) - vector format of the rho params (default uninitialized, overwrites previous rho params if given)
%           var_eta (2,1) spontaneous variance for the pair (default [0;0])
%           rho_eta - spontaneous correlation (default 0)
%           mu_eta (2,1) - spontaneous mean (default 0)
%           cont (1,num_stims) - contrast vector (overwrites num_stims if given)
%   OUTPUTS
%       spikes (2,num_stims,num_trials) - of spike counts/rates
%       D (2,num_stims,num_trials) - of the true denominators used to generate spikes
%       parameters (struct) -  same as input parameters with additional vector
%           pars (1,19) - Vector of the input parameters in proper format for use by other functions
%                           [Rmax(2) eps(2) alphan1 betan1 alphad1 betad1 alphan2 betan2 alphasd2 betad2 rhon rhod rhon1d2 rhon2d1 var_eta(2) rho_eta]
%%
    arguments
        % % Base Parameters
        num_trials (1,1) {mustBePositive} = 1e3;
        num_stims (1,1) {mustBeInteger,mustBePositive} = 10;

        % % Mean Parameters (for contrast response function)
        parameters.Rmax (2,1) = 5+45*rand(2,1);
        parameters.eps (2,1) =  10+90*rand(2,1); % = internal_lognorm(29,16)

        % % Variance Parameters
        parameters.alpha_n (2,1) = max(rand(2,1),0.1);
        parameters.beta_n (2,1) = 1+rand(2,1);
        parameters.alpha_d (2,1) = max(rand(2,1),0.1);
        parameters.beta_d (2,1) = 1+rand(2,1);

        % % Rho Parameters
        parameters.rho_n (1,1) = -0.5+rand(1);
        parameters.rho_d (1,1) = -0.5+rand(1);
        parameters.rho_n1d2 (1,1)  = 0;
        parameters.rho_n2d1 (1,1)  = 0;
        parameters.rho (4,1) {mustBeNumeric};

        % % Noise Parameters
        parameters.var_eta (2,1) = zeros(2,1);
        parameters.rho_eta (1,1) = 0;
        parameters.mu_eta (2,1) = zeros(2,1);

        % % Function Handle
        parameters.tuning function_handle = @contrastResponse;
        
        % % Contrast Space
        parameters.cont (1,:);

        % % Override Parameters
        parameters.fit_parameters (1,19);
    end
%% Preprocessing
    % % If 'cont' is provided, this overwrites num_stims to be the size of 'cont'
    if isfield(parameters,'cont')
        num_stims = size(parameters.cont,2);
    end

    % % If 'cont' is not provided, then 'cont' is set using 'num_stims'
    if ~isfield(parameters,'cont')
        parameters.cont = linspace(10,100,num_stims);
    end
    cont = parameters.cont;

    % % If 'fit_parameters' is provided, this will override all other
    % parameter settings (that are optimized)
    if isfield(parameters,'fit_parameters')
        parameters.Rmax = parameters.fit_parameters(1:2);
        parameters.eps = parameters.fit_parameters(3:4);
        parameters.alpha_n = parameters.fit_parameters([5 9]);
        parameters.beta_n  = parameters.fit_parameters([6 10]);
        parameters.alpha_d = parameters.fit_parameters([7 11]);
        parameters.beta_d = parameters.fit_parameters([8 12]);
        parameters.rho = parameters.fit_parameters(13:16);
        parameters.var_eta = parameters.fit_parameters(17:18);
        parameters.rho_eta = parameters.fit_parameters(19);
    end

    % % If 'rho' is provided, this overwrites the 'rho_X' parameters (besides 'rho_eta')
    if ~isfield(parameters,'rho')
        parameters.rho =[parameters.rho_n parameters.rho_d parameters.rho_n1d2 parameters.rho_n2d1];
    end

%% Mean and Covariance Setting
%     [parameters.mu_n,parameters.mu_d] = contrastResponse(parameters.Rmax,parameters.eps,cont);
    [parameters.mu_n,parameters.mu_d] = parameters.tuning(parameters.Rmax,parameters.eps,cont);
%     parameters.mu_n = mu_n;
%     parameters.mu_d = mu_d;

    [cov_ND,parameters.var_n,parameters.var_d] = covRogFull(...
            parameters.alpha_n,...
            parameters.beta_n,...
            parameters.mu_n,...
            parameters.alpha_d,...
            parameters.beta_d,...
            parameters.mu_d,...
            parameters.rho);

%% Noise Setting
    c_eta = parameters.rho_eta*sqrt(parameters.var_eta(1)*parameters.var_eta(2));
    sigma_eta = diag(parameters.var_eta)+diag(c_eta,1)+diag(c_eta,-1);
    parameters.sigma_eta = sigma_eta;
%%
%% Spike Generation
%%
    spikes = zeros(2,num_stims,num_trials);
    ND = zeros(4,num_stims,num_trials);
    for k=1:num_stims
            ND(:,k,:) = [parameters.mu_n(1,k); parameters.mu_d(1,k); parameters.mu_n(2,k);parameters.mu_d(2,k)].*ones(1,num_trials)...
                            + internal_chol(squeeze(cov_ND(k,:,:))).'*randn(4,num_trials);


        ND(:,k,:) = reshape(...
                    squeeze(ND(:,k,:)).*[ones(1,num_trials); (reshape(ND(2,k,:),1,[])>0); ones(1,num_trials); (reshape(ND(4,k,:),1,[])>0)],...
                    4,...
                    1,...
                    []);
        ND(ND==0) = NaN;

        spikes(:,k,:) = [squeeze(ND(1,k,:)./ND(2,k,:)).'; squeeze(ND(3,k,:)./ND(4,k,:)).']...
                            +mvnrnd(parameters.mu_eta,sigma_eta,num_trials).';
    end
    D = ND([2 4],:,:);

        parameters.pars = [parameters.Rmax(:);...
                        parameters.eps(:);...
                        parameters.alpha_n(1);...
                        parameters.beta_n(1);...
                        parameters.alpha_d(1);...
                        parameters.beta_d(1);...
                        parameters.alpha_n(2);...
                        parameters.beta_n(2);...
                        parameters.alpha_d(2);...
                        parameters.beta_d(2);...
                        parameters.rho(:);...
                        parameters.var_eta(:);...
                        parameters.rho_eta].';
end
function [U,isPosDef] = internal_chol(A)
%CHOLPROJ  Projected Cholesky factorization.
% cholproj(A) returns an upper triangular matrix U so that U'*U = A,
% provided A is symmetric positive semidefinite (sps).
%
% If A is not sps, then U will approximately satisfy U'*U = A.   
% This is useful when dealing with matrices that are affected
% by roundoff errors.  By multiplying U'*U you effectively round A to the 
% nearest sps matrix.
%
% [U,isPosDef] = cholproj(A) also returns whether A is positive definite.
% % %
% Taken from Lightspeed Toolbox (c) 2017 Tom Minka MIT License
    [rowsA,colsA] = size(A);
    U = zeros(rowsA,colsA);
    isPosDef = 1;
    for i = 1:colsA
        for j = i:rowsA
            k = 1:(i-1);
            s = A(i,j) - U(k,i)'*U(k,j);
            if i == j
                if s <= 0
                    isPosDef = 0;
                    U(i,i) = 0;
                else
                    U(i,i) = sqrt(s);
                end
            else
                if U(i,i) > 0
                    U(i,j) = s / U(i,i);
                else
                    U(i,j) = 0;
                end
            end
        end
    end
end