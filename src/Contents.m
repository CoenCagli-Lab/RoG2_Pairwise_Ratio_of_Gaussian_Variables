% SRC
%
% Files
%   contrastResponse                  - contrastResponse Tuning curve for standard contrast-response hyperbolic ratio function
%   covRogFull                        - covRogFull Creates a Covariance Matrix for a Ratio of Gaussians Model
%   covRogTaylorApprox                - covRogTaylorApprox Creates Covariance Matrix for a Ratio of Gaussians Model using the Taylor Series Approximation
%   generateRogCrf                    - generateRogCrf produces spike counts/rates for the pairwise Ratio of Gaussians model assuming a contrast response parametrization
%   mncovEmpirical                    - mncovEmpirical Computes the empirical mean and covariance matrices for pairwise neural response
%   mncovMgCrf                        - mncovMgCrf Returns the Mean/Covariance Approximations for the modulated Gaussians Model using the contrast response parametrization
%   mncovRogCrf                       - mncovRogCrf Returns the Mean/Covariance Approximations for the pairwise Ratio of Gaussians Model using the contrast response parametrization
%   normalizationSingleTrialInference - normalizationSingleTrialInference Finds the Maximum a Priori estimates for the normalization signal in a Ratio of Gaussians model
%   surroundSuppression               - surroundSuppresion Tuning curve for a pseudo-surround suppression model
%   varPowerLaw                       - varPowerLaw Generates variance according to a power relation
%   varRogTaylorApprox                - varRogTaylorApprox Variance for a Ratio of Univariate Gaussian Random Variables

% FIT
%
% Files
%   fitMgCrfNFoldCv         - fitMgCrfNFoldCv N-fold Cross Validated Fit for the Modulated Gaussian model using the contrast response function parametrization. Fits both the pairwise and independent models.
%   fitRogCrf               - 
%   fitRogCrfNFoldCv        - fitRogCrfNFoldCv N-fold Cross Validated Fit for the Ratio of Gaussians model using the contrast response function parametrization. Fits both the pairwise and independent models.
%   log_bvnpdf              - log_bvnpdf Log probability density for a bivariate Gaussian
%   mgLogLikeParam          - mgLogLikeParam Computes the parametrized negative Log-likelihood of a bivarate Gaussian model with Poisson-Gamma inspired moments.
%   mgOptimBounds           - 
%   modulatedPoissonLogLike - 
%   mvnneglogpdf            - 
%   rogNegLogLikeFull       - rogNegLogLikeFull Negative Log Likelihood for the pairwise ratio of Gaussians model
%   rogNegLogLikeParam      - NLL_GAUSS_PARAM Computes the negative log likelihood for the Ratio of Gaussians model
%   rogOptimBounds          - rogOptimBounds Produces the optimization bounds/starting points for the pairwise Ratio of Gaussians parametrized by the contrast response function
