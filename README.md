# Pairwise Ratio of Gaussians
MATLAB Toolbox and code to reproduce the results in our paper, Modeling the Diverse Effects of Divisive Normalization on Noise Correlations[^1].

# Quickstart
Run setup.m, which adds src/ to the MATLAB path for the current session. Then move to plots/ and try running makeFigureApproxVsTrue to produce Figure 1 in the paper. Functions in src/ have explanatory comments.

# Folder Layout
```
.
├── data
│   ├── processed
│   │   └── params_best_JN2019.mat
│   └── simulations
│       ├── ApproxVTrue.mat
│       ├── NoiseCorrVsNormalization.mat
│       ├── PairwiseVsIndependent.mat
│       ├── RhoFitVTrue.mat
│       └── SingleTrialNorm.mat
├── plots
│   ├── brewermap.m
│   ├── makeFigureApproxVsTrue.m
│   ├── makeFigureNoiseCorrVsNormalization.m
│   ├── makeFigureRhoFitVTrue.m
│   ├── makeFigureSingleTrial.m
│   └── nc_parametric.m
├── scripts
│   ├── genDataApproxVsTrue.m
│   ├── genDataNoiseCorrVsNormalization.m
│   ├── genDataSyntheticBootCIrand.m
│   └── gofVsTrials.m
├── setup.m
└── src
    ├── Contents.m
    ├── contrastResponse.m
    ├── covRogFull.m
    ├── covRogTaylorApprox.m
    ├── fit
    │   ├── fitMgCrfNFoldCv.m
    │   ├── fitRogCrf.m
    │   ├── fitRogCrfNFoldCv.m
    │   ├── log_bvnpdf.m
    │   ├── mgLogLikeParam.m
    │   ├── mgOptimBounds.m
    │   ├── modulatedPoissonLogLike.m
    │   ├── mvnneglogpdf.m
    │   ├── rogNegLogLikeFull.m
    │   ├── rogNegLogLikeParam.m
    │   └── rogOptimBounds.m
    ├── generateRogCrf.m
    ├── mncovEmpirical.m
    ├── mncovMgCrf.m
    ├── mncovRogCrf.m
    ├── normalizationSingleTrialInference.m
    ├── surroundSuppression.m
    ├── varPowerLaw.m
    └── varRogTaylorApprox.m
```
[^1]:Modeling the Diverse Effects of Divisive Normalization on Noise Correlations. Oren Weiss, Hayley A. Bounds, Hillel Adesnik, Ruben Coen-Cagli. bioRxiv 2022.06.08.495145; doi: https://doi.org/10.1101/2022.06.08.495145
