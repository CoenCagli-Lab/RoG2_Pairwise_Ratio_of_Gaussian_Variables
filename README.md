# Pairwise Ratio of Gaussians
MATLAB Toolbox and code to reproduce the results in our paper, Modeling the Diverse Effects of Divisive Normalization on Noise Correlations[^1].

# Quickstart
Run setup.m, which adds src/ to the MATLAB path for the current session. Then move to plots/ and try running makeFigureApproxVsTrue to produce Figure 1 in the paper. Functions in src/ have explanatory comments.

# Requirements
- Parallel Computing Toolbox: for generation of synthetic data and analysis of calcium imaging data
- Optimization Toolbox: used to fit the model to data using `fmincon` (note: it is possible to replace `fmincon` with another nonlinear optimization function that doesn't require this toolbox, e.g. [BADS](https://github.com/acerbilab/bads) or [minConf](https://www.cs.ubc.ca/~schmidtm/Software/minConf.html)

# Folder Descriptions
 - `src`: the main functions involved in using the model, including synthetic generation of data and optimization. See `src/Contents.m` for an overview, and the comments in the files for descriptions on functionality.
 - `scripts`: generate simulation data and to analyze neuronal data (this file is not yet available). These files take a long time to run, so recommended to run in batch mode on a server with the Parallel Computing Toolbox.
 - `data`: contains four subfolders
    - `raw (NOT YET AVAILABLE)`: raw calcium imaging data with preprocessing script (see the Methods section for more details)
    - `processed`: where calcium imaging data is saved after being processed (NOT YET AVAILABLE) and contains `params_best_JN2019.mat`, which contains a cell structure used in some of the simulations in `scripts` (see Methods for full description)
    - `simulations`: where simulation data generated `scripts` are saved. See `simulations/README.md` for a link to simulations used in the paper.
    - `analysis (NOT YET AVAILABLE)`: analysis of calcium imaging data (i.e. parameter fits)
- `plots`: plotting functions to generate Figures 1-4 in the paper using the simulation data. Scripts to generate Figures 5,6 and some of the Supplemental Figures are not yet available.

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
