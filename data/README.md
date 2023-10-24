Data files for the paper, Modeling the Diverse Effects of Divisive Normalization on Noise Correlations, available on [bioRxiv](https://www.biorxiv.org/content/10.1101/2022.06.08.495145v3.abstract)

Data files can be found on [Zenodo](https://zenodo.org/doi/10.5281/zenodo.8118187)
Please cite as:
```text
Modeling the Diverse Effects of Divisive Normalization on Noise Correlations
Oren Weiss, Hayley A. Bounds, Hillel Adesnik, Ruben Coen-Cagli
bioRxiv 2022.06.08.495145; doi: https://doi.org/10.1101/2022.06.08.495145
```
## Contents
| <center>Directory/File</center>                      |                                                                                                                                                                                                                                                                                                          <center>Contents </center> |
| :--------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------: |
| `data/raw`                                           | Raw calcium imaging recording sessions (files of the form `HBxx_x_xxx.mat`) and script `parseData.m`, a preprocessing script. This data contains calcium imaging data from primary visual cortex of mice viewing visual stimuli while performing a behavioral task. For more detail, see the publication and `data/raw/readme.txt`. |
| `data/processed`                                     |                                                                                                                                                                                                         Location of the output of `parseData.m` and of the parameters used for generating synthetic data (`params_best_JN2019.mat`) |
| `data/processed/PairwiseCalcium_NeuropilCoeff07.mat` |                                                                                                                                                                                           Responsive pairs of normalized fluorescence traces (with neuropil coefficient = 0.7) across all sessions. Main dataset used in this study |
| `data/simulations`                                   |                                                                                                                                                                        Stored results from simulation scripts (see `../scripts/` in the main [GitHub repository](https://github.com/oren-weiss/pairwiseRatioGaussians) for details) |
| `data/simulations/ApproxVTrue.mat`                   |                                                                                                                                                                                                                                     Simulation data used for Figure 1 in the paper (generated from `scripts/genDataApproxVsTrue.m`) |
| `data/simulations/NoiseCorrVsNormalization.mat`      |                                                                                                                                                                                                                        Simulation data used for Figures 2 in the paper (generated from `scripts/genDataNoiseCorrVsNormalization.m`) |
| `data/simulations/RhoFitVTrue.mat`                   |                                                                                                                                                                                                                              Simulation data used for Figure 3 in the paper (generated from `scripts/genDataSyntheticBootCIrand.m`) |
| `data/simulations/SingleTrialNorm.mat`               |                                                                                                                                                                                                                             Simulation data used for Figure 4 in the paper (generated from `scripts/genDataSingleTrialInference.m`) |
| `data/analysis`                                      |                                                                                                                                  Stored results from fitting the model to the processed calcium imaging data (`AnalysisCalciumCI.mat`). Corresponding analysis script is `scripts/analyzeCalciumDataCI.m`. Used for Figures 5 and 6 |

## File Tree
```text
data/
├── 📂 analysis
│   └── AnalysisCalciumCI.mat
├── 📂 processed
│   ├── PairwiseCalcium_NeuropilCoeff07.mat
│   └── params_best_JN2019.mat
├── 📂 raw
│   ├── HB67_2_209_190910_dataset_09012022.mat
│   ├── HB67_2_209_191014_dataset_09012022.mat
│   ├── HB67_2_209_191016_dataset_09012022.mat
│   ├── HB67_3_217_191008_dataset_09012022.mat
│   ├── HB67_3_217_191017_dataset_09012022.mat
│   ├── HB67_3_217_191025_dataset_09012022.mat
│   ├── HB67_3_217_191114_dataset_09012022.mat
│   ├── HB67_4_210_190926_dataset_09012022.mat
│   ├── HB67_4_210_191018_dataset_09012022.mat
│   ├── parseData.m
│   └── readme.txt
└── 📂 simulations
    ├── ApproxVTrue.mat
    ├── NoiseCorrVsNormalization.mat
    ├── PairwiseVsIndependent.mat
    ├── RhoFitVTrue.mat
    └── SingleTrialNorm.mat
```
