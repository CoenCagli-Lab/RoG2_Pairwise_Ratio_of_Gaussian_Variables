function [] = genDataSyntheticBootCIrand(nParallel,nEXP,nBoots,nTrials,Contrasts,p)
%   Generate synthetic data to investigate parametric inference in the model with respect to bootstrapping
%   Must be done in parallel in order to be completed in a reasonable timespan.
%%  INPUT
%       nParallel - number of cores (must be > 1)
%       nEXP - number of simulated pairs (default/max 11628)
%       nBoots - number of bootstrap resamples (default 1e3)
%       nTrials - number of trials for simulations (default 1e3)
%       Contrasts - set of contrasts (1-100%) (default 6.25,12.5,25,50,100)
%%% NAME-VALUE INPUTS
%       name - appended to the basename (which is YYYYMMDD_RhoFitVTrue)
%       seed - random seed for sampling
%       tuning - tuning function, unused
%       nfold - number of cross-validation folds (default 10, don't go above this or runtime will increase)
%%  OUTPUT File Contents
%       allCorr (nEXP, # of contrasts) - empirical noise correlations (computed from simulated activity)
%       trueParameters (nEXP, 19) - true parameters generating the ratio of Gaussians
%       fitParameters (nEXP, 19) - fit parameters for the ratio of Gaussians
%       fitParametersBoot (nEXP, nBoots, 19) - fit parameters for bootstrap resamples to the ratio of Gaussians
%       nllPairwise (nEXP) - cross-validated negative log-likelihood for the pairwise Ratio of Gaussians
%       gofPairwise (nEXP) - cross-validated goodness of fit for the pairwise ratio of Gaussians
%       gofIndependent (nEXP) - cross-validated goodness of fit for the independent Ratio of Gaussaisn
%       genParams - convenience struct containing input parameters
%%
arguments
    nParallel (1,1) {mustBeInteger,mustBePositive};
    nEXP (1,1) {mustBeInteger} = 0;
    nBoots (1,1) {mustBeInteger,mustBePositive} = 1e3;
    nTrials (1,1) {mustBePositive,mustBeInteger}= 1e3;
    Contrasts (1,:) {mustBeLessThanOrEqual(Contrasts, 100)} = [6.25 12.5 25 50 100];
    p.name string = "";  %names will be appended at the end of filename
    p.seed (1,1) {mustBeInteger,mustBePositive} = 96309;
    p.tuning function_handle = @contrastResponse;
    p.nfold (1,1) {mustBeInteger,mustBePositive} = 10;
end
%%
baseFolder = fileparts(cd);
fprintf("Current Folder: %s\nBase Folder: %s\n", pwd,baseFolder);
name = p.name;
seed = p.seed;
tuning = p.tuning;
nfold = p.nfold;
rng(seed); % Reproducibility
try
    fprintf([repelem('-',29) 'Loading crf parameters' repelem('-',29) newline])
    addpath(genpath(fullfile(baseFolder,'src')));
    load('../data/processed/params_best_JN2019.mat','parampairs');
catch
    error('Error loading data/adding to path. Check path is correct')
end

genParams = struct('tuning',tuning,'nParallel',nParallel,'nTrials',nTrials,'Contrasts',Contrasts,'nBoots',nBoots,'time',[]); %#ok<*PFOUS>

if nEXP == 0
    nEXP = size(parampairs,1);
elseif nEXP > 11628
    error("nEXP must be less than 11628");
else
    parampairs = parampairs(randsample(size(parampairs,1),nEXP),:);
end
nContrasts = numel(Contrasts);
baseName = "RhoFitVTrue";
nIterations = nEXP;
%% Setup Parallel Pool
fprintf([repelem('-',28) 'Setting up Parallel Pool' repelem('-',28) newline])
if nParallel > 1
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        poolobj = parpool('local',nParallel);
        addAttachedFiles(poolobj,regexp(genpath(fullfile(baseFolder,'src')),':', 'split'));
    end
else
    error('Must set nParallel > 1')
end
paramsFixed = parallel.pool.Constant(parampairs);
Cont = parallel.pool.Constant(Contrasts);
nFold = parallel.pool.Constant(nfold);
%% Figure Rho Identifiable
allCorr = NaN(nIterations,nContrasts);
trueParameters = NaN(nIterations,19);
gofPairwise = NaN(nIterations,1);
gofIndependent = NaN(nIterations,1);
fitParameters = NaN(nIterations,19);
fitParametersBoot = NaN(nIterations,nBoots,19);
nllPairwise = NaN(nIterations,1);
objSave = {'genParams',...
    'fitParameters',...
    'trueParameters',...
    'fitParametersBoot',...
    'allCorr',...
    'gofPairwise',...
    'gofIndependent',...
    'nllPairwise'};
fprintf('Saving Variables:');fprintf(' %s ',objSave{:}); fprintf('\n')
%% Setup File Names
Today = string(datetime('now','Format','yMMdd'));
saveprefix = strjoin([Today,baseName],"_");
if name ~= ""
    name = strjoin([saveprefix,name],"_");
else
    name = saveprefix;
end
name = name + ".mat";
fprintf('Saving in %s\n',name);
fprintf([repelem('-',80) newline]);
fprintf("Beginning Fitting\nTuning function: %s\n# of Iterations (Pairs): %d\n# of Bootstraps: %d:\n# of Trials: %d\n",functions(tuning).function,nIterations,nBoots,nTrials);
fitSTART = datetime('now','Format','eeee, MMMM d, yyyy h:mm a');
fprintf('Start time: %s\n',string(fitSTART));
allt = tic;

parfor iIteration=1:nIterations
    % for iIteration=1:nIterations
    %         [iEXP,iRhoN,iRhoD] = ind2sub(szLoop.Value,iIteration);
    %         rN = RhoNloop.Value(iRhoN);
    %         rD = RhoDloop.Value(iRhoD);
    pairLoop = paramsFixed.Value(iIteration,:);
    % pairLoop = parampairs(iIteration,:);
    C = Cont.Value;
    % C = Contrasts;
    %
    % pairLoop = parampairs(iIteration,:);
    % C = Contrasts;
    rN = -0.9 + 1.8*rand();
    rD = -0.9 + 1.8*rand();
    %% Generate Data from the Ratio of Gaussians
    [spikes,~,parameters] = generateRogCrf(nTrials,...
        nContrasts,...
        pairLoop{1:12},...
        'cont',C,...
        'rho_n',rN,...
        'rho_d',rD, ...
        'tuning',tuning);
    [~,~,allCorr(iIteration,:)] = mncovEmpirical(spikes);
    %% Fit Ratio of Gaussians to synthetic data
    [nllPairwise(iIteration),~,params,gofPairwise(iIteration),gofIndependent(iIteration),~] = ...
        fitRogCrfNFoldCv(nFold.Value,spikes,C,parameters.mu_eta,parameters.sigma_eta,"tuning",tuning,"kind","two-step");

    fitParameters(iIteration,:) = params;
    trueParameters(iIteration,:) = parameters.pars;
    %% Bootstrap Fitting
    for iBoot=1:nBoots
        bootResponse = spikes(:,:,randsample(nTrials,nTrials,true));
        paramfitBoot = ...
            fitRogCrfTwoStep(bootResponse,C,parameters.mu_eta,parameters.sigma_eta,"tuning",tuning);
        fitParametersBoot(iIteration,iBoot,:) = paramfitBoot;
    end
end
fitEND = datetime('now','Format','eeee, MMMM d, yyyy h:mm a');
fprintf('Finshed!\nEnd time: %s\n',string(fitEND));
fitDURATION = fitEND - fitSTART;
allt = toc(allt);
genParams.time = allt;
fprintf('Elapsed Time = %s\n',fitDURATION);
fprintf("\nSaving as %s in ../data/simulations...\n",name)
try
    save(fullfile(baseFolder,"data","simulations",name),objSave{:},'-v7.3');
    fprintf("Success!\n");
catch
    save("fail.mat");
    fprintf("Failed to Save!\n")
end
fprintf([repelem('-',80) newline]);
end