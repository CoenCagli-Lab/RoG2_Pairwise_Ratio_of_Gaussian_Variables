function [] = genDataSingleTrialInference(nParallel,nEXP,nTrials,Contrasts,RhoN,RhoD,p)
% genDataSingleTrialInference Generates synthetic data from the pairwise Ratio of Gaussians model using model fits to previously recorded data to measure single-trial normalization.
% 	Produces large matrices comparing the mean squared error and correlation
% 		between the true normalization and the MAP estimate for pairwise and independent
% 		Ratio of Gausians
%%  INPUT
%       nParallel - number of parallel cores (must be > 1)
%       nEXP - number of simulated pairs (default/max 11628)
%       nTrials - number of synthetic trials (default 1e3)
%       Contrasts - contrast levels used for synthesis (1-100%) (default 6.25, 12.5, 25, 50, 100)
%       RhoN, RhoD - space of rho_N, rho_D values used to generate data (default RhoN = -0.5:0.05:0.5, RhoD = -0.5:0.05:0.5)
%       name - name appended to the default name (
%%% NAME-VALUE INPUTS
%       name - appended to the basename (which is YYYYMMDD_RhoFitVTrue)
%       seed - random seed for sampling
%       tuning - tuning function, unused
%%  OUTPUT File Contents
%       genParams - saves input parameters (convenience structure)
%       CorrNormSignalPairwise - Corr between pairwise estimate of normalization and ground truth
% 	    CorrNormSignalIndependent - Corr between independent estimate of normalization and ground truth
% 		MSENormSignalPairwise - mean squared error between pairwise estimate of normalization and ground truth
% 		MSENormSignalIndependent - mean squared error between pairwise estimate of normalization and ground truth
% 		MSENormSignalPairwiseRand
% 		MSENormSignalIndependentRand
% 		CorrNormSignalPairwiseRand
%		CorrNormSignalIndependentRand
% Files are saved in data/simulations
%%
arguments
    nParallel (1,1) {mustBeInteger,mustBePositive};
    nEXP (1,1) {mustBeInteger,mustBePositive} = 0;
    nTrials (1,1) {mustBePositive,mustBeInteger}= 1e3;
    Contrasts (1,:) {mustBeLessThanOrEqual(Contrasts, 100)} = [6.25 12.5 25 50 100];
    RhoN (1,:) {mustBeLessThan(RhoN,1),mustBeGreaterThan(RhoN,-1)} = -0.5:0.05:0.5;
    RhoD (1,:) {mustBeLessThan(RhoD,1),mustBeGreaterThan(RhoD,-1)} = -0.5:0.05:0.5;
    p.name string = "";  %names will be appended at the end of filename
    p.seed (1,1) {mustBeInteger,mustBePositive} = 96309;
    p.tuning function_handle = @contrastResponse;
end
%% Load files, get loop size
name = p.name;
seed = p.seed;
tuning = p.tuning;
rng(seed) % Reproducibility
baseFolder = fileparts(cd);
fprintf("Current Folder: %s\n", pwd);
fprintf("Base Folder: %s\n",baseFolder);
fprintf([repelem('-',29) 'Loading crf parameters' repelem('-',29) newline])
try
    addpath(genpath(fullfile(baseFolder,'src')));
    load('../data/processed/params_best_JN2019.mat','parampairs');
catch
    disp('Error loading data/adding to path. Check path is correct')
end


if nEXP == 0
    nEXP = size(parampairs,1);
elseif nEXP > 11628
    error("nEXP must be less than 11628");
else
    parampairs = parampairs(randsample(size(parampairs,1),nEXP),:);
end
nRhoNs = length(RhoN);
nRhoDs = length(RhoD);
nContrasts = length(Contrasts);
szIterations = [nEXP nRhoNs nRhoDs];
nIterations = prod(szIterations);
baseName = "SingleTrialNorm";
%% Setup Parallel Pool
fprintf([repelem('-',28) 'Setting up Parallel Pool' repelem('-',28) newline])
poolobj = gcp('nocreate');
if isempty(poolobj)
    if nParallel > 1
        poolobj = parpool('local',nParallel);
        addAttachedFiles(poolobj,regexp(genpath(fullfile(baseFolder,'src')),':', 'split'));
    end
end
paramsFixed = parallel.pool.Constant(parampairs);
szSpace = parallel.pool.Constant(szIterations);
RhoNloop = parallel.pool.Constant(RhoN);
RhoDloop = parallel.pool.Constant(RhoD);
%% Ouput Structure
genParams = struct('nEXP',nEXP,'nTrials',nTrials,'Contrasts',Contrasts,'RhoN',RhoN,'RhoD',RhoD,'time',[]);
CorrNormSignalPairwise = NaN(nIterations,nContrasts);
CorrNormSignalIndependent = NaN(nIterations,nContrasts);
MSENormSignalPairwise = NaN(nIterations,nContrasts);
MSENormSignalIndependent = NaN(nIterations,nContrasts);


objSave = {'genParams',...
    'CorrNormSignalPairwise',...
    'CorrNormSignalIndependent',...
    'MSENormSignalPairwise',...
    'MSENormSignalIndependent'};
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
fprintf("Beginning Run\nTuning function: %s\n# of Iterations: %d\n# of Trials: %d\n",functions(tuning).function,nIterations,nTrials);

fitSTART = datetime('now','Format','eeee, MMMM d, yyyy h:mm a');
fprintf('Start time: %s\n',string(fitSTART));
parfor iIteration=1:nIterations
    [iEXP,iRhoN,iRhoD] = ind2sub(szSpace.Value,iIteration);
    pairLoop = paramsFixed.Value(iEXP,:);

    rN = RhoNloop.Value(iRhoN);
    rD = RhoDloop.Value(iRhoD);

    %% Generate Data from the Ratio of Gaussians
    [spikes,Dtrue,parameters] = generateRogCrf(nTrials,...
        nContrasts,...
        pairLoop{1:12},...
        'cont',Contrasts,...
        'rho_n',rN,...
        'rho_d',rD, ...
        'tuning',tuning);

    %% Compute empirical Noise Correlations
    [~,~,allCorr(iIteration,:)] = mncovEmpirical(spikes);

    [Dpair,Dind] = normalizationSingleTrialInference(spikes,parameters);
    MSENormSignalPairwise(iIteration,:) = mean((Dtrue-Dpair).^2,[1 3],'default','omitnan');
    MSENormSignalIndependent(iIteration,:) = mean((Dtrue-Dind).^2,[1 3],'default','omitnan');
    [CorrNormSignalPairwise(iIteration,:),CorrNormSignalIndependent(iIteration,:)] = corrTrueMap(Dtrue,Dpair,Dind,nContrasts);
end

fitEND = datetime('now','Format','eeee, MMMM d, yyyy h:mm a');
fprintf('Finshed!\nEnd time: %s\n',string(fitEND));
fitDURATION = fitEND - fitSTART;
fprintf('Elapsed Time = %s\n',fitDURATION);
genParams.time = fitDURATION;

allCorr = reshape(allCorr,[szIterations,nContrasts]); %#ok<*NASGU>
MSENormSignalPairwise = reshape(MSENormSignalPairwise,[szIterations,nContrasts]);
MSENormSignalIndependent = reshape(MSENormSignalIndependent,[szIterations,nContrasts]);
CorrNormSignalIndependent = reshape(CorrNormSignalIndependent,[szIterations,nContrasts]);
CorrNormSignalPairwise = reshape(CorrNormSignalPairwise,[szIterations,nContrasts]);

fprintf("\nSaving as %s in ../data/simulations...\n",name)
try
    save(fullfile(baseFolder,"data","simulations",name),objSave{:},'-v7.3');
    fprintf("Success!\n");
catch
    fprintf("Failed to Save!\nSaving everything as FAIL_SINGLETRIALNORM.mat\n")
    save("FAIL_SINGLETRIALNORM.mat",'-v7.3');
end
fprintf([repelem('-',80) newline]);
end

function [corrPair,corrInd] = corrTrueMap(Dt,DmapPair,DmapInd,nStims)
corrPair = NaN(2,nStims);
corrInd = NaN(2,nStims);
for iStim = 1:nStims
    corrPair(:,iStim) = diag(corr(permute(Dt(:,iStim,:),[3 1 2]),permute(DmapPair(:,iStim,:),[3 1 2]),'rows','complete')).';
    corrInd(:,iStim) = diag(corr(permute(Dt(:,iStim,:),[3 1 2]),permute(DmapInd(:,iStim,:),[3 1 2]),'rows','complete')).';
end

corrPair = mean(corrPair,1,'omitnan');
corrInd = mean(corrInd,1,'omitnan');
end
