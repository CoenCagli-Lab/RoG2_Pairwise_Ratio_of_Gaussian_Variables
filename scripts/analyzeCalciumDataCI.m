function [] = analyzeCalciumDataCI(nParallel,nBoots,p)
% analyzeCalciumDataCI Takes Processed Calcium data and fits pairwise Ratio of Gaussians model
%%  INPUT
%       nParallel ( must be > 1) - number of cores to run on parallel
%       nBoots - number of bootstrap resamples (default 1000)
%%% NAME-VALUE INPUTS
%       name (optional) - appended to the end of the file. Full output file name is YYYYMMDD_AnalysisCalciumCI.mat. If you add a name, the file name becomes YYYYMMDD_AnalysisCalciumCI_name.mat
%       seed (optional) - random seed for reproducibility
%       tuning (optional, unused) - tuning function
%%  OUTPUT File Contents
% 		nllPairwise (number of pairs, 1) - median negative log-likelihood for pairwise Ratio of Gaussians
% 		nllIndependent (number of pairs, 1) - median negative log-likelihood for independent Ratio of Gaussians
% 		gofPairwise (number of pairs, 1) - median goodness of fit (normalized negative log-likelihood) for pairwise Ratio of Gaussians
% 		gofIndependent (number of pairs, 1) - median goodness of fit (normalized negative log-likelihood) for independent Ratio of Gaussians
%       nllModGauss (number of pairs, 1) - median negative log-likelihood for the modulated Gaussian model
%       gofModGauss (numbr of pairs, 1) - median goodness of fit for the modulated Gaussian model
% 		parameterstRoG (number of pairs, 19) - best fit parameters to pairwise Ratio of Gaussians
%       parametersMG (number of pairs, 11) - best fit parameters to modulated Gaussian
%       parametersRoGBoot (nummber of pairs, nBoots, 19) - bootstrap fits for the pairwise Ratio of Gaussians
%       parametersMGBoot (number of pairs, nBoots, 11) - bootstrap fit for the modulated Gaussian
%       fitStruct struct (number of pairs) - structure containing fit moments for each pair (according to the index)
%           Contains:
%                   muRoG (2, number of stimuli) - fit mean for Ratio of Gaussians
%                   sigmaRoG (number of stimuli, 2,2) - noise covariance for Ratio of Gaussians
%                   corrRoG (1, number of stimuli) - noise correlations for Ratio of Gaussians
%                   muMG (2, number of stimuli) - fit mean for modulated Gaussian
%                   sigmaMG (number of stimuli, 2,2) - noise covariance for modulated Gaussian
%                   corrMG (1, number of stimuli) - noise correlations for modulated Gaussian
%                   muEmp (2, number of stimuli) - empirical mean
%                   sigmaEmp (number of stimuli, 2,2) empirical noise covariance
%                   corrEmp (1, number of stimuli) - empirical noise correlations
%       runSettings - convenience structure
% File saved in data/processed
%%
arguments
    nParallel (1,1) {mustBeInteger,mustBePositive};
    nBoots (1,1) {mustBeInteger,mustBePositive} = 1000;
    p.name string = "";  %names will be appended at the end of filename
    p.seed (1,1) {mustBeInteger,mustBePositive} = 96309;
    p.tuning function_handle = @contrastResponse;
end
baseFolder = fileparts(cd);
fprintf("Current Folder: %s\nBase Folder: %s\n", pwd,baseFolder);
name = p.name;
seed = p.seed;
tuning = p.tuning;
rng(seed); % Reproducibility
try
    fprintf([repelem('-',31) 'Loading Ca2+ data' repelem('-',31) newline])
    addpath(genpath(fullfile(baseFolder,"src")));
    load(...
        fullfile(baseFolder,"data","processed","PairwiseCalcium_NeuropilCoeff07.mat"), ...
        'Data');
catch
    error('Error loading data/adding to path. Check path is correct')
end
%% Setup Parallel Pool
fprintf([repelem('-',28) 'Setting up Parallel Pool' repelem('-',28) newline])
optsfmin = optimoptions(@fmincon,...
    'Algorithm','sqp',...
    'MaxIter',1000,...
    'Display','off',...
    'TolX',1e-6,...
    'TolFun',1e-6,...
    'TolCon',1e-3);

disp('Setting up parallel pool')
poolobj = gcp('nocreate'); %#ok<*PFOUS>
if isempty(poolobj)
    poolobj = parpool('local',nParallel); %#ok<NASGU>
end
Session = parallel.pool.Constant(Data);

%% Setup File Names
baseName = 'AnalysisCalciumCI';
Today = string(datetime('now','Format','yMMdd'));
saveprefix = strjoin([Today,baseName],"_");
if name ~= ""
    name = strjoin([saveprefix,name],"_");
else
    name = saveprefix;
end
name = name + ".mat";
fprintf('Saving in %s\n',name);

nPairs = length(Data);
% nPairs = 10;
nllPairwise = NaN(nPairs,1);
nllIndependent = NaN(nPairs,1);
gofPairwise = NaN(nPairs,1);
gofIndependent = NaN(nPairs,1);
nllModGauss = NaN(nPairs,1);
gofModGauss = NaN(nPairs,1);
parametersRoG = NaN(nPairs,19);
parametersMG = NaN(nPairs,11);
parametersRoGBoot = NaN(nPairs,nBoots,19);
parametersMGBoot = NaN(nPairs,nBoots,11);
fitStruct = struct(...
    "muRoG",[],...
    "sigmaRoG",[],...
    "corrRoG",[],...
    "muMG",[],...
    "sigmaMG",[],...
    "corrMG",[],...
    "muEmp",[],...
    "sigmaEmp",[],...
    "corrEmp",[]...
    );
fitStruct(nPairs).nllPairwise = [];
runSettings = struct("time",[],"nBoots",nBoots,"tuning",functions(tuning).function);
objSave = {...
    "nllPairwise",...
    'nllIndependent',...
    'nllModGauss',...
    'gofPairwise',...
    'gofIndependent', ...
    'gofModGauss', ...
    "parametersRoG",...
    "parametersMG",...
    "parametersRoGBoot",...
    "parametersMGBoot",...
    'fitStruct',...
    'runSettings'};
fprintf('Saving Variables:');fprintf(' %s ',objSave{:}); fprintf('\n')
fprintf([repelem('-',80) newline]);
fprintf("Beginning Fitting\nRunning %d (# of pairs) Iterations with %d bootstraps:\n",nPairs,nBoots);
fitSTART = datetime('now','Format','eeee, MMMM d, yyyy h:mm a');
fprintf('Start time: %s\n',string(fitSTART));
allt = tic;
parfor iPair=1:nPairs
    %     for iPair=1:nPairs
    loopSession = Session.Value(iPair);
    %         loopSession = Data(iPair);
    loopResponse = loopSession.resp_pair;
    loopMuEta = loopSession.spont_pair;
    loopVarEta = (loopSession.stdspont_pair).^2;
    loopCovEta = diag(loopVarEta); % % No spontaneous correlations
    Contrasts = loopSession.Contrasts;
    [nllPairwise(iPair),...
        nllIndependent(iPair),...
        nllModGauss(iPair),...
        gofPairwise(iPair),...
        gofIndependent(iPair),...
        gofModGauss(iPair),...
        loopRoG,...
        loopModG] =...
        fitRoGCalciumWrapperCV( ...
                                loopResponse,...
                                Contrasts, ...
                                loopMuEta,...
                                loopCovEta, ...
                                optsfmin,...
                                tuning);
    parametersRoG(iPair,:) = loopRoG;
    parametersMG(iPair,:) = loopModG;


    [fitStruct(iPair).muRoG,...
        fitStruct(iPair).sigmaRoG,...
        fitStruct(iPair).corrRoG,...
        fitStruct(iPair).muMG,...
        fitStruct(iPair).sigmaMG,...
        fitStruct(iPair).corrMG] = modelMoments(loopRoG,loopModG,Contrasts,loopMuEta,tuning);
    [loopRoGBoot, loopModGBoot] =...
                    fitCalciumBootstrap(...
                    loopResponse,...
                    Contrasts,...
                    loopMuEta,...
                    loopCovEta,...
                    nBoots,...
                    optsfmin,...
                    tuning);
    parametersRoGBoot(iPair,:,:) = loopRoGBoot;
    parametersMGBoot(iPair,:,:) = loopModGBoot;
    fitStruct(iPair).muEmp = loopSession.mn_resp;
    fitStruct(iPair).sigmaEmp = loopSession.cov_resp;
    fitStruct(iPair).corrEmp = loopSession.noiseCorr;
    fitStruct(iPair).pairResponse = loopSession.resp_pair;
end

fitEND = datetime('now','Format','eeee, MMMM d, yyyy h:mm a');
fprintf('Finshed!\nEnd time: %s\n',string(fitEND));
fitDURATION = fitEND - fitSTART;
t=toc(allt);
runSettings.time = t; %#ok<STRNU>
fprintf('Elapsed Time = %s\n',fitDURATION);
fprintf("\nSaving as %s in ../data/analysis...\n",name)

save(fullfile(baseFolder,"data","analysis",name),objSave{:},'-v7.3')
end
%%
%% ---------------------------------------------------------------------- %
%% AUXILIARY FUNCTIONS
%% ---------------------------------------------------------------------- %
%%
function [mu_rog, sigma_rog, corr_rog, mu_mg, sigma_mg, corr_mg] = modelMoments(paramRoG,paramModG,cont,mu_eta,tuning)
[mu_rog,sigma_rog,corr_rog] = mncovRogCrf(paramRoG,cont,mu_eta,tuning);
[mu_mg, sigma_mg, corr_mg] = mncovMgCrf(paramModG,cont,mu_eta,tuning);
end

function [nll_pair,nll_ind,nll_mg,gof_pair,gof_ind,gof_mg,pfit_rog,pfit_mg] = fitRoGCalciumWrapperCV(response,contrasts,mu_eta,cov_eta,opts,tuning)

var_eta = diag(cov_eta);

rho_eta = cov_eta(1,2)./(sqrt(prod(var_eta)));

nTrials = size(response,3);


nFold = nTrials;
cvid = 1 + mod((1:nTrials)',nFold);
cvid = cvid(randperm(nTrials));
nll_pair = NaN(nFold,1);
nll_ind = NaN(nFold,1);
nll_pair_norm = NaN(nFold,1);
nll_ind_norm = NaN(nFold,1);
nll_mg = NaN(nFold,1);
nll_mg_norm = NaN(nFold,1);
for iFold=1:nFold

    %% CV Setting
    cv_set = (cvid~=iFold);
    r_train = response(:,:,cv_set);
    r_test = response(:,:,~cv_set);
    [MU_HAT_TRAIN,SIGMA_HAT_TRAIN,~] = mncovEmpirical(r_train);
    [x_pair,x_ind] = fitRogCrfTwoStepInner(MU_HAT_TRAIN,SIGMA_HAT_TRAIN,contrasts,mu_eta,var_eta,rho_eta,tuning,opts);
    [x_mg,~] = fitMgCrfTwoStepInner(MU_HAT_TRAIN,SIGMA_HAT_TRAIN,contrasts,mu_eta,var_eta,rho_eta,tuning,opts);
    [MU_PAIR,SIGMA_PAIR] = mncovRogCrf(x_pair,contrasts,mu_eta,tuning);
    [MU_IND,SIGMA_IND] = mncovRogCrf(x_ind,contrasts,mu_eta,tuning);
    [MU_NB,SIGMA_NB] = mncovMgCrf(x_mg,contrasts,mu_eta,tuning);



    %% Computing Goodness of Fit Metrics
    % %  Null Model % %
    tmp = -log_bvnpdf(reshape(r_test,2,[]).',...
        mean(r_train,[2 3],'omitnan').',...
        diag(var(reshape(r_train,2,[]),[],2,'omitnan')));
    tmp(isinf(tmp)) = NaN;
    nllnull = mean(tmp,'omitnan');
    %% ---------------------------------------------------------------------- %
    % % 5/24/22 - These are set this way to reproduce the independent RoG
    % % oracle/null structure and to allow for a fair comparison between
    % % independent and pairwise models
    % ----------------------------------------------------------------------- %
    SIGMA_HAT_TRAIN(:,1,2) = 0;
    SIGMA_HAT_TRAIN(:,2,1) = 0;
    %%
    %%
    nlloracle = rogNegLogLikeFull(r_test,MU_HAT_TRAIN,SIGMA_HAT_TRAIN);
    nllpair = rogNegLogLikeFull(r_test,MU_PAIR,SIGMA_PAIR);
    nllind = rogNegLogLikeFull(r_test,MU_IND,SIGMA_IND);

    nll_pair_norm(iFold) = (nllpair-nllnull)./(nlloracle-nllnull);
    nll_pair(iFold) = nllpair;
    nll_ind(iFold) = nllind;
    nll_ind_norm(iFold) = (nllind-nllnull)./(nlloracle-nllnull);

    nllmg = rogNegLogLikeFull(r_test,MU_NB,SIGMA_NB);
    nll_mg_norm(iFold) = (nllmg-nllnull)./(nlloracle-nllnull);
    nll_mg(iFold) = nllmg;
end
pfit_rog = fitRogCrfTwoStep(response,contrasts,mu_eta,cov_eta,"optimopts",opts,"tuning",tuning);
pfit_mg = fitMgCrfTwoStep(response,contrasts,mu_eta,cov_eta,"optimopts",opts,"tuning",tuning);
nll_pair = median(nll_pair,1,'omitnan');
nll_ind = median(nll_ind,1,'omitnan');
gof_pair = median(nll_pair_norm,1,'omitnan');
gof_ind = median(nll_ind_norm,1,'omitnan');
nll_mg = median(nll_mg(:),1,'omitnan');
gof_mg = median(nll_mg_norm(:),'omitnan');

end

function [fitRoG,fitModG] = fitCalciumBootstrap(response,contrasts,mu_eta,cov_eta,nBoots,opts,tuning)
nT = size(response,3);
fitRoG = NaN(nBoots,19);
fitModG = NaN(nBoots,11);
for iBoot=1:nBoots
    r = response(:,:,randsample(nT,nT,true));
    fitRoG(iBoot,:) = fitRogCrfTwoStep(r,contrasts,mu_eta,cov_eta,"optimopts",opts,"tuning",tuning);
    fitModG(iBoot,:) = fitMgCrfTwoStep(r,contrasts,mu_eta,cov_eta,"optimopts",opts,"tuning",tuning);
end
end
