function [] = genDataNoiseCorrVsNormalization(nParallel, nEXP,nTrials,RhoSpace,Contrasts,mu_eta,sigma_eta,name)
%% genDataNoiseCorrVsNormalization simulates data from Ratio of Gaussians to compare how normalization
% effects noise correlations.
% 	Keep RhoSpace small for visualization purposes
%	OUTPUT File Contents
%		genParams (struct) - returns the generating parameters for the function
% 		epsNorm (nEXP,2,number of Rhos, number of Rhos) - product of epsilon parameters for
% 															the paiwise ratio of gaussians.
% 															Contrast independent measure of
% 															normalization strength
% 		allCorr (nEXP,number of contrasts, number of Rhos, number of Rhos) - all noise
% 															correlations
% Files are saved in data/simulations
%%
	arguments
		nParallel;
		nEXP = 1e5;
		nTrials = 1e3;
		RhoSpace = [-0.5 0 0.5];
		Contrasts = [6.25 12.5 25 50 100];
		mu_eta = [0;0];
		sigma_eta = [0 0;0 0];
		name = '';
	end
	rng(96309) % Reproducibility
	disp(['Current Folder: ',pwd])
	try
		addpath(genpath('../src/'));
	catch
		disp('Error loading data/adding to path. Check path is correct')
	end

	%% Setup Parallel Pool
	poolobj = gcp('nocreate'); %#ok<*PFOUS> 
	if isempty(poolobj)
		poolobj = parpool('local',nParallel);
	end

	%% Setup File Names
	Today = string(datetime('now','Format','yMMdd'));
	saveprefix = sprintf('%s_DataFigureNoiseCorrVsNormalization',Today);
	if ~isempty(name)
		name = char(sprintf('%s_%s.mat',saveprefix,name));
	else
		name = char([saveprefix,'.mat']);
	end

	fprintf('File will be saved as %s in data/simulations folder\n',name);
	fprintf('Begin Simulation\n');
	tic;
	[epsNorm,allCorr] = noiseCorrerlationNormalization_generateData(...
						nEXP,...
						nTrials,...
						RhoSpace,...
						Contrasts,...
						mu_eta,...
						diag(sigma_eta)); %#ok<*ASGLU> 

	t = toc;
	fprintf('Completed in %0.2f minutes\n',t/60)

	genParams = struct();
	genParams.NEXP = nEXP;
	genParams.nTrials = nTrials;
	genParams.RhoSpace = RhoSpace;
	genParams.Contrasts = Contrasts;
	genParams.mu_eta = mu_eta;
	genParams.sigma_eta = sigma_eta;
	objSave = {'epsNorm','allCorr','genParams'};
	fprintf('Saving as %s.mat in ../data/simulations\n',name)
	save(['../data/simulations/',name],objSave{:},'-v7.3');

end

function [GeoMeanNormalization,NoiseCorrelations] = noiseCorrerlationNormalization_generateData(nExperiments,nTrials,Rhos,Contrasts,spontmn,spontvar)
    nRho = numel(Rhos);
    nContrasts = numel(Contrasts);

    NoiseCorrelations=NaN(nExperiments,nContrasts,nRho,nRho);

    % GeoMeanNormalization = NaN(nExperiments,nContrasts,nRho,nRho);
    GeoMeanNormalization = NaN(nExperiments,2,nRho,nRho); %% - 5/30/22 - Epsilon instead of geomean(mu_D)
    parfor iExperiment = 1:nExperiments
    C = NaN(nContrasts,nRho,nRho);
        % Dtmp = NaN(nContrasts,nRho,nRho);
        Dtmp = NaN(2,nRho,nRho);
        for iRhoN = 1:nRho
            for iRhoD = 1:nRho
                [spikes,~,p] = generateRogCrf(nTrials,...
				nContrasts,...
				'rho_n',Rhos(iRhoN),...
				'rho_d',Rhos(iRhoD),...
				'cont',Contrasts,...
				'mu_eta',spontmn,...
				'var_eta',spontvar);
                % for iContrast=1:nContrasts
                %     tmp = corr(squeeze(spikes(:,iContrast,:))','rows','complete');
                %     C(iContrast,iRhoN,iRhoD) = tmp(1,2);
                % end
				C(:,iRhoN,iRhoD) = noisecorrEmpirical(spikes,nContrasts);
                % Dtmp(:,iRhoN,iRhoD) = prod(p.mu_d,1).^(0.5);
                Dtmp(:,iRhoN,iRhoD) = p.eps; %% - 5/30/22 - Epsilon instead of geomean(mu_D)
            end
        end
        NoiseCorrelations(iExperiment,:,:,:) = C;
        % GeoMeanNormalization(iExperiment,:,:,:) = Dtmp;
        GeoMeanNormalization(iExperiment,:,:,:) = Dtmp; %% - 5/30/22 - Epsilon instead of geomean(mu_D)

    end
end
