function [] = genDataApproxVsTrue(nParallel,nEXP,nTrials,name)
% genDataApproxVsTrue compares covariance/correlation from the Ratio of Gaussians distribution
%%  INPUT
%       nParallel - number of cores (must be > 1)
%       nEXP - number of simulated pairs (default 1e4)
%       nTrials - number of simulated trials (default 1e6)
%       name - name appended to end of default name (default name - YYYYMMDD_ApproxVsTrue
%%  OUTPUT File Contents
%		genParams (struct) - returns the generating parameters for the function
% 		ratioCovariance (nEXP,1) - covariance between Ratio of Gaussians using full distribution
% 		ratioCorrelation (nEXP,1) - correlation between Ratio of Gaussians using full distribution
% 		taylorCovariance (nEXP,1) - covariance between Ratio of Gaussians using Taylor Approximation
% 		taylorCorrelation (nEXP,1) - correlation between Ratio of Gaussians using Taylor Approximation
%       genParams - convenience struct containing input parameters
% File saved in data/simulations
%%
	arguments
		nParallel;
		nEXP = 1e4;
		nTrials = 1e6;
		name char = ''; % appended to end of file name
	end
	rng(96309) % Reproducibility
	disp(['Current Folder: ',pwd])
	try
		addpath(genpath('../src/'));
	catch
		disp('Error loading data/adding to path. Check path is correct')
	end

	ratioCovariance = NaN(nEXP,1); %#ok<*PFOUS> 
	taylorCovariance = NaN(nEXP,1);
	ratioCorrelation = NaN(nEXP,1);
	taylorCorrelation = NaN(nEXP,1);

	%% Setup Parallel Pool
	poolobj = gcp('nocreate');
	if isempty(poolobj)
		poolobj = parpool('local',nParallel);
	end

	%% Setup File Names
	Today = string(datetime('now','Format','yMMdd'));
	savename = sprintf('%s_ApproxVsTrue',Today);
	if ~isempty(name)
		savename=[savename,'_',name];
	end

	fprintf('File will be saved as %s\n',savename);
	fprintf('Begin Generation...\n')
	tic;
	parfor iEXP = 1:nEXP
		mu_n = 100*rand(2,1);
		mu_d = 0.5+rand(2,1);
		mu_eta = zeros(2,1);

		alpha_n = ones(2,1);
		alpha_d = 0.001*ones(2,1);

		beta_n = 1+0.5*rand(2,1);
		beta_d = 1+0.5*rand(2,1);

		rho_n = -0.5+rand();
		rho_d = -0.5+rand();
		rho_n2d1 = 0;
		rho_n1d2 = 0;

		rho = [rho_n rho_d rho_n1d2 rho_n2d1];
		var_eta = 0.1*mu_n./mu_d;
		rho_eta = 0;
		c_eta = rho_eta*sqrt(var_eta(1)*var_eta(2));
        cov_eta = diag(var_eta)+diag(c_eta,1)+diag(c_eta,-1);

        loopCovND = covRogFull(alpha_n,beta_n,mu_n,alpha_d,beta_d,mu_d,rho);
    	[loopCovR,taylorCorrelation(iEXP)] = covRogTaylorApprox(alpha_n,beta_n,mu_n,alpha_d,beta_d,mu_d,rho,cov_eta); %#ok<*PFOUS> 
		ND = [mu_n(1);mu_d(1);mu_n(2);mu_d(2)].*ones(1,nTrials)...
				+ chol(squeeze(loopCovND)).'*randn(4,nTrials);
		spikes = [squeeze(ND(1,:)./ND(2,:)); squeeze(ND(3,:)./ND(4,:))] + mvnrnd(mu_eta,cov_eta,nTrials).';
		loopCovFull = cov(spikes.','omitrows');
		ratioCovariance(iEXP) = loopCovFull(1,2);
		ratioCorrelation(iEXP) = loopCovFull(1,2)/sqrt(loopCovFull(1,1)*loopCovFull(2,2));
        taylorCovariance(iEXP) = loopCovR(1,1,2);
	end

	t=toc;
	fprintf('Total = %.2f minutes\n',t/60)
	genParams = struct('nEXP',nEXP,'nTrials',nTrials,'time',seconds(t)); %#ok<*NASGU> 
	objSave  = {'ratioCovariance','ratioCorrelation','taylorCovariance','taylorCorrelation','genParams'};
	save(['../data/simulations/',savename,'.mat'],objSave{:},'-v7.3')
end
