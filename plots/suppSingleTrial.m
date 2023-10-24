nSimulations = 10000;
nContrasts = 10;
nTrials = 100;
contLower = 20;
contUpper = 50;
Contrasts = linspace(contLower,contUpper,nContrasts);



Dtrue = NaN(nSimulations,2,nContrasts,nTrials);
Dmap = NaN(size(Dtrue));
VarN = NaN(nSimulations,2,nContrasts);
VarD = NaN(size(VarN));
corrD = NaN(nSimulations,2);
for iSim = 1:nSimulations
    Rmax = 10+90*rand(2,1);
    eps = 10+15*rand(2,1);
    beta_n = 1 + 0.5*rand(2,1);
    beta_d = 1 + 0.5*rand(2,1);


    [spikes,D,parameters] = generateRogCrf(nTrials,nContrasts, ...
        "Rmax",Rmax, ...
        "eps",eps, ...
        "beta_n", beta_n, ...
        "beta_d",beta_d, ...
        "mu_eta",zeros(2,1), ...
        "var_eta",zeros(2,1));
    Dtrue(iSim,:,:,:) = D;
    
    [Dfit,~] = normalizationSingleTrialInference(spikes,parameters);

    corrD(iSim,:) = diag(corr( ...
        permute(reshape(zscore(D,0,3),2,[]),[2 1]), ...
        permute(reshape(zscore(Dfit,0,3),2,[]),[2 1]) ...
        ));
    Dmap(iSim,:,:,:) = Dfit;
    VarN(iSim,:,:) = parameters.var_n;
    VarD(iSim,:,:) = parameters.var_d;


end