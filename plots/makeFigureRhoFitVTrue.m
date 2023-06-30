
%% load data
% load('../data/simulations/RhoFitVTrue.mat');
%% Seeting up Data 
%  Rho_N,Rho_D
rhoBoot = fitParametersBoot(:,:,13:14); 
rhoFit = fitParametersAll(:,13:14);
rhoTrue = trueParameters(:,13:14);
medCorr = abs(median(allCorr,2,'omitnan'));
[nPairs,nBoots,~] = size(rhoBoot);
rhoBootCI = prctile(rhoBoot,[2.5,97.5],2);
rhoBootMean = squeeze(mean(rhoBoot,2,'omitnan')); % Mean of bootstrap sample

%% Index Filtering
% IndSig - 95% CI excludes 0
INDSIG = squeeze(~and(rhoBootCI(:,1,:)<0,rhoBootCI(:,2,:)>0));
% IndGof - gof for pairwise model is better than independent
INDGOF = and(gofPairwise>gofIndependent,gofPairwise > 0.5);
% IND - the indices that meet both criteria
% INDop - the indices that are not significant but still meet gof criterion
IND = logical(INDSIG.*INDGOF);
INDop = logical(~INDSIG.*INDGOF);

%% RhoFit vs. RhoTrue Scatter
% Colors
cols = NaN(3,3,2);
cols(:,:,1) = brewermap(3,'Blues');
cols(:,:,2) = brewermap(3,'Reds');
cols = cols([1 3],:,:);
%
figureDirect = figure("Name","Fit v True");
tiledlayout(2,1)
for i=1:2
    colFig = squeeze(cols(:,:,i));
    INDSCATER = INDSIG(:,i);
    nexttile
    hold on
    s = scatter(rhoTrue(INDSCATER,i),rhoFit(INDSCATER,i),36,colFig(2,:),'filled','MarkerFaceAlpha',0.5);
    scatter(rhoTrue(~INDSCATER,i),rhoFit(~INDSCATER,i),36,colFig(1,:),'filled','MarkerFaceAlpha',0.5)
    uistack(s,'top');
    hold off
    xlabel('\rho True')
    ylabel('\rho Fit')
    xlim([-1 1])
    ylim([-1 1])
    [r1,p1] = corr(rhoTrue(INDSCATER,i),rhoFit(INDSCATER,i),"rows","complete",'type','Pearson');
    title(sprintf('# Cases Significant: %d/%d\nPearson Correlation = %0.2f, P-value = %0.4e\nProp. cases same sign = %d/%d',...
    sum(INDSCATER),...
    nPairs,...
    r1,...
    p1,...
    sum(sign(rhoTrue(INDSCATER,i))==sign(rhoFit(INDSCATER,i)),"all"),...
    sum(INDSCATER)));
    axis square
end
%% RhoFit v. RhoTrue Bootstrap MSE/Bias/Variance
tileArray = @(oldArray)(repmat(reshape(oldArray,nPairs,[],2),1,nBoots,1));
MSE = squeeze(mean((tileArray(rhoTrue) - rhoBoot).^2,2,'omitnan'));
BIAS = (rhoBootMean - rhoTrue).^2;
VAREST = squeeze(var(rhoBoot,[],2,'omitnan'));
figureMSE = figure("Name","MSE");
names = ["\rho_N","\rho_D"];
tileHandle = tiledlayout(2,3);
title(tileHandle,'Bootstrapped Estimator Stats for \rho_N and \rho_D');
sameSign = NaN(2,1);
for i=1:2
    colFig = squeeze(cols(:,:,i));
%         colSig = (1-INDSIG(:,i)).*colFig(1,:) + INDSIG(:,i) .* colFig(2,:);
%         aData = 0.5 .*INDSIG(:,i) + 0.25 .* (1-INDSIG(:,i));
    loopInd = IND(:,i);
    loopIndOp = INDop(:,i);
    nexttile
%         scatter(medCorr,MSE(:,i),36,colSig,'filled','MarkerFaceAlpha',0.5);
    hold on
    s = scatter(medCorr(loopInd),MSE(loopInd,i),36,colFig(2,:),'filled','MarkerFaceAlpha',0.5);
    scatter(medCorr(loopIndOp),MSE(loopIndOp,i),36,colFig(1,:),'filled','MarkerFaceAlpha',0.5);
    uistack(s,'top');
    hold off
    xlabel('Median Noise Correlation')
    ylabel('MSE')
%         title(names(i),'FontSize',14,'FontWeight','bold');
    title('MSE as Function of Noise Correlation','FontSize',11,'FontWeight','normal');
    ylim([0 2]);
    nexttile
    hold on
%         scatter(medCorr,(BIAS(:,i)),36,colSig,'filled','MarkerFaceAlpha',0.5);
    s = scatter(medCorr(loopInd),BIAS(loopInd,i),36,colFig(2,:),'filled','MarkerFaceAlpha',0.5);
    scatter(medCorr(loopIndOp),BIAS(loopIndOp,i),36,colFig(1,:),'filled','MarkerFaceAlpha',0.5);
    uistack(s,'top');
    hold off
    xlabel('Median Noise Correlation')
    ylabel('Bias of Estimator ^2')
    title('Bias ^2 as Function of Noise Correlation','FontSize',11,'FontWeight','normal');
        ylim([0 2]);

    nexttile
    hold on
%         scatter(medCorr,VAREST(:,i),36,colSig,'filled','MarkerFaceAlpha',0.5);
    s = scatter(medCorr(loopInd),VAREST(loopInd,i),36,colFig(2,:),'filled','MarkerFaceAlpha',0.5);
    scatter(medCorr(loopIndOp),VAREST(loopIndOp,i),36,colFig(1,:),'filled','MarkerFaceAlpha',0.5);
    uistack(s,'top');
    hold off
    xlabel('Median Noise Correlation')
    ylabel('Variance of Estimator')
    title('Variance of Estimator as Function of Noise Correlation','FontSize',11,'FontWeight','normal');
            ylim([0 2]);
    sameSign(i) = sum(sign(rhoTrue(loopInd,i))==sign(rhoFit(loopInd,i)));
end
propCases = struct('SIGGOF',sum(IND,1),'NONSIGGOF',sum(INDop,1),'GOF',sum(INDGOF),'SIG',sum(INDSIG,1),'sameSign',sameSign);

% function newArray = tileArray(oldArray)
%     newArray = repmat(reshape(oldArray,nPairs,[],2),1,nBoots,1);
% end
