% load('../data/simulations/SingleTrialNorm.mat')
%% 
% MSENormSignalPairwise - MSE between true denominator and pairwise
%           estimator across trials, varying rho_n,rho_d systematically
% MSENormSignalIndependent - MSE between true denominator and independent
%           estimator across trials, varying rho_n,rho_d systematically
% CorrNormSignalPairwise - Correlation between true denominator and
%   pairwise estimator across trials, varying rho_n,rho_d 
% CorrNormSignalPairwise - Correlation between true denominator and
%   independent estimator across trials, varying rho_n,rho_d 
%% Compute relative %difference
mseDIFF = 100*(-(MSENormSignalPairwise-MSENormSignalIndependent)./MSENormSignalIndependent);
corrDIFF = 100*(CorrNormSignalPairwise-CorrNormSignalIndependent)./CorrNormSignalIndependent;
medianMSE = NaN([size(mseDIFF,2:3) 3]);
medianCorr = NaN([size(corrDIFF,2:3) 3]);
%% Plot Setup
cmap = brewermap(100,'Blues'); %'YlOrBr'
Fig = figure('Color','w','Units','inches','Position',[0.5625    0.5573    9.0938    5.9114]);
t = tiledlayout(2,3);
FigAxes = gobjects(2,3);
colormap(cmap);
nRho = size(mseDIFF,2);
tickLoc = 1:5:nRho;
tickLabels = linspace(-0.5,0.5,numel(tickLoc));
% aData = ones(nRho); aData(floor(nRho/2),floor(nRho/2)) = 0; 
    %%
for iFig = 1:3
    nexttile
    ax=gca();
    imageFig = squeeze(median(mseDIFF(:,:,:,3),1,'omitnan'));
    I = imagesc(ax,imageFig);
    set(ax,'YDir','normal','XTick',tickLoc,'YTick',tickLoc,'XTickLabel',tickLabels,'YTickLabel',tickLabels,'XTickLabelRotation',45,'TickDir','out')
    axis square
    box off
    FigAxes(1,iFig) = ax;
    medianMSE(:,:,iFig) = imageFig;
    title("MSE")
end
%%
for iFig=1:3
    clim(FigAxes(1,iFig),[min(medianMSE(:)),max(medianMSE(:))])
end
c1 = colorbar(FigAxes(1,3));
c1.YLabel.String = 'Improvement over Independent (%)';
%%
for iFig = 1:3
    nexttile
    ax=gca();
    imageFig = squeeze(median(corrDIFF(:,:,:,3),1,'omitnan'));
    I = imagesc(ax,imageFig);
    set(ax,'YDir','normal','XTick',tickLoc,'YTick',tickLoc,'XTickLabel',tickLabels,'YTickLabel',tickLabels,'XTickLabelRotation',45,'TickDir','out')
    axis square
    box off
    FigAxes(2,iFig) = ax;
    medianCorr(:,:,iFig) = imageFig;
end
%%
for iFig=1:3
    clim(FigAxes(2,iFig),[min(medianCorr(:)),max(medianCorr(:))])
end
c2 = colorbar(FigAxes(2,3));
c2.YLabel.String = 'Improvement over Independent (%)';
%%
xlabel(t,{'','\rho_D'},'FontSize',16);
ylabel(t,{'\rho_N',''},'FontSize',16);
ylabel(FigAxes(1,1),{'MSE',''})
ylabel(FigAxes(2,1),{'Corr',''})
title(FigAxes(1,1),'Low Contrast')
title(FigAxes(1,2),'Medium Contrast')
title(FigAxes(1,3),'High Contrast')