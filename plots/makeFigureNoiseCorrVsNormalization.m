load('../data/simulations/NoiseCorrVsNormalization.mat')
%% - Rho Values are -0.5,0,0.5, so there are 9 possible combinations of rho_N,rho_D.
% % rho vector is -0.5:0.5:0.5, and the pairs below are [rho_N,rho_D]
% % locations, so [2,3] => rho_N = 0, rho_D = 0.5 and similarly
rhoValues = struct();
rhoValues.Positive = {[2 3],[3 2], [3 3]};
rhoValues.Negative = {[2 1],[1 2], [1 1]};
rhoValues.Mixed = {[1 3],[3 1]};
prodEps = squeeze(prod(epsNorm,2));
%% - Main Figure - Positive Rho Values
[figPositive,axPositive] = plotFigureNoiseCorrNorm(allCorr,prodEps,genParams,'ContLevels',[1 2 3 5],'FigRho',rhoValues.Positive);
%% - Supplement Figures - Negative Rho Values
[figNegative,axNegative] = plotFigureNoiseCorrNorm(allCorr,prodEps,genParams,'ContLevels',[1 2 3 5],'FigRho',rhoValues.Negative);
%% - Supplement Figures - Mixed Positive and Negative
[figMixed,axMixed] = plotFigureNoiseCorrNorm(allCorr,prodEps,genParams,'ContLevels',[1 2 3 5],'FigRho',rhoValues.Mixed,'Position',[0.5625    2.0833    8.3906    4.3854
]);
yline(axMixed(1),0,'--','LineWidth',1.5)
yline(axMixed(2),0,'--','LineWidth',1.5)
Leg = figMixed.Children(1).Children(1).String;
figMixed.Children(1).Children(1).String = Leg(1:end-1);
%%
function [figureNC,FigAxes] = plotFigureNoiseCorrNorm(allCorr,eps,genParams,p)
    arguments
        allCorr;
        eps;
        genParams;
        p.ContLevels = [1 2 3 5];
        p.FigRho = {[2 3],[3 2],[3 3]};
        p.MarkerSize = 26;
        p.EdgeWidth = 0.1;
        p.MarkerEdgeColor = 'none';
        p.Legend = 1;
        p.Position (1,4) = [0.5625    2.0833   11.6406    4.3854];
        p.Units char = 'inches';
    end
    ContLevels =p.ContLevels;
    FigRho = p.FigRho;
    nContrasts = numel(genParams.Contrasts);
    rhos = genParams.RhoSpace;
    nRhos = numel(rhos);
%     NoiseCorrMedian = cell(nContrasts,nRhos,nRhos);
%     CenterBins = cell(nContrasts,nRhos,nRhos);
%     NoiseCorrCI = cell(nContrasts,nRhos,nRhos);
    szLoop = [nContrasts,nRhos,nRhos];
    nIterations = prod(szLoop);
    NoiseCorrMedian = cell(nIterations,1);
    CenterBins = cell(nIterations,1);
    NoiseCorrCI = cell(nIterations,1);

    % tic;
%     for iContrast = 1:nContrasts
%         for iRhoN = 1:nRhos
%             for iRhoD = 1:nRhos
    for iIteration=1:nIterations
        [iContrast,iRhoN,iRhoD] = ind2sub(szLoop,iIteration);
                % % Bin by product of epsilon, exclude unused bins, compute
                % % median in this bin, exclude bins with <100 values
                [N,EDGES,BINIDS] = histcounts(squeeze(eps(:,iRhoN,iRhoD)));
                BINEX = or(isnan(BINIDS),BINIDS == 0);
                BINIDS(BINEX) = [];
                NoiseCorrBin = splitapply(@(x){x},allCorr(~BINEX,iContrast,iRhoN,iRhoD),BINIDS);
                NoiseCorrBin = NoiseCorrBin(N>100);
                NoiseCorrBinMedian = cell2mat(cellfun(@(x)(mean(x,1,'omitnan')),NoiseCorrBin,'UniformOutput',0));
                % NoiseCorrBinCI = cell2mat(cellfun(@(x)(bootci(10000,{@(y)(mean(y,1,'omitnan')),x})),NoiseCorrBin,'UniformOutput',0).');
%                 NoiseCorrBin = accumarray(BINIDS,squeeze(allCorr(~BINEX,iContrast,iRhoN,iRhoD)),[numel(EDGES)-1 1],@(x)(median(x,'omitnan')));
                CENTERS = 1/2*(EDGES(1:end-1)+EDGES(2:end));
%                 NoiseCorrMedian{iContrast,iRhoN,iRhoD} = NoiseCorrBinMedian;
%                 NoiseCorrCI{iContrast,iRhoN,iRhoD} = NoiseCorrBinCI;
%                 CenterBins{iContrast,iRhoN,iRhoD} = CENTERS(N>100);
                NoiseCorrMedian{iIteration} = NoiseCorrBinMedian;
                % NoiseCorrCI{iIteration} = NoiseCorrBinCI;
                CenterBins{iIteration} = CENTERS(N>100);
    end
    % toc

    nContPlot = numel(ContLevels);
    cmap = brewermap(nContPlot+1,'YlOrRd');

    figureNC = figure('Color','w');
    figureNC.Units = p.Units;
    figureNC.Position = p.Position;
    nPlots = numel(p.FigRho);
    t = tiledlayout(1,nPlots); %t.TileSpacing = 'compact'; t.Padding = 'compact';
    ContNames = string(genParams.Contrasts(ContLevels));
    FigAxes = gobjects(1,nPlots);
    hobj = gobjects(1,nContPlot);
    CenterBins = reshape(CenterBins,szLoop);
    NoiseCorrMedian = reshape(NoiseCorrMedian,szLoop);
    % NoiseCorrCI = reshape(NoiseCorrCI,szLoop);
    for iFig = 1:nPlots
        nexttile
        ax = gca;
        axis(ax,'square');
        hold(ax,'on')
        ND = FigRho{iFig};

       for iContrast = 1:nContPlot
           hobj(iContrast) = scatter(ax,CenterBins{ContLevels(iContrast),ND(1),ND(2)},NoiseCorrMedian{ContLevels(iContrast),ND(1),ND(2)},p.MarkerSize,'filled','MarkerFaceColor',cmap(iContrast+1,:),'MarkerEdgeColor',p.MarkerEdgeColor,'LineWidth',p.EdgeWidth);
            % errorbar(gca,...
            %     CenterBins{ContLevels(iContrast),ND(1),ND(2)},...
            %     NoiseCorrMedian{ContLevels(iContrast),ND(1),ND(2)},...
            %     NoiseCorrCI{ContLevels(iContrast),ND(1),ND(2)}(1,:),...
            %     NoiseCorrCI{ContLevels(iContrast),ND(1),ND(2)}(2,:),...
            %     'LineStyle','none','Marker','o',...
            %     'MarkerEdgeColor',cmap(iContrast+1,:),...
            %     'MarkerFaceColor',cmap(iContrast+1,:),...
            %     'MarkerSize',10,...
            %     'CapSize',0,...
            %     'LineWidth',2,...
            %     'Color',cmap(iContrast+1,:));
       end
       set(ax,'FontName','Arial','LineWidth',1.5)
       yt=get(ax,'YTick');
       TMP = (get(ax,'YLim'));
       yt(1) = TMP(1); yt(end) = TMP(2);
       set(ax,'YTick',[yt(1) 0.5*(yt(1)+yt(end)) yt(end)],'XTick',[0 5000 10000])
       ax.YAxis.TickLabelFormat = '%.2f';
       set(ax,'TickDir','out')

        title(ax,{sprintf('\\rho_N = %0.1f, \\rho_D = %0.1f',rhos(ND(1)),rhos(ND(2))),''},'FontWeight','bold','Color','k','FontName','Arial')
        FigAxes(1,iFig) = ax;
    end
    if p.Legend
        hCopy = copyobj(hobj,ax);
        set(hCopy,'XData', NaN', 'YData', NaN)

        L = legend(hCopy(end:-1:1),ContNames(end:-1:1));
        L.Layout.Tile = 'east'; L.Box = 'off';
    end
    xlabel(t,{'','','Normalization Strength (\epsilon_1\times\epsilon_2)'},'FontWeight','bold','Color','k','FontName','Arial')
    ylabel(t,{'','','Noise Correlations'},'FontWeight','bold','Color','k','FontName','Arial')
end