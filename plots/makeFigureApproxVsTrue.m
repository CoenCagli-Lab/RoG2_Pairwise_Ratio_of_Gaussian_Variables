load('../data/simulations/ApproxVTrue.mat');
%%
figureApprox = figure('Color','w','Units','inches','Position',[1.9010    0.6094    7.7708    5.6771],'PaperPositionMode','auto');

labelFontSize = 12;
scatterMarkerSize = 50;
t = tiledlayout(2,2); %t.TileSpacing = 'tight'; t.Padding = 'compact';
axApprox = gobjects(2,2);
nExperiments = genParams.nEXP;
IND = 1:10:nExperiments; %For improved speed of scatter plotting, subsample
% IND = 1:nExperiments;
covValues = [ratioCovariance taylorCovariance];
corrValues = [ratioCorrelation taylorCorrelation];
%% - Covariance Scatter
nexttile
ax = gca;
scatter(ax,covValues(IND,1),covValues(IND,2),scatterMarkerSize,'MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.1)
axis(ax,'square')
box off
xlabel(ax,'True','FontSize',labelFontSize,'FontWeight','bold');
ylabel(ax,'Taylor approx.','FontSize',labelFontSize,'FontWeight','bold');
set(ax,'XLim',get(ax,'YLim'));
set(ax,'XTick',get(ax,'YTick'),'XTickLabelRotation',45)
axApprox(1,1) = ax;
title(ax,'Covariance','FontSize',labelFontSize,'FontWeight','bold')
%% - Correlation Scatter
nexttile
ax = gca;
scatter(ax,corrValues(IND,1),corrValues(IND,2),scatterMarkerSize,'MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.1)
axis(ax,'square')
box off
xlabel(ax,'True','FontSize',labelFontSize,'FontWeight','bold');
ylabel(ax,'Taylor approx.','FontSize',labelFontSize,'FontWeight','bold');
set(ax,'XLim',get(ax,'YLim'));
set(ax,'XTick',get(ax,'YTick'),'XTickLabelRotation',45)
axApprox(1,2) = ax;
title(ax,'Correlation','FontSize',labelFontSize,'FontWeight','bold')
%% - Covariance Histogram 
nexttile
ax = gca;
hold on
covHist = histogram(ax,(covValues(:,2)-covValues(:,1))./covValues(:,1).*100,-12.5:1:12.5,'FaceColor',[0.75 0.75 0.75],'EdgeColor','none','Normalization','probability');
axis square
box off
xline(0,'k--','LineWidth',1.5);
scatter(ax,median(covHist.Data,'omitnan'),0.9,scatterMarkerSize,'r','filled','Marker','v')
xlabel(ax,'Approx - True (%)','FontSize',labelFontSize,'FontWeight','bold')
ylabel(ax,'Prop. cases','FontSize',labelFontSize,'FontWeight','bold')
set(ax,'YScale','log')
axApprox(2,1) = ax;
%% - Correlation Histogram
nexttile
ax = gca;
hold on
corrHist = histogram(ax,(corrValues(:,2)-corrValues(:,1))./corrValues(:,1).*100,-12.5:1:12.5,'FaceColor',[0.75 0.75 0.75],'EdgeColor','none','Normalization','probability');
axis square
box off
xline(0,'k--','LineWidth',1.5);
scatter(ax,median(corrHist.Data,'omitnan'),0.9,scatterMarkerSize,'r','filled','Marker','v')
xlabel(ax,'Approx - True (%)','FontSize',labelFontSize,'FontWeight','bold')
ylabel(ax,'Prop. cases','FontSize',labelFontSize,'FontWeight','bold')
set(ax,'YScale','log')
axApprox(2,2) = ax;
%%
for i=1:4; axApprox(i).TickDir='out'; axApprox(i).LineWidth=1.5; end