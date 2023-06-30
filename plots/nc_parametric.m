function [contData,rhonData,rhodData,rmData,epsData,contID,YLimLine,YLimMesh,colors,cmap] = nc_parametric(rho_n,rho_d,cont,flag)
% This generates Supplementary Figure S4
    arguments
        rho_n (1,1) {mustBeNumeric};
        rho_d (1,1) {mustBeNumeric};
        cont (1,:) {mustBeNumeric} = 1:100;
        flag.type string = "diag";
        flag.common logical = 1;
        flag.save logical =1;
        flag.visible string = "on";
        flag.exportVars logical = 1;
    end
    params ={'cont',cont,'Rmax',[25 50],'eps',[10 25],'rho_n',rho_n,'rho_d',rho_d,'alpha_n',1*ones(2,1),'alpha_d',1*ones(2,1),'beta_n',1.5*ones(2,1),'beta_d',1.5*ones(2,1)};
    [spikes,~,p] = generateRogCrf(1000,params{:});
    nContrasts = size(spikes,2);
    % rhoSpace = -0.9:0.1:0.9;
    rhoSpace = linspace(-1,1,nContrasts);
    eps = linspace(5,50,nContrasts);
    rm = linspace(5,50,nContrasts);
    [~,~,nc_cont] = mncovRogCrf(p.pars,p.cont,p.mu_eta);
    % nc_cont = noisecorrEmpirical(spikes,nContrasts);
    nc_rhon = NaN(numel(rhoSpace),nContrasts);
    nc_rhod = NaN(numel(rhoSpace),nContrasts);
    nc_rm = NaN(nContrasts,nContrasts,nContrasts);
    % rm_strength = NaN(numel(rm),1);
    % rm2 = NaN(numel(rm),2);
    nc_eps = NaN(nContrasts,nContrasts,nContrasts);
    % nc_eps = NaN(numel(eps),nContrasts);
%     eps_strength = NaN(numel(eps),1);
%     eps2 = NaN(numel(eps),2);
%     alpha = linspace(1,10,nContrasts);
%     nc_alphan = NaN(nContrasts,nContrasts);
%     nc_alphad = NaN(nContrasts,nContrasts);
%     beta = linspace(1,2,nContrasts);
%     nc_betan = NaN(nContrasts,nContrasts);
%     nc_betad = NaN(nContrasts,nContrasts);
%     ID = floor(nContrasts/30);
%     IDscale = [1,5,25];
    contID = find(ismember(cont,[3,15,75]));
    ncontIDs = numel(contID);
    cmap = {brewermap(1000,'Blues'),brewermap(1000,'Greens'),brewermap(1000,'Reds')};
    colors = cell2mat(arrayfun(@(x)(cmap{x}(600,:)),1:3,'UniformOutput',false).');
    for i=1:nContrasts
        pLoop = parametersRogCrf(1000,params{:},'rho_n',rhoSpace(i));
        nc_rhon(i,:) = noiseCorrRogCrf(pLoop);
        pLoop = parametersRogCrf(1000,params{:},'rho_d',rhoSpace(i));
        nc_rhod(i,:) = noiseCorrRogCrf(pLoop);
        for j=1:nContrasts
    %     loopRm = rm(i)+5*rand(2,1);
    %     loopRm = min(rm(i)*rand(2,1),rm(i)/1.5);
    %     loopRm = rm(i)*ones(2,1);
    %     pLoop = parametersRogCrf(1000,params{:},'Rmax',loopRm);
    %     [s,~,~] = generateRogCrf(1000,params{:},'Rmax',loopRm);
            pLoop = parametersRogCrf(1000,params{:},'Rmax',[rm(i); rm(j)]);
            nc_rm(i,j,:) =  noiseCorrRogCrf(pLoop);
    %     nc_rm(i,:) = noisecorrEmpirical(s,nContrasts);
    %     rm_strength(i) = sqrt(prod(loopRm(:)));
    %     rm2(i,:) = loopRm;
    %     loopEps = eps(i)+5*rand(2,1);
    %     loopEps = min(eps(i)*rand(2,1),eps(i)/1.5);
    %     loopEps = eps(i)*ones(2,1);
    %     eps2(i,:) = loopEps;
    %     loopEps = eps(i)*ones(2,1);
    %     pLoop = parametersRogCrf(1000,params{:},'eps',loopEps);
            pLoop = parametersRogCrf(1000,params{:},'eps',[eps(i); eps(j)]);
            nc_eps(i,j,:) = noiseCorrRogCrf(pLoop);
        end
    %     [s,~,~] = generateRogCrf(1000,params{:},'Rmax',loopRm);
    %     nc_eps(i,:) = noisecorrEmpirical(s,nContrasts);
    %     eps_strength(i) = sqrt(prod(loopEps(:)));
    
    %     loopAn = alpha(i)*ones(2,1);
    %     pLoop = parametersRogCrf(1000,params{:},'alpha_n',loopAn);
    %     nc_alphan(i,:) = noiseCorrRogCrf(pLoop);
    %     loopAd = alpha(i)*ones(2,1);
    %     pLoop = parametersRogCrf(1000,params{:},'alpha_n',loopAd);
    %     nc_alphad(i,:) = noiseCorrRogCrf(pLoop);
    end
%%
F = figure('Position',[20          98        1675         903],'Name',sprintf('RhoN = %0.2f, RhoD = %0.2f',p.rho_n,p.rho_d),'Visible',flag.visible);
contData = struct(); rhonData = struct(); rhodData = struct(); rmData = struct(); epsData = struct();
contData.contrast = cont; contData.noisecorr = nc_cont;
rhonData.rho = rhoSpace; rhonData.noisecorr = nc_rhon;
rhodData.rho = rhoSpace; rhodData.noisecorr = nc_rhod;
rmData.rmax = rm; rmData.noisecorr = nc_rm;
epsData.eps = eps; epsData.noisecorr = nc_eps;
contScatter = NaN(3,2);
if flag.type=="diag"
    t = tiledlayout(1,5);
    nexttile
    hold on
    plot(cont,nc_cont,'LineWidth',2,'Color','k');
    for i=1:ncontIDs
        plot(p.cont(contID(i)),nc_cont(contID(i)),'o','MarkerSize',10,'Color',colors(i,:),'MarkerFaceColor',colors(i,:));
        contScatter(i,:) = [p.cont(contID(i)) nc_cont(contID(i))]; 
    end
    % plot(p.cont(5*ID),nc_cont(5*ID),'o','MarkerSize',10,'Color',blue,'MarkerFaceColor',blue);
    % plot(p.cont(15*ID),nc_cont(15*ID),'o','MarkerSize',10,'Color',red,'MarkerFaceColor',red);
    xlabel('Contrast')
    ylabel('noise corr.')
    axis square
    % ID = floor(nContrasts/4);
    mesh_obj = gobjects(2,3);
    line_obj = gobjects(2,3);
    rax = NaN(2,2,3);
    cax = NaN(2,2,3);
        nexttile
        hold on
        for k=1:ncontIDs
            loopID = contID(k);
            plot(rhoSpace,nc_rhon(:,loopID).','LineWidth',2,'Color',colors(k,:),'DisplayName',sprintf('Contrast = %d',p.cont(loopID)));
            rax(:,1,k) = get(gca,'YLim');
            line_obj(1,k) = gca;
        end
        xlabel('\rho_N')
    %     ylabel('noise corr.')
        axis square
    %     S = sprintf('\\fontsize{18}\\bf Contrast = %d',p.cont(loopID));
    %     ylabel({S; '\fontsize{12}\rm noise corr.'});
        nexttile
        hold on
        for k=1:ncontIDs
            loopID = contID(k);
            plot(rhoSpace,nc_rhod(:,loopID).','LineWidth',2,'Color',colors(k,:),'DisplayName',sprintf('Contrast = %d',p.cont(loopID)));
            rax(:,2,k) = get(gca,'YLim');
            line_obj(2,k) = gca;
        end
        xlabel('\rho_D');
        axis square
        nexttile
        hold on
        for k=1:ncontIDs
            loopID = contID(k);
            plot(rm,diag(nc_rm(:,:,loopID)),'LineWidth',2,'Color',colors(k,:),'DisplayName',sprintf('Contrast = %d',p.cont(loopID)));
            cax(:,1,k) = get(gca,'YLim');
            mesh_obj(1,k) = gca;
        end
        xlabel('R^{max}_1=R^{max}_2');

        axis square
        nexttile
        hold on
        for k=1:ncontIDs
            loopID = contID(k);
            plot(eps,diag(nc_eps(:,:,loopID)),'LineWidth',2,'Color',colors(k,:),'DisplayName',sprintf('Contrast = %d',p.cont(loopID)));
            cax(:,2,k) = get(gca,'Ylim');
            mesh_obj(2,k) = gca;
        end
        xlabel('\epsilon_1=\epsilon_2');
    
        % [eps,I] = sort(eps_strength);
        % plot(eps,median(nc_eps(I,:),2).');
        % xlabel('\epsilon_1 \times \epsilon_2');
        axis square
        if flag.common
            set(mesh_obj(:),'Ylim',[min(cax(:)) max(cax(:))])
            set(line_obj(:),'Ylim',[min(rax(:)) max(rax(:))])
        end
    lg = legend('NumColumns',1,'FontSize',24);
    lg.Layout.Tile = 'East';
elseif flag.type=="mesh"
%%
    t = tiledlayout(ncontIDs,5);
    nexttile([ncontIDs 1])
    hold on
    plot(p.cont,nc_cont,'LineWidth',2,'Color','k');
    for i=1:ncontIDs
        plot(p.cont(contID(i)),nc_cont(contID(i)),'o','MarkerSize',10,'Color',colors(i,:),'MarkerFaceColor',colors(i,:));
        contScatter(i,:) = [p.cont(contID(i)) nc_cont(contID(i))]; 
    end
    % plot(p.cont(5*ID),nc_cont(5*ID),'o','MarkerSize',10,'Color',blue,'MarkerFaceColor',blue);
    % plot(p.cont(15*ID),nc_cont(15*ID),'o','MarkerSize',10,'Color',red,'MarkerFaceColor',red);
    xlabel('Contrast')
    ylabel('noise corr.')
    axis square
    % ID = floor(nContrasts/4);
    mesh_obj = gobjects(2,3);
    line_obj = gobjects(2,3);
    rax = NaN(2,2,3);
    cax = NaN(2,2,3);
    for k=1:ncontIDs
        loopID = contID(k);
        nexttile
        plot(rhoSpace,nc_rhon(:,loopID).','LineWidth',2,'Color',colors(k,:));
        xlabel('\rho_N')
        rax(:,1,k) = get(gca,'YLim');
        axis square
        S = sprintf('\\fontsize{18}\\bf Contrast = %d',p.cont(loopID));
        ylabel({S; '\fontsize{12}\rm noise corr.'});
        line_obj(1,k) = gca;
        nexttile
        plot(rhoSpace,nc_rhod(:,loopID).','LineWidth',2,'Color',colors(k,:));
        xlabel('\rho_D');
        axis square
        rax(:,2,k) = get(gca,'YLim');
        line_obj(2,k) = gca;
        nexttile
        [X,Y] = meshgrid(rm);
        contourf(X,Y,nc_rm(:,:,loopID));
        colormap(gca,cmap{k});
        xlabel('R^{max}_1');
        ylabel('R^{max}_2');
        colorbar
        % [rm,I] = sort(rm_strength);
        % plot(rm,median(nc_rm(I,:),2).');
        % xlabel('R^{max}_1 \times R^{max}_2');
        axis square
        cax(:,1,k) = clim(gca);
        mesh_obj(1,k) = gca;
        nexttile
        [X,Y] = meshgrid(eps);
        contourf(X,Y,nc_eps(:,:,loopID));
        colormap(gca,cmap{k});
        xlabel('\epsilon_1');
        ylabel('\epsilon_2');
        cax(:,2,k) = clim(gca);
        mesh_obj(2,k) = gca;
        colorbar
        % [eps,I] = sort(eps_strength);
        % plot(eps,median(nc_eps(I,:),2).');
        % xlabel('\epsilon_1 \times \epsilon_2');
        axis square
    end
    if flag.common
        set(mesh_obj(:),'CLim',[min(cax(:)) max(cax(:))]);
        set(line_obj(:),'Ylim',[min(rax(:)) max(rax(:))]);
    end
else
    error("Must be 'diag' or 'mesh'")
end
title(t,sprintf('Base Parameters\nR^{max} = [%d, %d], \\epsilon = [%d, %d], \\rho_N = %0.2f, \\rho_D = %0.2f',p.Rmax,p.eps,p.rho_n,p.rho_d),'FontSize',24,'FontWeight','bold')
YLimMesh = [min(cax(:)) max(cax(:))];
YLimLine = [min(rax(:)) max(rax(:))];
if flag.save
    exportgraphics(F,sprintf('noisecorrtuning_nContrasts=%d_rhon=%s_rhod=%s_type=%s_common=%s.png',nContrasts,erase(string(rho_n),'.'),erase(string(rho_d),'.'),flag.type,string(flag.common)));
end

% if flag.exportVars
%     fileIDsBase = [...
%         sprintf("cont_tuning_rm=[%d,%d]_eps=[%d,%d]_rhon=%s_rhod=%s.dat",p.Rmax,p.eps,erase(string(rho_n),'.'),erase(string(rho_d),'.')),...
%         sprintf("rhon_tuning_rm=[%d,%d]_eps=[%d,%d]_rhon=%s_rhod=%s.dat",p.Rmax,p.eps,erase(string(rho_n),'.'),erase(string(rho_d),'.')),...
%         sprintf("rhod_tuning_rm=[%d,%d]_eps=[%d,%d]_rhon=%s_rhod=%s.dat",p.Rmax,p.eps,erase(string(rho_n),'.'),erase(string(rho_d),'.')),...
%         ];
%     for k=1:numel(fileIDsBase)
%         fID = fopen(fileIDsBase(k),"w");
% end

end
%%

function nc = noiseCorrRogCrf(p)
    [~,~,nc] = mncovRogCrf(p.pars,p.cont,p.mu_eta);
end

function p = parametersRogCrf(num_trials,num_stims,parameters)
    arguments
            num_trials (1,1) {mustBePositive} = 1e3;
            num_stims (1,1) {mustBeInteger,mustBePositive} = 10;
            parameters.Rmax (2,1) = 5+45*rand(2,1);
            parameters.eps (2,1) =  10+90*rand(2,1); % = internal_lognorm(29,16)
    
            % % Variance Parameters
            parameters.alpha_n (2,1) = max(rand(2,1),0.1);
            parameters.beta_n (2,1) = 1+rand(2,1);
            parameters.alpha_d (2,1) = max(rand(2,1),0.1);
            parameters.beta_d (2,1) = 1+rand(2,1);
    
            % % Rho Parameters
            parameters.rho_n (1,1) = -0.5+rand(1);
            parameters.rho_d (1,1) = -0.5+rand(1);
            parameters.rho_n1d2 (1,1)  = 0;
            parameters.rho_n2d1 (1,1)  = 0;
            parameters.rho (4,1) {mustBeNumeric};
    
            % % Noise Parameters
            parameters.var_eta (2,1) = zeros(2,1);
            parameters.rho_eta (1,1) = 0;
            parameters.mu_eta (2,1) = zeros(2,1);
    
            % % Function Handle
            parameters.tuning function_handle = @contrastResponse;
            
            % % Contrast Space
            parameters.cont (1,:) = 1:100;
    
            % % Override Parameters
            parameters.fit_parameters (1,19);
    end
        if isfield(parameters,'cont')
            num_stims = size(parameters.cont,2);
        end

    % % If 'cont' is not provided, then 'cont' is set using 'num_stims'
    if ~isfield(parameters,'cont')
        parameters.cont = linspace(10,100,num_stims);
    end
%     cont = parameters.cont;

    % % If 'fit_parameters' is provided, this will override all other
    % parameter settings (that are optimized)
    if isfield(parameters,'fit_parameters')
        parameters.Rmax = parameters.fit_parameters(1:2);
        parameters.eps = parameters.fit_parameters(3:4);
        parameters.alpha_n = parameters.fit_parameters([5 9]);
        parameters.beta_n  = parameters.fit_parameters([6 10]);
        parameters.alpha_d = parameters.fit_parameters([7 11]);
        parameters.beta_d = parameters.fit_parameters([8 12]);
        parameters.rho = parameters.fit_parameters(13:16);
        parameters.var_eta = parameters.fit_parameters(17:18);
        parameters.rho_eta = parameters.fit_parameters(19);
    end

    % % If 'rho' is provided, this overwrites the 'rho_X' parameters (besides 'rho_eta')
    if ~isfield(parameters,'rho')
        parameters.rho =[parameters.rho_n parameters.rho_d parameters.rho_n1d2 parameters.rho_n2d1];
    end
        p.cont = parameters.cont;
        p.mu_eta = parameters.mu_eta;
        p.pars = [parameters.Rmax(:);...
                        parameters.eps(:);...
                        parameters.alpha_n(1);...
                        parameters.beta_n(1);...
                        parameters.alpha_d(1);...
                        parameters.beta_d(1);...
                        parameters.alpha_n(2);...
                        parameters.beta_n(2);...
                        parameters.alpha_d(2);...
                        parameters.beta_d(2);...
                        parameters.rho(:);...
                        parameters.var_eta(:);...
                        parameters.rho_eta].';

end

% for i=1:numel(contID)
% %     figure("Visible","off");
%     [X,Y] = meshgrid(rmData.rmax);
%     Z = rmData.noisecorr(:,:,contID(i));
% %     X = X(1:10,1:10); Y = Y(1:10,1:10);
% %     Z = Z(1:10,1:10);
%     data = [X(:) Y(:) Z(:)];
%     data = array2table(data,'VariableNames',{'X','Y','Z'});
%     writetable(data,sprintf('rm%d.dat',contID(i)),'Delimiter','\t')
%     [X,Y] = meshgrid(epsData.eps);
%     Z = epsData.noisecorr(:,:,contID(i));
% %     X = X(1:10,1:10); Y = Y(1:10,1:10);
% %     Z = Z(1:10,1:10);
%     data = [X(:) Y(:) Z(:)];
%     data = array2table(data,'VariableNames',{'X','Y','Z'});
%     writetable(data,sprintf('eps%d.dat',contID(i)),'Delimiter','\t')
% %     meshRows = size(X,2);
% %     [C,B] = contourf(X,Y,Z);
% %     save P.dat data -ASCII
% %     levels = B.LevelList;
% end