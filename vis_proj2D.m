clc
close all

bd = 'EXPORT/';
sn = 'wrap0523_7162_07272023';

frames2plot = 1:2:20;
% frames2plot = 1:1:20;
nf = numel(frames2plot);

if nf == 20
    ya = 10;
    ny = 100 + ya;
    nx = (ny-ya) * (nf);
    gap = [0 0]; 
    mh  = [0.05 0.005]; % bot/top
    mw  = [0.015 0.005]; % left/right
elseif nf == 10
    ya = 10;
    ny = 100 + ya;
    nx = (ny-ya) * (nf);
    gap = [0 0]; 
    mh  = [0.06 0.005]; % bot/top
    mw  = [0.025 0.005]; % left/right
end
f = figure;
f.Position = [400 600 nx 3*ny];
ha = tight_subplot(3,nf,gap,mh,mw);

rnames  = {'SC','V4','CA'};
tnames = {'SC (C2)', 'V4-CA', 'CA-V3'};
vscales = [0.2 0.8 0.9];

for r = 1:3
    
    rn = rnames{r};
    tn = tnames{r};
    vscale = vscales(r);
    fn = [bd rn '-' sn '.mat'];
    load(fn,'proj2D');
    
    vmin = min(proj2D.tvel(:));
    vmax = max(proj2D.tvel(:));
    
    for k = 1:nf
        
        axind = (r-1)*nf + k;
        axes(ha(axind));
        
        imagesc(proj2D.tvel(:,:,frames2plot(k)));
        axis image
        colormap(gca,'gray')
        hold on
        visboundaries(proj2D.bseg,'Color','r','LineWidth',0.75);
        set(gca,'XTick',[],'YTick',[]);
        
        % bottom row frame labels
        if r == 3
            xlabel(sprintf('Frame %d',frames2plot(k)),'FontSize',12);
        end
        
        % row labels
        if k == 1
            ylabel(tn,'FontWeight','bold','FontSize',12);
        end
        
        clim(vscale * [vmin vmax]);
    end
end

saveas(gcf, ['EXPORT/proj2D-SC-V4-CA.png']);




% clc
% close all
% bd = 'EXPORT/';
% sn = 'wrap0523_7162_07272023';
% % rn = 'CA'; vscale = 1;
% rn = 'SC'; vscale = 0.2;
% fn = [bd rn '-' sn '.mat'];
% 
% load(fn,'proj2D');
% 
% frames2plot = [1 3 5 7 9 11 13 15 17 19];
% % frames2plot = [3 7 9 13 17];
% % frames2plot = [3 6 9 12 15 18];
% % frames2plot = [2 4 6 8 10 12 14 16 18];
% % frames2plot = [4 9 14 19];
% nf = numel(frames2plot);
% 
% ny = 225;
% nx = (ny-35) * (nf);
% 
% f = figure;
% f.Position = [400 600 nx ny];
% gap = [0.005 0.005]; 
% mh = [0.01 0.01]; 
% mw = [0.01 0.01];
% ha = tight_subplot(1,nf,gap, mh, mw);
% 
% % axes(ha(1))
% % imagesc(proj2D.patch);
% % axis image
% % colormap(gca,'gray')
% % hold on
% % visboundaries(proj2D.bseg,'Color','r','LineWidth',0.75);
% % set(gca,'XTick',[],'YTick',[]);   % hide ticks, keep axes
% % xlabel('CUBE x VSTD','FontSize',8);
% vmin = min(proj2D.tvel(:));
% vmax = max(proj2D.tvel(:));
% for k = 1:nf
%     axes(ha(k)); %#ok<*LAXES>
%     imagesc(proj2D.tvel(:,:,frames2plot(k)));
%     axis image
%     colormap(gca,'gray')
%     hold on
%     visboundaries(proj2D.bseg,'Color','r','LineWidth',0.75);
%     set(gca,'XTick',[],'YTick',[]);   % hide ticks, keep axes
%     xlabel(sprintf('Frame %d',frames2plot(k)),'FontSize',12);
%     clim(vscale * [vmin 1.0*vmax]);
% end
% 
% % exportgraphics(gcf, ['EXPORT/proj2D-' rn '.png']);
% saveas(gcf, ['EXPORT/proj2D-' rn '.png']);
% 
% %%
% 
% clc
% close all
% 
% bd = 'EXPORT/';
% sn = 'wrap0523_7162_07272023';
% 
% frames2plot = [1 3 5 7 9 11 13 15 17 19];
% nf = numel(frames2plot);
% 
% ya = 25;
% ny = 200 + ya;
% nx = (ny-ya) * (nf);
% 
% f = figure;
% f.Position = [400 600 nx 2*ny];
% 
% gap = [0.01 0.0025]; 
% mh  = [0.01 0.01]; 
% mw  = [0.01 0.01];
% ha = tight_subplot(2,nf,gap,mh,mw);
% 
% rnames = {'SC','CA'};
% vscales = [0.2 0.95];
% 
% for r = 1:2
% 
%     rn = rnames{r};
%     vscale = vscales(r);
%     fn = [bd rn '-' sn '.mat'];
%     load(fn,'proj2D');
% 
%     vmin = min(proj2D.tvel(:));
%     vmax = max(proj2D.tvel(:));
% 
%     for k = 1:nf
% 
%         axind = (r-1)*nf + k;
%         axes(ha(axind));
% 
%         imagesc(proj2D.tvel(:,:,frames2plot(k)));
%         axis image
%         colormap(gca,'gray')
%         hold on
%         visboundaries(proj2D.bseg,'Color','r','LineWidth',0.75);
%         set(gca,'XTick',[],'YTick',[]);
% 
%         if r == 2
%             xlabel(sprintf('Frame %d',frames2plot(k)),'FontSize',12);
%         end
% 
%         clim(vscale * [vmin vmax]);
%     end
% end
% 
% saveas(gcf, ['EXPORT/proj2D-SC-CA.png']);
% 
% %%
% %%
% clc
% close all
% 
% f = figure;
% f.Position = [400 600 nx ny];
% f.Color = 'w';
% 
% % --- SC (left axis) ---
% yyaxis left
% rn = 'SC'; 
% fn = [bd rn '-' sn '.mat']; 
% load(fn,'proj2D');
% h1 = plot(proj2D.flow, 'LineWidth', 1.8);
% ylabel('SC flow rate (ml/min)');
% ax = gca;
% ax.YColor = [0 0.4470 0.7410];      % blue
% h1.Color  = [0 0.4470 0.7410];      % match line to axis
% 
% % --- CA (right axis) ---
% yyaxis right
% rn = 'CA'; 
% fn = [bd rn '-' sn '.mat']; 
% load(fn,'proj2D');
% h2 = plot(proj2D.flow, '--', 'LineWidth', 1.8);
% ylabel('CA flow rate (ml/min)');
% ax.YColor = [0.8500 0.3250 0.0980]; % orange
% h2.Color  = [0.8500 0.3250 0.0980];  % match line to axis
% 
% % --- Styling ---
% grid on
% box off
% set(gca, 'TickDir', 'out');
% axis tight
% xlabel('Frame')
% 
% legend([h2 h1], ...
%     {'CA flow rate (ml/min)', 'SC flow rate (ml/min)'}, ...
%     'Location','best', ...
%     'Box','off');
% 
% exportgraphics(gcf, ['EXPORT/projWFx2.png'], 'BackgroundColor','white');
% 
