clc
close all
load export 

%% --- FIGURE 1: cube, dist, rms, patch ---
h1 = 250;              % base height
w1 = 4 * h1;           % width = 4:1 aspect
yTop = 900;            % top position (start point)

rms_clim = [0, 0.25 * max(export.rms(:))];

figure(1);
ha = tight_subplot(1,4,[0.01 0.025],[0.05 0.05],[0.05 0.05]);

axes(ha(1));
imagesc(imrotate(export.cube,-90)); axis image off; colormap(gca,gray);
title('Cube');

axes(ha(2));
imagesc(imrotate(export.dist,-90)); axis image off; colormap(gca,"parula");
title('Distance');

axes(ha(3));
imagesc(imrotate(export.rms,-90)); axis image off; colormap(gca,gray);
title('RMS');

axes(ha(4));
imagesc(export.patch);
% imagesc(imrotate(export.patch,-90));
axis image off; colormap(gca,gray);
hold on; 
visboundaries(gca, export.bseg, 'Color', 'r', 'LineWidth', 1.5);
title('Patch');

% Place figure 1
set(gcf, 'Units', 'pixels', 'Position', [100, yTop-h1, w1, h1]);


%% --- FIGURE 2: proj2D frames 1:20 ---
frames = 1:20; 
h2 = 0.6 * h1;   % second figure height
figure(2);
set(gcf, 'Units', 'pixels', 'Position', [100, yTop-h1-h2, w1, h2])

ha = tight_subplot(2,10,[0.01 0.01],[0.05 0.05],[0.05 0.05]);

for i = 1:length(frames)
    axes(ha(i));
    imagesc(export.proj2D(:,:,frames(i)));
    % imagesc(imrotate(export.proj2D(:,:,frames(i)),-90));
    axis image off; colormap(gca,gray);
    hold on;
    visboundaries(gca, export.bseg, 'Color', 'r', 'LineWidth', 1.0);
    title('Patch');
    title(['f=' num2str(frames(i))],'FontSize',8);
end


%% --- FIGURE 3: waveforms with dual y-axes ---
h3 = 0.65 * h1;   % third figure height
figure(3);

yyaxis left
plot(export.flow,'-o','DisplayName','Flow');
ylabel('Flow');
ylim([min(export.flow) max(export.flow)]);

yyaxis right
plot(export.pc1,'-s','DisplayName','PC1');
ylabel('PC1');
ylim([min(export.pc1) max(export.pc1)]);

xlabel('Frame');
grid on;
title('Figure 3: Flow & PC1 waveforms');
legend('Flow','PC1','Location','best');

% Place figure 3 directly below figure 2
set(gcf, 'Units', 'pixels', 'Position', [100, yTop-h1-h2-h3, w1, h3])
