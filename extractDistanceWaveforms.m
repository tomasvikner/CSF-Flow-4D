%% Batch save ROI velocities for all subjects
clc; clear; close all;

parentFolder = '/Users/TXV016/Desktop/WRAP2/SUBJ/PDT';

% Get subject folders
subjectDirs = dir(parentFolder);
subjectDirs = subjectDirs([subjectDirs.isdir]); % only folders
subjectDirs = subjectDirs(~ismember({subjectDirs.name}, {'.','..'})); % remove . and ..

ilist = 1:length(subjectDirs);

% Loop over subjects
for s = ilist

    subjFolder = fullfile(parentFolder, subjectDirs(s).name);
    fprintf('\nProcessing subject: %s\n', subjectDirs(s).name);

    [foldername, subjname] = fileparts(subjFolder);
    baseparent = fileparts(foldername);

    % Locate dist and velocity data
    dmatfile = fullfile(baseparent, 'CSFmasks', subjname, 'D.mat');
    load(dmatfile, 'D');
    dist = flip(D, 2);

    cubefile = fullfile(baseparent, 'PDTreg', subjname, 'RT2.nii');
    cube = MRIread(cubefile).vol;
    cube = imrotate(cube, -90);
    cube = cube/max(cube(:));

    broi = dist > 0; 
    rvelsfolder = fullfile(subjFolder, 'rvels');
    load([rvelsfolder, '/rx.mat'], 'rx');
    load([rvelsfolder, '/ry.mat'], 'ry');
    load([rvelsfolder, '/rz.mat'], 'rz');
    rrms = sqrt(mean(double(rx).^2 + double(ry).^2 + double(rz).^2, 2));
    vrms = zeros(size(broi));
    inds = find(broi);
    vrms(inds) = rrms;
    vrms = vrms/max(vrms(:));
    vrms = vrms.*cube;

    slicemask = zeros(size(dist));
    slicemask(3:end-2, 3:end-2, 3:end-2) = 1;
    vmask = zeros(size(vrms));
    minvol = sum(broi(:)) * 0.05; 
    thresh = 1; 
    while sum(vmask(:)) < minvol % iterate 
        vmask = vrms > thresh;
        vmask = largestcomp(vmask);
        thresh = thresh - 0.01;
    end

    vinds = find(vmask);
    nvox = numel(broi);
    nframes = 20;
    vx = zeros(nvox, nframes);
    vy = zeros(nvox, nframes);
    vz = zeros(nvox, nframes);
    vx(inds, :) = rx;
    vy(inds, :) = ry;
    vz(inds, :) = rz;

    % Find inds in broi and map back to rx/ry/rz
    nx = vx(vinds, :);
    ny = vy(vinds, :);
    nz = vz(vinds, :);

    dlist = dist(vinds);

    pc1 = zeros(20, numel(dlist));
    X = [nx', ny', nz'];
    [coeff, PC1] = pca(X, 'NumComponents', 1);
    for i = 1:numel(dlist)
        ix = nx(i, :)';
        iy = ny(i, :)';
        iz = nz(i, :)';
        x = [ix, iy, iz];
        [coeff, score] = pca(x, 'NumComponents', 1);
        y = zscore(score(:, 1));
        if corr(y, PC1) < 1
            y = -y;
        end
        pc1(:, i) = y;
    end

    [sdlist, order] = sort(dlist);
    pc1 = pc1(:, order);

    % Number of bins
    nBins = 100;
    
    % Make sure distance list is sorted (you already did this)
    % dlist = distances sorted ascending
    % pc1 = [frames x voxels] matrix sorted by distance
    
    % Create bin edges
    binEdges = linspace(min(sdlist), max(sdlist), nBins+1);
    
    % Bin index for each voxel
    [~,~,binIdx] = histcounts(sdlist, binEdges);
    
    % Preallocate binned result
    nFrames = size(pc1, 1);
    bpc1 = nan(nFrames, nBins); % (frames x bins)
    bdist = nan(nBins, 1);
    for b = 1:nBins
        % Find all voxels in this bin
        voxelsInBin = (binIdx == b);
        
        if any(voxelsInBin)
            % Average across voxels in this bin
            bpc1(:, b) = zscore(mean(pc1(:, voxelsInBin), 2));
            bdist(b) = mean(sdlist(voxelsInBin));
        end
    end

    RMDEL = true;
    if RMDEL
        meanWave = zscore(PC1);
        rmse = sqrt(mean((bpc1 - meanWave).^2, 1, 'omitnan'));
        zrmse = zscore(rmse);
        outlierThresh = 3; % bins >3 SD away from mean are outliers
        outlierThresh = 2; % edit 
        validBins = zrmse <= outlierThresh;
        invalidBins = zrmse > outlierThresh;
        disp(['n-invalid: ' num2str(sum(invalidBins))]);
        bpc1_clean = bpc1(:, validBins);
        bdist_clean = bdist(validBins);
    else
        bpc1_clean = bpc1; %#ok<*UNRCH>
        bdist_clean = bdist;
    end
    
    figure;
    hold on;
    
    cmap = parula(sum(validBins));
    cmap = flipud(cmap);
    
    for b = 1:size(bpc1_clean,2)
        if ~any(isnan(bpc1_clean(:,b)))
            plot(1:nFrames, bpc1_clean(:,b), 'Color', cmap(b,:), 'LineWidth', 1.5);
        end
    end
    
    colormap(cmap);
    cb = colorbar;
    cb.Ticks = linspace(0,1,5);
    cb.TickLabels = round(linspace(min(bdist_clean), max(bdist_clean), 5),2);
    ylabel(cb,'Distance','FontSize',12);
    
    xlabel('Frame','FontSize',12);
    ylabel('Z-scored PC1','FontSize',12);
    title('Waveforms Along Distance (Outliers Removed)','FontSize',14);
    set(gca,'FontSize',12,'Box','on');

    drawnow();

    save(fullfile(subjFolder, 'bpc1.mat'), 'bpc1')
    save(fullfile(subjFolder, 'bdist.mat'), 'bdist')

end
