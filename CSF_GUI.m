function CSF_GUI_V6

addpath('Freesurfer Tools')

delete(findall(0, 'Type', 'figure', 'Name', '4D CSF Flow Viewer'));

ncols = 4;
[fig, g, ax, s_frame, s_slice] = setupGUI(@loadData, @updateDisplays, @keyControl, ncols);

% Add separate slider for ax(1,4) (Flow Plane Viewer)
s_frame_plane = uislider(g, ...
    'Limits', [1 10], ...
    'Value', 1, ...
    'MajorTicks', [], ...
    'Orientation', 'vertical');
s_frame_plane.Layout.Row = 3;
s_frame_plane.Layout.Column = ncols+1;  % adjust layout if needed

addlistener(s_frame_plane, 'ValueChanged', @(src, evt) updateFlowPlane());

% Create dropdown and position it over ax(1,4)
toggle = uidropdown(fig, ...
    'Items', {'proj', 'vx', 'vy', 'vz'}, ...
    'Value', 'proj', ...
    'Position', [830, 865, 55, 30]);

szDropdown = uidropdown(fig, ...
    'Items', {'15', '25', '35', '45', '55', '65', '75'}, ...
    'Value', '25', ...
    'Position', [830, 570, 55, 30], ...
    'Tooltip', 'Patch size');
thrDropdown = uidropdown(fig, ...
    'Items', {'20', '30', '40', '50', '60', '70', '80'}, ...
    'Value', '50', ...
    'Position', [1065, 570, 55, 30], ...
    'Tooltip', 'Patch size');
thrDropdown.ValueChangedFcn = @(src, evt) updateFlowPlane();
szDropdown.ValueChangedFcn = @(src, evt) onPatchOrThreshChange();
thrDropdown.ValueChangedFcn = @(src, evt) onPatchOrThreshChange();
function onPatchOrThreshChange()
    updateDisplays();
    updateFlowPlane();
    % Update waveforms only if clicked point exists
    if ~isempty(clickedX) && ~isempty(clickedY)
        updateWaveformsFromCoords(clickedX, clickedY, round(s_slice.Value));
    end
end

playBtn = uibutton(fig, ...
    'Text', 'Loop', ...
    'Position', [1065, 870, 55, 30]);
playBtn.ButtonPushedFcn = @(btn, evt) togglePlay();

% Add Set Dir / Reset Dir toggle button
setDirBtn = uibutton(fig, ...
    'Text', 'Set Dir', ...
    'Position', [555, 570, 80, 30], ...
    'ButtonPushedFcn', @(btn,event) toggleSetDir());

direction = [];     % current flow direction (from PCA or manual)
manualdir = [];  % Global manual direction override
manualMode = false; % Flag manual direction active

ROIPCA = false; % Flag to start direction estimation based on voxel PCA
VOLPCA = false; % Flag to decide whether to sample voxel or volume (eg LVs) PCA waveform
volumeCoords = []; % Store location about currently sampled volume PCA

% Add voxel/ROI PCA toggle button for direction estimation
setDirBtn = uibutton(fig, ...
    'Text', 'ROI PCA', ...
    'Position', [750, 585, 75, 30], ...
    'ButtonPushedFcn', @(btn,event) toggleDirVoxROI());

setPCABtn = uibutton(fig, ...
    'Text', 'Vol PCA', ...
    'Position', [750, 615, 75, 30], ...
    'ButtonPushedFcn', @(btn,event) togglePCAVoxVol());

resetGUIbutton = uibutton(fig, ...
    'Text', 'Reset GUI', ...
    'Position', [10 570 100 30], ...
    'ButtonPushedFcn', @(btn, evt) resetGUI());

popWFsBtn = uibutton(fig, ...
    'Text', 'Pop WFs', ...
    'Position', [750, 645, 75, 30], ...
    'ButtonPushedFcn', @(btn,event) togglePopWFs());

exportFiguresBtn = uibutton(fig, ...
    'Text', 'Export', ...
    'Position', [10 540 100 30], ...
    'ButtonPushedFcn', @(btn, evt) exportFigures());

% Create button panel for saving waveforms (row 3, col 4)
btnPanel = uipanel(g);
btnPanel.Layout.Row = 3;
btnPanel.Layout.Column = 4;

btnGrid = uigridlayout(btnPanel, [8,1]);
btnGrid.RowHeight = repmat({'1x'}, 1, 8);
btnGrid.ColumnWidth = {'1x'};

btnNames = {'Left LV', 'Right LV', 'Left FMo', 'Right FMo', '3rd Vent', 'C. Aqueduct', '4th Vent', 'Spinal Canal'};
fileNames = {'LLV.mat', 'RLV.mat', 'LFMo.mat', 'RFMo.mat', 'V3.mat', 'CA.mat', 'V4.mat', 'SC.mat'};

btnHandles = gobjects(length(btnNames), 1);
for k = 1:length(btnNames)
    btnHandles(k) = uibutton(btnGrid, 'Text', btnNames{k}, ...
        'ButtonPushedFcn', @(btn,event) saveWaveforms(fileNames{k}));
end

    function keyControl(~, event)
        step = 1;

        switch event.Key
            case 'rightarrow'
                if s_frame.Value < s_frame.Limits(2)
                    s_frame.Value = s_frame.Value + step;
                    updateDisplays();
                end
            case 'leftarrow'
                if s_frame.Value > s_frame.Limits(1)
                    s_frame.Value = s_frame.Value - step;
                    updateDisplays();
                end
            case 'uparrow'
                if s_slice.Value < s_slice.Limits(2)
                    s_slice.Value = s_slice.Value + step;
                    updateDisplays();
                end
            case 'downarrow'
                if s_slice.Value > s_slice.Limits(1)
                    s_slice.Value = s_slice.Value - step;
                    updateDisplays();
                end
        end
    end

addlistener(s_frame,'ValueChanged',@(src,evt) updateDisplays());
addlistener(s_slice,'ValueChanged',@(src,evt) updateDisplays());

data = struct();

isPlaying = false;
playTimer = [];

marker = [];
flowArrow = [];

clickedX = [];
clickedY = [];

viewplane = 'sagittal';

vel_clim = [];
rms_clim = [];

subjectFolder = ''; % To store folder path from loadData
savefolder = '';
WF = [];

% *** LOAD DATA ***
    function loadData()
        basefolder = uigetdir();
        if basefolder == 0, return; end

        subjectFolder = basefolder;  % Save basefolder for saving waveforms
        savefolder = fullfile(subjectFolder, 'Waveforms');
        if ~exist(savefolder, 'dir')
            mkdir(savefolder);
        end

        % Extract folders for RT2 and D.mat
        [foldername, subjname] = fileparts(subjectFolder);
        baseparent = fileparts(foldername);

        cubefile = fullfile(baseparent, 'PDTreg', subjname, 'RT2.nii');
        dmatfile = fullfile(baseparent, 'CSFmasks', subjname, 'D.mat');
        load(dmatfile, 'D');

        files = dir(fullfile(basefolder, '**', '*'));
        vxfile = ''; vyfile = ''; vzfile = ''; vbgfile = ''; magfile = '';
        for i = 1:length(files) 
            if files(i).isdir, continue; end
            fname = files(i).name;
            fpath = fullfile(files(i).folder, fname);

            if contains(fname, 'CSF_FLOW_-_LR') && endsWith(fname, '.nii.gz')
                vxfile = fpath;
            elseif contains(fname, 'CSF_FLOW_-_AP') && endsWith(fname, '.nii.gz')
                vyfile = fpath;
            elseif contains(fname, 'CSF_FLOW_-_SI') && endsWith(fname, '.nii.gz')
                vzfile = fpath;
            elseif strcmp(fname, 'VBG.nii')
                vbgfile = fpath;
            elseif strcmp(fname, 'MAG.nii')
                magfile = fpath;
            end
        end

        if isempty(vxfile) || isempty(vyfile) || isempty(vzfile)
            uialert(fig, 'Missing required NIfTI files.', 'Load Error');
            return;
        end

        data.vx = MRIread(vxfile).vol;
        data.vy = MRIread(vyfile).vol;
        data.vz = MRIread(vzfile).vol;
        data.mag = MRIread(magfile).vol;
        data.cube = MRIread(cubefile).vol;

        if ~isempty(vbgfile)
            disp('VBG file found')
            vbg = double(niftiread(vbgfile));
            if size(vbg, 1) == 3
                disp('Applying BG correction')
                data.vx = data.vx - squeeze(vbg(1,:,:,:));
                data.vy = data.vy - squeeze(vbg(2,:,:,:));
                data.vz = data.vz - squeeze(vbg(3,:,:,:));
            else
                disp('VBG matrix wrong size')
            end
        end

        data.vx = imrotate(data.vx, -90);
        data.vy = imrotate(data.vy, -90);
        data.vz = imrotate(data.vz, -90);
        data.mag = imrotate(data.mag, -90);
        data.cube = imrotate(data.cube, -90);
        data.dist = flip(D, 2);

        data.rms = sqrt(mean(double(data.vx).^2 + double(data.vy).^2 + double(data.vz).^2, 4));
        data.mixed = data.rms.*data.cube; 

        % Automatically extract reference WF in bottom of SC 
        nt = size(data.vx, 4); nvox = numel(data.vx)/nt;
        rx = reshape(data.vx, nvox, nt);
        ry = reshape(data.vy, nvox, nt);
        rz = reshape(data.vz, nvox, nt);
        SCinds = extractSCinit(data.dist);
        X = [rx(SCinds, :); ry(SCinds, :); rz(SCinds, :)]';
        [COEFF, SCORE] = pca(X); %#ok<*ASGLU> 
        data.SCWF = SCORE(:, 1);

        CROPON = true; % reduce from ~21 million to ~9 million voxels 
        if CROPON
            data = cropdata(data);
        end
    
        % Set display limits
        vmax = 0.8 * max([max(abs(data.vx(:))), max(abs(data.vy(:))), max(abs(data.vz(:)))]);
        vel_clim = [-vmax, vmax];
        rms_clim = [0, 0.25 * max(data.rms(:))];

        [~, ~, zres, nframes] = size(data.vx);
        s_frame.Limits = [1 nframes];
        s_slice.Limits = [1 zres];
        s_frame.Value = 1;
        s_slice.Value = round(zres/2);

        s_frame_plane.Limits = [1 nframes];
        s_frame_plane.Value = 1;

        updateDisplays();
        updateFlowPlane(); 

        disp('LoadData complete.')

    end

    function updateDisplays()
        if isempty(fieldnames(data)), return; end
        frame = round(s_frame.Value);
        slice = round(s_slice.Value);

        plotVelocitySlice(ax(1,1), data.vx, slice, frame, vel_clim, 'vx');
        plotVelocitySlice(ax(1,2), data.vy, slice, frame, vel_clim, 'vy');
        plotVelocitySlice(ax(1,3), data.vz, slice, frame, vel_clim, 'vz');
        imagesc(ax(2,1), squeeze(data.dist(:,:,slice))'); axis(ax(2,1), 'image'); title(ax(2,1), 'CSF dist.'); colormap(ax(2,1), hot);
        imagesc(ax(2,2), squeeze(data.cube(:,:,slice))'); axis(ax(2,2), 'image'); title(ax(2,2), 'T2 cube')

        if ~isempty(marker)
            % Only delete valid handles that still exist
            validMarkers = marker(ishandle(marker));
            if ~isempty(validMarkers)
                delete(validMarkers);
            end
            marker = [];
        end

        prevLim = axis(ax(2,3));  % Save zoom state
        cla(ax(2,3));
        imagesc(ax(2,3), squeeze(data.rms(:,:,slice))');
        if ~isequal(prevLim, [0 1 0 1])  % If zoomed, restore view
            axis(ax(2,3), prevLim);
        else
            axis(ax(2,3), 'image');
        end

        % Allow pointer to follow slice scroll
        if ~isempty(clickedX) && ~isempty(clickedY)
            updateWaveformsFromCoords(clickedX, clickedY, slice);
        end

        colormap(ax(2,3), gray);
        title(ax(2,3), 'RMS velocity');
        caxis(ax(2,3), rms_clim);

        % Initiate ax(3,3) as clickable
        ax(3,3);

        clickableAxes = [ax(1,1), ax(1,2), ax(1,3), ax(2,1), ax(2,2), ax(2,3), ax(3,3)];
        for k = 1:numel(clickableAxes) %#ok<*FXUP> 
            axx = clickableAxes(k);
            axx.PickableParts = 'all';
            axChildren = axx.Children;
            if ~isempty(axChildren)
                axChildren(end).HitTest = 'on';
                axChildren(end).ButtonDownFcn = @(src, event) updateWaveformsFromClick(event, slice);
            end
        end

        fig.Pointer = 'crosshair';

        % Show time-resolved PCA-aligned 2D plane (proj2D) in ax(1,4)
        if isfield(data, 'proj2D') && ~isempty(data.proj2D)
            frame = round(s_frame.Value);
            imagesc(ax(1,4), data.proj2D(:,:,frame));
            axis(ax(1,4), 'image');
            colormap(ax(1,4), gray);
            d2d = data.proj2D .* data.bseg; 
            caxis(ax(1,4), [0.8 * min(d2d(:)), 0.8 * max(d2d(:))]);
            title(ax(1,4), sprintf('Flow plane (%s), frame %d', toggle.Value, frame));
            hold(ax(1,4), 'on'); 
            visboundaries(ax(1,4), data.bseg, 'Color', 'r', 'LineWidth', 1.5); 
            hold(ax(1,4), 'off');
        end

    end

    function updateFlowPlane()
        if ~isfield(data, 'proj2D') || isempty(data.proj2D), return; end
        frame = round(s_frame_plane.Value);

        switch toggle.Value
            case 'proj'
                img = data.proj2D(:,:,frame);
            case 'vx'
                img = data.velx2D(:,:,frame);
            case 'vy'
                img = data.vely2D(:,:,frame);
            case 'vz'
                img = data.velz2D(:,:,frame);
        end
    
        imagesc(ax(1,4), img);
        axis(ax(1,4), 'image');
        colormap(ax(1,4), gray);
        title(ax(1,4), sprintf('Flow plane (%s), frame %d', toggle.Value, frame));
        hold(ax(1,4), 'on'); 
        visboundaries(ax(1,4), data.bseg, 'Color', 'r', 'LineWidth', 1.5); 
        hold(ax(1,4), 'off');
    end

    function saveWaveforms(filename)
        if isempty(subjectFolder)
            uialert(fig, 'No subject folder loaded. Load data first.', 'Save Error');
            return;
        end
        if isempty(marker)
            uialert(fig, 'Select a point on RMS image first.', 'Save Error');
            return;
        end

        if isempty(clickedX) || isempty(clickedY)
            uialert(fig, 'No click selected for waveforms.', 'Save Error');
            return;
        end

        x = clickedX;
        y = clickedY;
        z = round(s_slice.Value);

        [xres, yres, zres, ~] = size(data.vx);
        if any([x y z] < 1) || x > xres || y > yres || z > zres
            uialert(fig, 'Invalid position for waveforms.', 'Save Error');
            return;
        end

        % vx_t = squeeze(data.vx(x, y, z, :));
        % vy_t = squeeze(data.vy(x, y, z, :));
        % vz_t = squeeze(data.vz(x, y, z, :));
        % WF.vel.x = vx_t;
        % WF.vel.y = vy_t;
        % WF.vel.z = vz_t;

        WF.vel.x = data.velx2D;
        WF.vel.y = data.vely2D;
        WF.vel.z = data.velz2D;
        WF.proj = data.proj2D;
        WF.flow = data.flow;
        WF.pc1 = data.pc1;
        if isfield(data, 'PC1')
            WF.PC1 = data.PC1;
        end

        WF.coords.x = x; 
        WF.coords.y = y;
        WF.coords.z = z;
        WF.dist = data.dist(x, y, z);

        % TODO: save ROIs

        save(fullfile(savefolder, filename), 'WF');
        disp(['Saved waveforms to ', fullfile(savefolder, filename)]);
    end

    % TODO: when clicking in coronal
    function updateWaveformsFromClick(event, slice)
        pos = round(event.IntersectionPoint(1:2));

        % Currently NOT working, seems annoying to and low prio
        if strcmp(viewplane, 'cor') % gotta modify s_slice? 
            % clickedY = pos(2);
            % slice = pos(1);
            % s_slice = slice;
            % y = clickedY;
            % disp(' *** TEMP: updating from coronal *** ')
        else % else assume sagittal 
            clickedX = pos(1);
            clickedY = pos(2);
            x = clickedX;
            y = clickedY;
        end
        updateWaveformsFromCoords(x, y, slice);
    end

    function updateVolumeBounds(ai1, ai2, seg2D)
        if VOLPCA && ~isempty(volumeCoords)

            % Delete old boundary if it exists
            old = findobj(ax(ai1, ai2), 'Tag', 'segBoundary');
            delete(old);
    
            % Draw new boundary
            hb = visboundaries(ax(ai1, ai2), volumeCoords.(seg2D), 'Color', 'r', 'LineWidth', 1.5);
            set(hb, 'Tag', 'segBoundary');
        end
    end

    function updateWaveformsFromCoords(x, y, z)

        [xres, yres, zres, nframes] = size(data.vx);
        if any([x y z] < 1) || x > xres || y > yres || z > zres
            return;
        end
        
        % Extract time series
        vx_t = squeeze(data.vx(x, y, z, :));
        vy_t = squeeze(data.vy(x, y, z, :));
        vz_t = squeeze(data.vz(x, y, z, :));

        plotWaveforms(ax(3,1), vx_t, vy_t, vz_t);

        % Perform PCA to get CSF flow direction
        V = [vx_t vy_t vz_t];  % (frames x 3)
        [coeff, score, ~] = pca(V, 'NumComponents', 1);
        pc1 = score(:,1);  % Time series
        if manualMode && ~isempty(manualdir)
            direction = manualdir(:);  % use manualdir for plane extraction etc.
        else
            direction = coeff(:,1);    % PCA direction loadings
            manualdir = [];            % clear manualdir if any
        end

        % patch_width = 17; 
        patch_width = str2double(szDropdown.Value);
        local_thresh = str2double(thrDropdown.Value);
        % MODE = 'rms';
        % MODE = 'cube';
        MODE = 'mixed';
        [flow, ~, patch, bseg, projV, proj2D] = extractThroughPlaneFlow_interp2(data, [x, y, z], direction, patch_width, MODE, local_thresh);
        data.proj2D = proj2D.proj;
        data.velx2D = proj2D.velx;
        data.vely2D = proj2D.vely;
        data.velz2D = proj2D.velz;
        data.patch = patch; % this is actually patch_interp
        data.bseg = bseg;

        % Edit for ROI based direction estimation; only if proj2D.velx2D exists for coord 
        if ~manualMode && ROIPCA 
            inds = find(data.bseg);
            [nx, ny, nt] = size(data.velx2D);
            rx = reshape(data.velx2D, nx*ny, nt);
            ry = reshape(data.vely2D, nx*ny, nt);
            rz = reshape(data.velz2D, nx*ny, nt);
            rx = mean(rx(inds, :));
            ry = mean(ry(inds, :));
            rz = mean(rz(inds, :));
            X = [rx' ry' rz'];  % (frames x 3 x Voxels in ROI)
            [coeff, score, ~] = pca(X, 'NumComponents', 1);
            pc1 = score(:,1);  % Time series
            direction = coeff;
            [flow, patch, patch, bseg, projV, proj2D] = extractThroughPlaneFlow_interp2(data, [x, y, z], direction, patch_width, MODE, local_thresh);
            data.proj2D = proj2D.proj;
            data.velx2D = proj2D.velx;
            data.vely2D = proj2D.vely;
            data.velz2D = proj2D.velz;
            data.bseg = bseg;
        end

        data.flow = flow;
        data.pc1 = pc1; 

        % Clear old arrow if it exists
        if isgraphics(flowArrow)
            delete(flowArrow);
            flowArrow = [];
        end
        dir_xy = -direction([2, 3]);  % vy (x), vz (z)
        arrowLength = 15;
        if norm(dir_xy) > 0
            % dir_xy = dir_xy / norm(dir_xy) * arrowLength;
            dir_xy = dir_xy / norm(direction) * arrowLength;
        else
            dir_xy = [0 arrowLength];
        end
        if ~isempty(marker)
            delete(marker(ishandle(marker)));
        end
        marker = [];
        for j = 1:3, hold(ax(2,j), 'on'); end
        marker = [ ...
            quiver(ax(2,1), x, y, dir_xy(1), dir_xy(2), 'Color', 'b', 'LineWidth', 2.5, 'AutoScale', 'off', 'MaxHeadSize', 2.5); 
            quiver(ax(2,2), x, y, dir_xy(1), dir_xy(2), 'Color', 'b', 'LineWidth', 2.5, 'AutoScale', 'off', 'MaxHeadSize', 2.5); 
            quiver(ax(2,3), x, y, dir_xy(1), dir_xy(2), 'Color', 'b', 'LineWidth', 2.5, 'AutoScale', 'off', 'MaxHeadSize', 2.5)];
        for j = 1:3, hold(ax(2,j), 'off'); end

        % Local CS 
        imagesc(ax(2,4), patch);
        axis(ax(2,4), 'image'); 
        title(ax(2,4), ['PCA-aligned plane: ' MODE]);
        colormap(ax(2,4), gray); 
        hold(ax(2,4), 'on'); 
        visboundaries(ax(2,4), data.bseg, 'Color', 'r', 'LineWidth', 1.5); 
        hold(ax(2,4), 'off');

        % Clear and prepare left axis 
        cla(ax(3,2), 'reset');
        yyaxis(ax(3,2), 'left');
        hold(ax(3,2), 'on');
        plot(ax(3,2), zscore(data.SCWF), 'k--', 'LineWidth', 1.5);

        % Plot Volume or Voxel PCA depending on mode
        ndils = 8; % 4 mm (10 slow, might need for LVs tho)
        nrodes = 2; % 1 mm 
        if VOLPCA
            [PC1, volumeCoords] = extractPCAFromVols(data, x, y, z, ndils, nrodes);
            for ai = 1:3
                updateVolumeBounds(2, ai, 'seg2D');
            end
            if corr(PC1, data.SCWF) < 0
                PC1 = - PC1;
            end
            data.PC1 = PC1;
            plot(ax(3,2), zscore(PC1), 'k-', 'LineWidth', 1.5);
            ylabel(ax(3,2), 'Volume PCA (z-score)');
        else
            if corr(pc1, data.SCWF) < 0
                pc1 = - pc1;
            end
            data.pc1 = pc1;
            plot(ax(3,2), zscore(pc1), 'k-', 'LineWidth', 1.5); %#ok<*UNRCH> 
            ylabel(ax(3,2), 'Voxel PCA (z-score)');
        end

        % Right Y-axis: Flow through plane
        yyaxis(ax(3,2), 'right');
        hold(ax(3,2), 'on'); % need to reapply hold on for each axis? 
        plot(ax(3,2), flow, 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5);
        ylabel(ax(3,2), 'Flow rate (ml/s)');

        % Final touches
        title(ax(3,2), sprintf('PC1 dir: [%.2f  %.2f  %.2f]', direction(1), direction(2), direction(3)));
        xlabel(ax(3,2), 'Frame');
        if VOLPCA
            legend(ax(3, 2), 'SCinit PCA', 'Volume PCA', 'CS flow rate')
        else
            legend(ax(3, 2), 'SCinit PCA', 'Voxel PCA', 'CS flow rate')
        end

        % Plot coronal RMS view in ax(3,3)
        cla(ax(3,3), 'reset');
        cor = squeeze(data.cube(clickedX, :, :));  
        % cor = squeeze(data.mixed(clickedX, :, :));  
        nz = size(cor, 2);
        ys = clickedY-nz/2; ye = clickedY+nz/2;
        if ye > size(data.cube, 2)
            ye = size(data.cube, 2);
        end
        cor = cor(ys:ye, :);
        imagesc(ax(3,3), cor);
        axis(ax(3,3), 'image');
        colormap(ax(3,3), gray);
        title(ax(3,3), 'CUBE \times RMS (coronal)');
        caxis(ax(3,3), [0, 0.9 * max(data.cube(:))]); 
        hold(ax(3,3), 'on');

        % Coronal view boundary 
        if VOLPCA
            cvslice = squeeze(volumeCoords.cvol(clickedX, :, :));
            cvbound = cvslice(ys:ye, :);
            old = findobj(ax(3, 3), 'Tag', 'segBoundary');
            delete(old);
            hb = visboundaries(ax(3, 3), cvbound, 'Color', 'r', 'LineWidth', 1.5); % Draw new boundary
            set(hb, 'Tag', 'segBoundary');
        end
        
        % Add velocity direction as quiver
        dir_xz = direction([1, 3]);  % vx and vz components
        dir_xz = dir_xz / norm(direction) * 10;  % scale for visibility
        
        quiver(ax(3,3), z, (nz/2) - 1, dir_xz(1), -dir_xz(2), ...
            'Color', 'b', 'LineWidth', 2.5, 'AutoScale', 'off', 'MaxHeadSize', 2.5);
        
        hold(ax(3,3), 'off');

    end

    function togglePlay()
        isPlaying = ~isPlaying;
        if isPlaying
            playBtn.Text = 'Stop';
            playTimer = timer( ...
                'ExecutionMode', 'fixedRate', ...
                'Period', 0.2, ...
                'TimerFcn', @(~,~) nextFrame());
            start(playTimer);
        else
            playBtn.Text = 'Play';
            if isvalid(playTimer)
                stop(playTimer);
                delete(playTimer);
            end
            playTimer = [];
        end
    end

    function nextFrame()
        if isempty(data) || ~isfield(data, 'proj2D') || isempty(data.proj2D), return; end
        curr = round(s_frame_plane.Value);
        maxFrame = s_frame_plane.Limits(2);
        next = mod(curr, maxFrame) + 1;  % Loop back to 1
        s_frame_plane.Value = next;
        updateFlowPlane();
    end

    fig.CloseRequestFcn = @(src, event) closeGUI();
    function closeGUI()
        if isvalid(fig)
            delete(fig);
        end
        if ~isempty(playTimer) && isvalid(playTimer)
            stop(playTimer);
            delete(playTimer);
        end
        clear global data;  % optional
        clc;
    end

    function resetGUI()
        closeGUI();    % Call the same cleanup
        clear all; close all force; %#ok<CLALL> 
        pack;  % force memory defragmentation (may take time)
        CSF_GUI_V6();  % Restart GUI
    end

    function exportFigures()
        export = data;
        export = rmfield(export, 'vx');
        export = rmfield(export, 'vy');
        export = rmfield(export, 'vz');
        slice = round(s_slice.Value);
        export.mag = data.mag(:, :, slice);
        export.cube = data.cube(:, :, slice);
        export.rms = data.rms(:, :, slice);
        export.mixed = data.mixed(:, :, slice);
        export.dist = data.dist(:, :, slice);
        export.proj2D = data.proj2D;
        export.velx2D = data.velx2D;
        export.velx2D = data.vely2D;
        export.velx2D = data.velx2D;
        export.patch = data.patch;
        export.patch_interp = data.patch_interp;
        save('export.mat', 'export');
    end

    function toggleSetDir()
        manualMode = ~manualMode;
        if manualMode
            setDirBtn.Text = 'Reset Dir';
            manualdir = direction;  % Initialize manualdir from current PCA direction
            openDirectionPopup();
            updateDirection(manualdir);
        else
            setDirBtn.Text = 'Set Dir';
            closeDirectionPopup();
            updateDirection(direction);  % revert to PCA direction
        end
    end
    
    popupFigDir = [];  % persistent handle for popup window

    popupFigWFs = [];
    scatterHandles = [];
    
    function togglePopWFs()

        if ~isempty(popupFigWFs) && isvalid(popupFigWFs)
            % Delete scatter points
            if ~isempty(scatterHandles) && all(ishandle(scatterHandles))
                delete(scatterHandles);
                scatterHandles = [];
            end
            % Close the figure
            delete(popupFigWFs);
            popupFigWFs = [];
            return;
        end
    
        % Load waveform files
        d = dir(fullfile(savefolder, '*.mat'));
        popupFigWFs = figure('Name','Waveforms', 'Position', [590 535 340 200]);
        hold on;
    
        y_coords = zeros(1, numel(d));
        z_coords = zeros(1, numel(d));
    
        for i = 1:numel(d)
            wn = fullfile(savefolder, d(i).name);
            load(wn, 'WF');
            plot(WF.pc1); % TODO
            % plot(WF.flow); % TODO
    
            y_coords(i) = WF.coords.x;
            z_coords(i) = WF.coords.y;
        end
    
        legend({d.name}, 'Interpreter', 'none');
    
        % Plot scatter on ax(2,2) without clearing image
        axes(ax(2,2)); 
        hold(ax(2,2), 'on');
        scatterHandles = scatter(ax(2,2), y_coords, z_coords, 40, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r');
    
        % Set close callback on popup figure to remove scatter points
        set(popupFigWFs, 'CloseRequestFcn', @onPopupClose);
    
        function onPopupClose(src, ~)
            if ~isempty(scatterHandles) && all(ishandle(scatterHandles))
                delete(scatterHandles);
                scatterHandles = [];
            end
            delete(src);
            popupFigWFs = [];
        end
    end
    
    function openDirectionPopup()
        % Ensure popupFigDir is shared (define as persistent or outside this function in your GUI)
        % persistent popupFigDir
        if ~isempty(popupFigDir) && isvalid(popupFigDir)
            figure(popupFigDir); % bring to front
            return;
        end
    
        popupFigDir = uifigure('Name','Manual Direction Control', 'Position', [590 735 340 200]);
    
        % Rotate direction buttons
        btnLabels = {'-x', '+x', '-y', '+y', '-z', '+z'};
        btnPositions = [20 140 60 40; 100 140 60 40;
                        20 80 60 40; 100 80 60 40;
                        20 20 60 40; 100 20 60 40];
    
        for k = 1:6
            uibutton(popupFigDir, 'Text', btnLabels{k}, 'Position', btnPositions(k,:), ...
                'ButtonPushedFcn', @(btn,~) rotateManualDir(btn.Text));
        end
    
        % Set-direction buttons
        setLabels = {'X', '-X', 'Y', '-Y', 'Z', '-Z'};
        setDirs = [1 0 0; -1 0 0; 0 1 0; 0 -1 0; 0 0 1; 0 0 -1];
        setPositions = [180 140 60 40; 260 140 60 40;
                        180 80 60 40; 260 80 60 40;
                        180 20 60 40; 260 20 60 40];
    
        for k = 1:6
            dirVec = setDirs(k,:); % preserve loop variable
            uibutton(popupFigDir, 'Text', setLabels{k}, 'Position', setPositions(k,:), ...
                'ButtonPushedFcn', @(~,~) setManualDirection(dirVec));
        end
    end

    function closeDirectionPopup()
        if ~isempty(popupFigDir) && isvalid(popupFigDir)
            close(popupFigDir);
            popupFigDir = [];
        end
    end

    function rotateManualDir(directionStr)
        angleRad = deg2rad(5);
        switch directionStr
            case '-x', R = rotX(-angleRad);
            case '+x', R = rotX(angleRad);
            case '-y', R = rotY(-angleRad);
            case '+y', R = rotY(angleRad);
            case '-z', R = rotZ(-angleRad);
            case '+z', R = rotZ(angleRad);
            otherwise, R = eye(3);
        end
        manualdir = manualdir(:);
        manualdir = R * manualdir; % Rotate 
        manualdir = manualdir / norm(manualdir); 
        updateDirection(manualdir);  % Update and trigger downstream updates (waveforms, flowplane)
    end

    function setManualDirection(vec)
        vec = vec(:);
        if norm(vec) == 0
            return;
        end
        manualdir = vec / norm(vec);
        manualMode = true;
        updateDirection(manualdir);
    end
    
    function updateDirection(newdir)
        direction = newdir(:);
        direction = direction / norm(direction);
    
        % If you have a valid clicked position, update waveforms & quivers
        if ~isempty(clickedX) && ~isempty(clickedY) && ~isempty(direction)
            % Round slice index
            zslice = round(s_slice.Value);
            % This will update quivers, patch, waveforms, etc.
            updateWaveformsFromCoords(clickedX, clickedY, zslice); 
        end
    
        % Update flow plane view as well
        updateFlowPlane();
    end

    function R = rotX(theta)
        R = [1 0 0;
             0 cos(theta) -sin(theta);
             0 sin(theta) cos(theta)];
    end
    
    function R = rotY(theta)
        R = [cos(theta) 0 sin(theta);
             0 1 0;
             -sin(theta) 0 cos(theta)];
    end
    
    function R = rotZ(theta)
        R = [cos(theta) -sin(theta) 0;
             sin(theta) cos(theta) 0;
             0 0 1];
    end

    function toggleDirVoxROI()
        ROIPCA = ~ROIPCA;
        if ROIPCA
            setDirBtn.Text = 'Vox PCA';
        else
            setDirBtn.Text = 'ROI PCA';
        end
        zslice = round(s_slice.Value);
        updateWaveformsFromCoords(clickedX, clickedY, zslice); 
    end

    function togglePCAVoxVol()
        VOLPCA = ~VOLPCA;
        if VOLPCA
            setPCABtn.Text = 'Vox PCA';
        else
            setPCABtn.Text = 'Vol PCA';
        end
        zslice = round(s_slice.Value);
        updateWaveformsFromCoords(clickedX, clickedY, zslice); 
    end

% Link the axes limits for synchronized zoom/pan
% linkaxes([ax(1,1), ax(1,2), ax(1,3), ax(2,1), ax(2,2), ax(2,3), ax(2,4)], 'xy');

end
