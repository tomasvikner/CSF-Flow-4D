% Load local vels in brainmask if LOCAL=true 
% Save to server subject folders (BV_Waveforms)
% Circularity constraint for CA (extractCircular)
% - Should adapt this for larger ROIs, e.g., SC 

function CSF_GUI_brain

% USER SETTINGS: EDIT BEFORE RUN
LOCAL = true; 
LOCALVELS = '/Users/txv016/Documents/BRAINVELS'; 

% Directory for saving waveforms in remote subject folder 
% OUTFOLDER = 'BV_Waveforms'; % WRAP 2025 U/L CA 
OUTFOLDER = 'TEMP';

% PATH ON REMOTE SERVER
BASEPATH = '/Volumes/groups/CVMRIGroup/Users/txv016/WRAP2/niis/CURRENT/';

currentFrame = 1; % Track current frame
maxFrame = 20;
scatterHandles = [];

% TEMP 
addpath('/Users/txv016/Documents/MATLAB/Freesurfer Tools')

delete(findall(0, 'Type', 'figure', 'Name', '4D CSF Flow Viewer'));

ncols = 4;
fig = uifigure('Name','4D CSF Flow Viewer','Position',[100 100 300*ncols 900]);
g = uigridlayout(fig,[3, ncols+1]);
g.RowHeight = {'1x','1x','1x'};
widths = repmat({sprintf('%gx', 1)}, 1, ncols);
g.ColumnWidth = [widths, {80}];
g.RowSpacing = 2.5;       % Space between rows
g.ColumnSpacing = 2.5;    % Space between columns
g.Padding = [2.5 2.5 2.5 2.5];  % Space around entire grid

ax = gobjects(3, ncols);
for i = 1:3
    for j = 1:ncols
        if i+j > 2
            ax(i,j) = uiaxes(g);
            ax(i,j).Layout.Row = i;
            ax(i,j).Layout.Column = j;
            colormap(ax(i,j), gray);
        end
    end
end

for i = 1:2
    for j = 1:ncols
        if i+j > 2
            ax(i,j).XTick = [];
            ax(i,j).YTick = [];
        end
    end
end

% Right buttons X pos and height/width
RBX = 1120;
RBH = 30; 
RBW = 55;

LOADMODE = 'MASK';
mvBtn = uibutton(fig, ...
    'Text', 'Load Mode: Masked', ...
    'Position', [15, 835, 120, 30]);
mvBtn.ButtonPushedFcn = @(btn, evt) toggleLOADMODE();

uibutton(fig, 'Text','Load 4D data', ...
    'Position',[15 865 120 30], ...
    'ButtonPushedFcn', @(btn,event) loadData());

% --- SAFE CLOSE BUTTON ---
closeBtn = uibutton(fig, ...
    'Text', 'Close GUI', ...
    'Position', [15, 805, 120, 30], ...
    'ButtonPushedFcn', @(btn, evt) closeApp());

s_slice = uislider(g,'Limits',[1 10],'MajorTicks',[],'Orientation','vertical');
s_slice.Layout.Row = 2;
s_slice.Layout.Column = ncols+1;

% Link slider callbacks to update function
addlistener(s_slice,'ValueChanged',@(src,evt) updateFcn());

% Keyboard arrow control
fig.WindowKeyPressFcn = @keyControl;

% Initiate segmentation buttons and set y-positions
y1 = 730; 
y2 = y1 - 2*RBH;
y3 = y2 - 2*RBH;
szDropdown = uidropdown(fig, ...
    'Items', {'15', '25', '35', '45', '55', '65', '75', '85', '95', '105', '115', '125'}, ...
    'Value', '15', ...
    'Position', [RBX, y1, RBW, RBH], ...
    'Tooltip', 'Set Patch Width (voxels before interp.)');

thrDropdown = uidropdown(fig, ...
    'Items', {'10', '20', '30', '40', '50', '60', '70'}, ...
    'Value', '20', ...
    'Position', [RBX, y2, RBW, RBH], ...
    'Tooltip', 'Set Local Threshold (%)');

cscDropdown = uidropdown(fig, ...
    'Items', {'10', '20', '30', '40', '50', '60', '70', '80', '90', '100'}, ...
    'Value', '80', ...
    'Position', [RBX, y3, RBW, RBH], ...
    'Tooltip', 'Set Velocity Range (%)');

% Add labels for the dropdowns
szLabel = uilabel(fig, ...
    'Position', [RBX, y1+RBH, 55, 20], ...
    'Text', 'SegWidth');

thrLabel = uilabel(fig, ...
    'Position', [RBX, y2+RBH, 55, 20], ...
    'Text', 'Threshold');

cscLabel = uilabel(fig, ...
    'Position', [RBX, y3+RBH, 55, 20], ...
    'Text', 'VelScale');

thrDropdown.ValueChangedFcn = @(src, evt) updateFlowPlane();
szDropdown.ValueChangedFcn = @(src, evt) onPatchOrThreshChange();
thrDropdown.ValueChangedFcn = @(src, evt) onPatchOrThreshChange();
cscDropdown.ValueChangedFcn = @(src, evt) togglePlay();

% Create dropdown and position it over ax(1,4)
toggle = uidropdown(fig, ...
    'Items', {'proj', 'vx', 'vy', 'vz'}, ...
    'Value', 'proj', ...
    'Position', [830, 865, RBW, RBH]);

% Play/Loop button
playBtn = uibutton(fig, ...
    'Text', 'Loop', ...
    'Position', [830, 835, 55, 30]);
playBtn.ButtonPushedFcn = @(btn, evt) togglePlay();

savedBtn = uibutton(fig, ... % show saved
    'Text', 'Show Saved', ...
    'Position', [RBX, 270, 75, 30]);
savedBtn.ButtonPushedFcn = @(btn, evt) toggleSaved();

printDelayBtn = uibutton(fig, ... % print delay
    'Text', 'Print Delay', ...
    'Position', [RBX, 240, 75, 30]);
printDelayBtn.ButtonPushedFcn = @(btn, evt) printDelay();

direction = []; % PCA and manual
direction.pca = [];
direction.man = [0, 0, 1]; %
rotation = [0, 0, 1]; % [x, y, z] 
DIRMODE = 'pca'; 
dirBtn = uibutton(fig, ...
    'Text', 'DIRECTION: PCA', ...
    'Position', [140, 835, 120, 30]); % xpos swap 900 -> 120
dirBtn.ButtonPushedFcn = @(btn, evt) toggleDIRMODE();

SHAPEMODE = 'CIRC';
shapeBtn = uibutton(fig, ...
    'Text', 'SHAPE: Circular', ...
    'Position', [140, 865, 120, 30]); % xpos swap 900 -> 120
shapeBtn.ButtonPushedFcn = @(btn, evt) toggleSHAPEMODE();

WFSHOWMODE = 'ALL';
wfsBtn = uibutton(fig, ...
    'Text', 'WF vis: Flow + PCA', ...
    'Position', [140, 805, 120, 30]);
wfsBtn.ButtonPushedFcn = @(btn, evt) toggleWFSHOWMODE();

% X rotation
yrpos = 845;
hrotbox = 25;
wrotbox = 55;
xRotBox = uieditfield(fig,'numeric', ...
    'Limits', [-180 180], ...
    'Value', rotation(1), ...
    'Position', [RBX, yrpos, wrotbox, hrotbox], ...
    'ValueChangedFcn', @(src, evt) setRotation(1, src.Value)); %#ok<*NASGU>

% Y rotation
yRotBox = uieditfield(fig,'numeric', ...
    'Limits', [-180 180], ...
    'Value', rotation(2), ...
    'Position', [RBX, yrpos-hrotbox, wrotbox, hrotbox], ...
    'ValueChangedFcn', @(src, evt) setRotation(2, src.Value));

% Z rotation
zRotBox = uieditfield(fig,'numeric', ...
    'Limits', [-180 180], ...
    'Value', rotation(3), ...
    'Position', [RBX, yrpos-2*hrotbox, wrotbox, hrotbox], ...
    'ValueChangedFcn', @(src, evt) setRotation(3, src.Value));

rotLabel = uilabel(fig, ...
    'Position', [RBX, yrpos+hrotbox, wrotbox, hrotbox], ...
    'Text', 'Set X/Y/Z');

% Create button panel for saving waveforms (row 3, col 4)
btnPanel = uipanel(g);
btnPanel.Layout.Row = 3;
btnPanel.Layout.Column = 4;

btnGrid = uigridlayout(btnPanel, [8,1]);
btnGrid.RowHeight = repmat({'1x'}, 1, 8);
btnGrid.ColumnWidth = {'1x'};

% TEMP edit: skip using LV for this anyway; add 2x CA + pCSF instead  
% btnNames = {'Left LV', 'Right LV', 'Left FMo', 'Right FMo', '3rd Vent', 'C. Aqueduct', '4th Vent', 'Spinal Canal'};
% fileNames = {'LLV.mat', 'RLV.mat', 'LFMo.mat', 'RFMo.mat', 'V3.mat', 'CA.mat', 'V4.mat', 'SC.mat'};

btnNames = {'Left FMo', 'Right FMo', '3rd Vent', 'L. Aqueduct', 'U. Aqueduct', '4th Vent', 'Spinal Canal', 'Perivasc.'};
fileNames = {'LFMo.mat', 'RFMo.mat', 'V3.mat', 'CAL.mat', 'CAU.mat', 'V4.mat', 'SC.mat', 'pCSF.mat'};

btnHandles = gobjects(length(btnNames), 1);
for k = 1:length(btnNames)
    btnHandles(k) = uibutton(btnGrid, 'Text', btnNames{k}, ...
        'ButtonPushedFcn', @(btn,event) saveWaveforms(fileNames{k}));
end

function onPatchOrThreshChange()
    updateDisplays();
    updateFlowPlane();
    % Update waveforms only if clicked point exists
    if ~isempty(clickedX) && ~isempty(clickedY)
        updateWaveformsFromCoords(clickedX, clickedY, round(s_slice.Value));
    end
end

function keyControl(~, event)
    step = 1;
    switch event.Key
        case 'rightarrow'  % Scroll forward through slices
            if s_slice.Value < s_slice.Limits(2)
                s_slice.Value = s_slice.Value + step;
                updateDisplays();
            end
        case 'leftarrow'  % Scroll backward through slices
            if s_slice.Value > s_slice.Limits(1)
                s_slice.Value = s_slice.Value - step;
                updateDisplays();
            end
    end
end

% Define cleanup function
function closeApp(~, ~)
    try
        % Stop and delete timer if active
        if ~isempty(playTimer) && isvalid(playTimer)
            stop(playTimer);
            delete(playTimer);
        end

        % Remove listeners safely
        delete(findall(0, 'Type', 'listener'));

        % Delete figure
        if isvalid(fig)
            delete(fig);
        end

        % Force garbage collection
        drawnow;
        pause(0.05);
        clearvars -except LOCAL group users user wrap2 BASEPATH
        fprintf('GUI closed and memory cleared.\n');
    catch ME
        warning('Error during GUI close: %s', ME.message);
    end
end

addlistener(s_slice,'ValueChanged',@(src,evt) updateDisplays());

data = struct();

isPlaying = false;
playTimer = [];

marker = [];
flowArrow = [];

clickedX = [];
clickedY = [];

% rms_clim = [];
mix_clim = [];

subjectFolder = ''; % To store folder path from loadData
savefolder = '';
WF = [];

clickableAxes = [];

% *** LOAD DATA ***
    function loadData()
        basefolder = uigetdir();
        if basefolder == 0, return; end

        subjectFolder = basefolder;  % Save basefolder for saving waveforms

        % Extract folders for RT2 and D.mat
        [foldername, subjname] = fileparts(subjectFolder);
        baseparent = fileparts(foldername);
        fig.Name = ['4D CSF Flow Viewer | ' subjname];

        % TODO load from SERVER - can we save D to server first also? 
        subjectFolder = fullfile(BASEPATH, subjname, 'processed');
        cubefile = fullfile(subjectFolder, 'cr2mag', 'r4dT2.nii'); 
        dmatfile = fullfile(subjectFolder, 'matproc', 'D.mat'); % NOT USED ATM
        if exist(dmatfile, 'file')
            load(dmatfile, 'D');
        end
        distMatFile = fullfile(baseparent, 'CSFmasks', subjname, 'distMat.mat'); % NOT USED ATM
        if exist(distMatFile, 'file')
            load(distMatFile, 'distMat'); % distance over branchMat branches
            data.distMat = distMat;
        end
        magfile = fullfile(subjectFolder, 'MAG.nii');

        savefolder = fullfile(subjectFolder, OUTFOLDER);
        if ~exist(savefolder, 'dir')
            mkdir(savefolder);
        end

        if LOCAL
            fvelsfolder = fullfile(LOCALVELS, subjname, 'fullvels');
            rvelsfolder = fullfile(LOCALVELS, subjname, 'brainvels'); % PRE MASKED TO REDUCE LOADING SIZE  
            disp('Loading from LOCAL brainvels folder')
        else
            fvelsfolder = fullfile(subjectFolder, 'fullvels');
            rvelsfolder = fullfile(subjectFolder, 'brainvels'); % PRE MASKED TO REDUCE LOADING SIZE 
            disp('Loading from REMOTE brainvels folder')
        end

        data.mag = MRIread(magfile).vol;
        data.cube = MRIread(cubefile).vol;

        if strcmp(LOADMODE, 'MASK')
            load([rvelsfolder, '/rx.mat'], 'rx');
            load([rvelsfolder, '/ry.mat'], 'ry');
            load([rvelsfolder, '/rz.mat'], 'rz');
            load([rvelsfolder, '/roi.mat'], 'roi');
            data.vx = rx;
            data.vy = ry;
            data.vz = rz;
        elseif strcmp(LOADMODE, 'FULL') % TEMP ZAYNAB: THIS LOAD IF NOT PREMASKED VELS (if .nii, might need to rotate, etc.)
            load([fvelsfolder, '/vx.mat'], 'vx');
            load([fvelsfolder, '/vy.mat'], 'vy');
            load([fvelsfolder, '/vz.mat'], 'vz');
            [nx, ny, nz, nt] = size(vx); % CORRECT DIMS?  
            data.vx = reshape(vx, nx*ny*nz, nt);
            data.vy = reshape(vy, nx*ny*nz, nt);
            data.vz = reshape(vz, nx*ny*nz, nt);
            data.roi = ones(size(data.mag));
        end

        % TEMP: adding FS seg instead of distance for CA definition? 
        aafile = fullfile(subjectFolder, 'FSproc', 'aa_nn4d.nii.gz');
        try
            data.aa = MRIread(aafile).vol;
            data.aa = imrotate(data.aa, -90);
        catch
            data.aa = zeros(size(data.mag));
        end

        data.mag = imrotate(data.mag, -90);
        data.cube = imrotate(data.cube, -90);
        data.dist = flip(D, 2); % NOT USED ATM
        data.roi = roi;

        % Global ROI and within-ROI velocities
        % data.groi = data.dist > 0; % uu, data.dist is modified from roi
        data.groi = data.roi > 0;
        data.ginds = find(data.groi(:));
        data.imap = zeros(size(data.groi));
        data.imap(data.ginds) = 1:numel(data.ginds);
        data.imap = flip(data.imap, 2);

        % TODO: filter before rms; also maybe do vstd, reduce BG errors 
        data.rms = zeros(size(data.imap));
        data.rms(data.ginds) = sqrt(mean(double(data.vx).^2 + double(data.vy).^2 + double(data.vz).^2, 2));
        data.rms = flip(data.rms, 2);
        data.mixed = data.rms.*data.cube;

        % TEMP: show vents instead of CSF dist when segmenting CA 
        VENTS = ismember(data.aa, [14 15 4 5 44 45]);
        data.dist = VENTS; 

        % Set display limits
        rms_clim = [0, 0.25 * max(data.rms(:))];
        % mix_clim = [0, 0.25 * max(data.mixed(:))];
        mix_clim = [0, 0.15 * max(data.mixed(:))];

        [~, ~, zres] = size(data.rms);
        s_slice.Limits = [1 zres];
        s_slice.Value = round(zres/2);

        currentFrame = 1; % Start at first frame

        % TEMP: taking from updateDisplays to here
        % clickableAxes = [ax(1,1), ax(1,2), ax(1,3), ax(2,1), ax(2,2), ax(2,3), ax(1,3)];
        clickableAxes = [ax(1,2), ax(1,3), ax(2,1), ax(2,2), ax(2,3), ax(1,3)];
        linkaxes([ax(2,1), ax(2,2), ax(2,3)], 'xy');

        updateDisplays();
        updateFlowPlane();

        disp('LoadData complete.')

    end

    function updateDisplays()
        if isempty(fieldnames(data)), return; end
        slice = round(s_slice.Value);

        if ~isempty(marker)
            % Only delete valid handles that still exist
            validMarkers = marker(ishandle(marker));
            if ~isempty(validMarkers)
                delete(validMarkers);
            end
            marker = [];
        end

        if isfield(data, 'distMat') && ~isempty(data.distMat)
            imagesc(ax(1,1), squeeze(data.distMat(:,:,slice))');
            axis(ax(1,1), 'image');
            title(ax(1,1), 'Centerline');
        end

        % TEMP: move to load function 
        % linkaxes([ax(2,1), ax(2,2), ax(2,3)], 'xy');

        % TEMP: no zoom state preserve
        prevLim = axis(ax(2,1));  % Save zoom state
        cla(ax(2,1));
        imagesc(ax(2,1), squeeze(data.dist(:,:,slice))');
        if ~isequal(prevLim, [0 1 0 1])  % If zoomed, restore view
            axis(ax(2,1), prevLim);
        else
            axis(ax(2,1), 'image');
        end
        title(ax(2,1), 'FS ROIs'); 

        prevLim = axis(ax(2,2));  % Save zoom state
        cla(ax(2,2));
        imagesc(ax(2,2), squeeze(data.cube(:,:,slice))');
        if ~isequal(prevLim, [0 1 0 1])  % If zoomed, restore view
            axis(ax(2,2), prevLim);
        else
            axis(ax(2,2), 'image');
        end
        title(ax(2,2), 'T2 CUBE'); 

        prevLim = axis(ax(2,3));  % Save zoom state
        cla(ax(2,3));
        imagesc(ax(2,3), squeeze(data.mixed(:,:,slice))'); % TEMP 
        if ~isequal(prevLim, [0 1 0 1])  % If zoomed, restore view
            axis(ax(2,3), prevLim);
        else
            axis(ax(2,3), 'image');
        end

        % Allow pointer to follow slice scroll
        if ~isempty(clickedX) && ~isempty(clickedY)
            updateWaveformsFromCoords(clickedX, clickedY, slice);
        end
        title(ax(2,3), 'CUBE x RMS velocity'); 
        clim(ax(2,3), mix_clim); 

        % TEMP: move to load function 
        % clickableAxes = [ax(1,1), ax(1,2), ax(1,3), ax(2,1), ax(2,2), ax(2,3), ax(1,3)];
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
            frame = currentFrame;
            imagesc(ax(1,4), data.proj2D(:,:,frame));
            % axis(ax(1,4), 'image');
            colormap(ax(1,4), gray);
            d2d = data.proj2D .* data.bseg;
            CSCALE = 0.01*str2double(cscDropdown.Value);
            clim(ax(1,4), [CSCALE * min(d2d(:)), CSCALE * max(d2d(:))]);
            % title(ax(1,4), sprintf('Flow plane (%s), frame %d', toggle.Value, frame));
            hold(ax(1,4), 'on');
            visboundaries(ax(1,4), data.bseg, 'Color', 'r', 'LineWidth', 1.5);
            hold(ax(1,4), 'off');
        end

    end

    function updateFlowPlane()
        if ~isfield(data, 'proj2D') || isempty(data.proj2D), return; end
        frame = currentFrame; % use current frame tracker

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
        d2d = data.proj2D .* data.bseg;
        CSCALE = 0.01*str2double(cscDropdown.Value);
        clim(ax(1,4), [CSCALE * min(d2d(:)), CSCALE * max(d2d(:))]);
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

        [xres, yres, zres, ~] = size(data.rms);
        if any([x y z] < 1) || x > xres || y > yres || z > zres
            uialert(fig, 'Invalid position for waveforms.', 'Save Error');
            return;
        end

        WF.vel.x = data.velx2D;
        WF.vel.y = data.vely2D;
        WF.vel.z = data.velz2D;
        WF.proj = data.proj2D;
        WF.flow = data.flow;
        WF.pc1 = data.pc1;
        WF.PC1 = data.PC1;
        if isfield(data, 'PC1')
            WF.PC1 = data.PC1;
        end

        WF.coords.x = x;
        WF.coords.y = y;
        WF.coords.z = z;
        WF.dist = data.dist(x, y, z);
        WF.patch = data.patch;
        WF.bseg = data.bseg; 
        WF.size = str2num(szDropdown.Value);
        WF.thr = str2num(thrDropdown.Value);

        save(fullfile(savefolder, filename), 'WF');
        disp(['Saved waveforms to ', fullfile(savefolder, filename)]);
    end

    function updateWaveformsFromClick(event, slice)
        pos = round(event.IntersectionPoint(1:2));
        clickedX = pos(1);
        clickedY = pos(2);
        x = clickedX;
        y = clickedY;
        updateWaveformsFromCoords(x, y, slice);
    end

    % This needs to be called for within ROI data
    function updateWaveformsFromCoords(x, y, z)

        ind = data.imap(x, y, z); % ind within ROI

        if ind == 0
            disp('Point outside ROI')
            return;
        end

        vx_t = squeeze(data.vx(ind, :))';
        vy_t = squeeze(data.vy(ind, :))';
        vz_t = squeeze(data.vz(ind, :))';

        % Volume PCA (shouldnt help alot when regularization is high)
        dilmap = zeros(size(data.imap));
        dilmap(x, y, z) = 1;
        dilmap = imdilate(dilmap, strel('sphere', 2));
        inds = data.imap(find(dilmap));
        inds(inds==0) = [];

        VX_t = squeeze(data.vx(inds, :))';
        VY_t = squeeze(data.vy(inds, :))';
        VZ_t = squeeze(data.vz(inds, :))';

        plotWaveforms(ax(3,1), vx_t, vy_t, vz_t);

        % Perform PCA to get CSF flow direction
        V = [vx_t vy_t vz_t]; % (frames x 3)
        [coeff, score, ~] = pca(V, 'NumComponents', 1);
        pc1 = score(:,1); % Time series
        direction.pca = coeff(:,1); % PCA direction loadings

        % Perform PCA to get CSF flow direction
        V = [VX_t VY_t VZ_t]; % (frames x 3)
        [~, score, ~] = pca(V, 'NumComponents', 1);
        PC1 = score(:,1); % Time series

        patch_width = str2double(szDropdown.Value);
        local_thresh = str2double(thrDropdown.Value);

        % *** SPECIFY IN TOP OF GUI *** 
        if strcmp(SHAPEMODE, 'RECT')
            SEGMODE = 'cube'; % alternative, rms or cube
            [flow, cube_patch, patch, bseg, ~, proj2D] = ... 
                extractThroughPlaneFlow_V3D(data, [x, y, z], direction, patch_width, SEGMODE, local_thresh, DIRMODE, SHAPEMODE);
        elseif strcmp(SHAPEMODE, 'CIRC')
            SEGMODE = 'mixed'; % alternative, rms or cube
            [flow, cube_patch, patch, bseg, ~, proj2D] = ... 
                extractThroughPlaneFlow_V3D(data, [x, y, z], direction, patch_width, SEGMODE, local_thresh, DIRMODE, SHAPEMODE);
        elseif strcmp(SHAPEMODE, 'ANY') % TEMP: use mixed for any 
            SEGMODE = 'mixed'; % alternative, rms or cube
            [flow, cube_patch, patch, bseg, ~, proj2D] = ... 
                extractThroughPlaneFlow_V3D(data, [x, y, z], direction, patch_width, SEGMODE, local_thresh, DIRMODE, SHAPEMODE);
        end

        data.proj2D = proj2D.proj;
        data.velx2D = proj2D.velx;
        data.vely2D = proj2D.vely;
        data.velz2D = proj2D.velz;
        data.patch = patch; % this is actually patch_interp
        data.bseg = bseg;
        data.flow = flow;
        data.pc1 = pc1;
        data.PC1 = PC1; 

        % Clear old arrow if it exists
        if isgraphics(flowArrow)
            delete(flowArrow);
            flowArrow = [];
        end
        dir_xy = -direction.(DIRMODE)([2, 3]);  % vy (x), vz (z)
        arrowLength = 15;
        if norm(dir_xy) > 0
            dir_xy = dir_xy / norm(direction.(DIRMODE)) * arrowLength;
        else
            dir_xy = [0 arrowLength];
        end
        if ~isempty(marker)
            delete(marker(ishandle(marker)));
        end
        marker = [];
        for j = 1:3, hold(ax(2,j), 'on'); end
        marker = [ ...
            quiver(ax(2,1), x, y, dir_xy(1), dir_xy(2), 'Color', 'r', 'LineWidth', 2.5, 'AutoScale', 'off', 'MaxHeadSize', 2.5);
            quiver(ax(2,2), x, y, dir_xy(1), dir_xy(2), 'Color', 'r', 'LineWidth', 2.5, 'AutoScale', 'off', 'MaxHeadSize', 2.5);
            quiver(ax(2,3), x, y, dir_xy(1), dir_xy(2), 'Color', 'r', 'LineWidth', 2.5, 'AutoScale', 'off', 'MaxHeadSize', 2.5)];
        for j = 1:3, hold(ax(2,j), 'off'); end

        % Local CS
        imagesc(ax(1,2), cube_patch);
        axis(ax(1,2), 'image');
        title(ax(1,2), 'Flow plane: T2 CUBE');
        colormap(ax(1,2), gray);
        hold(ax(1,2), 'on');
        visboundaries(ax(1,2), data.bseg, 'Color', 'r', 'LineWidth', 1.5);
        hold(ax(1,2), 'off');
        ax(1,2).XTick = [];
        ax(1,2).YTick = [];

        % Local CS
        imagesc(ax(1,3), patch);
        axis(ax(1,3), 'image');
        title(ax(1,3), ['Flow plane: ' SEGMODE]);
        colormap(ax(1,3), gray);
        hold(ax(1,3), 'on');
        visboundaries(ax(1,3), data.bseg, 'Color', 'r', 'LineWidth', 1.5);
        hold(ax(1,3), 'off');
        ax(1,3).XTick = [];
        ax(1,3).YTick = [];

        % Clear and prepare left axis
        cla(ax(3,2), 'reset');

        if strcmp(WFSHOWMODE, 'ALL')
            yyaxis(ax(3,2), 'left');
            hold(ax(3,2), 'on');
    
            plot(ax(3,2), zscore(pc1), 'c--', 'LineWidth', 1.5); %#ok<*UNRCH> % voxel PCA 
            plot(ax(3,2), zscore(PC1), 'b-', 'LineWidth', 1.5); %#ok<*UNRCH> % volume PCA 
            ylabel(ax(3,2), 'PCA vel. (z-score)');
    
            % Right Y-axis: Flow through plane
            yyaxis(ax(3,2), 'right');
            hold(ax(3,2), 'on'); % need to reapply hold on for each axis?
            plot(ax(3,2), flow, 'Color', 'r', 'LineWidth', 1.5);
            ylabel(ax(3,2), 'Flow rate (ml/s)');
    
            % Final touches
            title(ax(3,2), sprintf('PC-1 dir: [%.2f  %.2f  %.2f]', direction.pca(1), direction.pca(2), direction.pca(3)));
            % xlabel(ax(3,2), 'Frame');
            legend(ax(3,2), 'Voxel PCA', 'Volume PCA', 'Through-plane flow')
        elseif strcmp(WFSHOWMODE, 'FLOW') % Just show the through-plane flow rate waveforms
            hold(ax(3,2), 'on'); % need to reapply hold on for each axis?
            plot(ax(3,2), flow, 'Color', 'r', 'LineWidth', 1.5);
            ylabel(ax(3,2), 'Flow rate (ml/s)');
            title(ax(3,2), sprintf('PC-1 dir: [%.2f  %.2f  %.2f]', direction.pca(1), direction.pca(2), direction.pca(3)));
            legend(ax(3,2), 'Through-plane flow')
        end

        % Plot coronal RMS view in ax(1,3)
        cla(ax(2,4), 'reset');
        cor = squeeze(data.cube(clickedX, :, :)); % could use data.mixed
        nz = size(cor, 2);
        ys = clickedY-nz/2; ye = clickedY+nz/2;
        if ye > size(data.cube, 2)
            ye = size(data.cube, 2);
        end
        cor = cor(ys:ye, :);
        imagesc(ax(2,4), cor);
        axis(ax(2,4), 'image');
        colormap(ax(2,4), gray);
        title(ax(2,4), 'CUBE \times RMS (coronal)');
        clim(ax(2,4), [0, 1.0 * max(data.cube(:))]);
        hold(ax(2,4), 'on');

        % Add velocity direction as quiver
        dir_xz = direction.(DIRMODE)([1, 3]);  % vx and vz components
        dir_xz = dir_xz / norm(direction.pca) * 10;  % scale for visibility

        quiver(ax(2,4), z, (nz/2) - 1, dir_xz(1), -dir_xz(2), ...
            'Color', 'r', 'LineWidth', 1.5, 'AutoScale', 'off', 'MaxHeadSize', 1.5);

        hold(ax(1,3), 'off');
        ax(2,4).XTick = [];
        ax(2,4).YTick = [];

    end

    % TEMP voxel PCA or flow rate? 
    function printDelay()
        d = dir(fullfile(savefolder, '*.mat'));
        nwf = numel(d);
        sc = [];
        try 
            load(fullfile(savefolder, 'SC.mat'), 'WF');
            % sc = WF.pc1;
            sc = WF.flow';
            sc = zscore(sc);
        catch
            disp('No point in SC');
        end
        for i = 1:nwf
            wn = fullfile(savefolder, d(i).name);
            load(wn, 'WF');
            if contains(d(i).name, 'LV')
                wf = WF.PC1;
            else 
                % wf = WF.pc1;
                wf = WF.flow';
            end
            if ~isempty(sc)
                if corr(sc, wf) < 0
                    wf = -wf;
                end
            end
            
            % wf = smoothdata(wf);
            % wf = zscore(wf);

            maxlag = 5; % 5 frames ~ 250 ms 
            [rmax, mlag] = waveformCoupling_V3D(sc', wf', maxlag);
            disp(['ROI: ' d(i).name ' | rmax: ' num2str(rmax) ' | mlag: ' num2str(mlag)])
        end
    end

    % TODO this clears correctly? 
    function toggleSaved()
        d = dir(fullfile(savefolder, '*.mat'));
        cla(ax(3,3), 'reset');
        hold(ax(3,3), 'off');

        y_coords = zeros(numel(d), 1);
        z_coords = zeros(numel(d), 1);
        plotHandles = gobjects(numel(d), 1); % Preallocate handles

        % Try to load SC reference waveform
        sc = [];
        try
            load(fullfile(savefolder, 'SC.mat'), 'WF');
            sc = WF.pc1;
        catch
            disp('No point in SC');
        end

        for i = 1:numel(d)
            wn = fullfile(savefolder, d(i).name);
            load(wn, 'WF');

            if contains(d(i).name, 'LV')
                wf = WF.PC1;
            else
                % wf = WF.pc1;
                wf = WF.flow';
            end

            % Flip sign if negatively correlated with SC
            if ~isempty(sc)
                if corr(sc, wf) < 0
                    wf = -wf;
                end
            end

            % Plot and store handle
            hold(ax(3,3), 'on');

            % wf = smoothdata(wf);
            x = 1:numel(wf);
            xq = linspace(1, 20, 100);
            wf = interp1(x,wf,xq, "pchip");
            wf = zscore(wf);
            plotHandles(i) = plot(ax(3,3), wf, 'LineWidth', 1.5);

            % Save coords
            y_coords(i) = WF.coords.x;
            z_coords(i) = WF.coords.y;
        end

        % Create legend with handles + filenames
        legend(ax(3,3), plotHandles, {d.name}, 'Interpreter', 'none');

        % Clear and re-scatter points
        delete(scatterHandles);
        updateDisplays();
        for j = 1:3
            hold(ax(2,j), 'on');
            scatterHandles = scatter(ax(2,j), y_coords, z_coords, 10, ...
                'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r');
        end
    end

    function togglePlay()
        % Create timer if it doesn't exist
        if isempty(playTimer) || ~isvalid(playTimer)
            playTimer = timer( ...
                'ExecutionMode', 'fixedRate', ...
                'Period', 0.2, ...
                'TimerFcn', @(~,~) nextFrame());
        end

        % Toggle playback state
        if ~isPlaying
            isPlaying = true;
            playBtn.Text = 'Stop';
            start(playTimer);
        else
            isPlaying = false;
            playBtn.Text = 'Play';
            stop(playTimer);
        end
    end

    function nextFrame()
        if isempty(data) || ~isfield(data, 'proj2D') || isempty(data.proj2D), return; end
        if isempty(maxFrame) || maxFrame < 1, return; end

        currentFrame = currentFrame + 1;
        if currentFrame > maxFrame
            currentFrame = 1;
        end
        updateFlowPlane();
    end

    function toggleLOADMODE()
        if strcmp(LOADMODE, 'MASK')
            LOADMODE = 'FULL';
            mvBtn.Text = 'Load Mode: Full';
        elseif strcmp(LOADMODE, 'FULL')
            LOADMODE = 'MASK';
            mvBtn.Text = 'Load Mode: Masked';
        end
    end

    function toggleDIRMODE()
        if strcmp(DIRMODE, 'pca')
            DIRMODE = 'man';
            dirBtn.Text = 'DIRECTION: Manual';
        else
            DIRMODE = 'pca';
            dirBtn.Text = 'DIRECTION: PCA';
        end

        % Update flow-plane and waveforms if a point is selected
        if ~isempty(clickedX) && ~isempty(clickedY)
            updateWaveformsFromCoords(clickedX, clickedY, round(s_slice.Value));
        end
    end

    function toggleSHAPEMODE()
        if strcmp(SHAPEMODE, 'CIRC')
            SHAPEMODE = 'RECT';
            shapeBtn.Text = 'SHAPE: Rect.';
        elseif strcmp(SHAPEMODE, 'RECT')
            SHAPEMODE = 'ANY';
            shapeBtn.Text = 'SHAPE: Any';
        elseif strcmp(SHAPEMODE, 'ANY')
            SHAPEMODE = 'CIRC';
            shapeBtn.Text = 'SHAPE: Circular';
        end

        % Update flow-plane and waveforms if a point is selected
        if ~isempty(clickedX) && ~isempty(clickedY)
            updateWaveformsFromCoords(clickedX, clickedY, round(s_slice.Value));
        end
    end

    function toggleWFSHOWMODE()
        if strcmp(WFSHOWMODE, 'ALL')
            WFSHOWMODE = 'FLOW';
            wfsBtn.Text = 'WF vis: Flow only';
        elseif strcmp(WFSHOWMODE, 'FLOW')
            WFSHOWMODE = 'ALL';
            wfsBtn.Text = 'WF vis: Flow + PCA';
        end

        % Update flow-plane and waveforms if a point is selected
        if ~isempty(clickedX) && ~isempty(clickedY)
            updateWaveformsFromCoords(clickedX, clickedY, round(s_slice.Value));
        end
    end

    function setRotation(idx, val)
        rotation(idx) = val;  % update rotation array

        % Normalize vector to length 1
        direction.man = rotation / norm(rotation);
        direction.man = direction.man(:);
        % disp(['Manual direction set: ' num2str(direction.man')])

        % If a point is selected, update flow-plane
        if ~isempty(clickedX) && ~isempty(clickedY)
            updateWaveformsFromCoords(clickedX, clickedY, round(s_slice.Value));
        end
    end

end
