% 4D CSF flow analysis tool 

% TODO: export all frames from current flow-proj with button press 

function CSF_GUI

% ----------------------------SETTINGS---------------------------------- %

USER = 'TV';
TV_OUT = 'TEMP';

% USER = 'ZY';
ZY_OUT = 'Z:\CVMRIGroup\Users\zsy001\anti-FLAIR\quant\AF2026';

% USER = '';

if strcmp(USER, 'ZY')
    LOCAL = false;
    BASEPATH = 'Z:/CVMRIGroup/Users/txv016/WRAP2/niis/CURRENT/'; % Where all subj-folders are 
    LOCALVELS = ''; 
    OUTFOLDER = ZY_OUT;
    addpath('Z:\CVMRIGroup\Users\zsy001\anti-FLAIR\quant\Freesurfer Tools')
elseif strcmp(USER, 'TV')
    LOCAL = true; % Note: LOCAL might fail with load mode full if not rsynced
    BASEPATH = '/Volumes/groups/CVMRIGroup/Users/txv016/WRAP2/niis/CURRENT/'; 
    LOCALVELS = '/Users/txv016/Documents/BRAINVELS'; % This is only used if LOCAL is true 
    OUTFOLDER = TV_OUT;
    addpath('/Users/txv016/Documents/MATLAB/Freesurfer Tools') % TODO: better use niftiread
else
    disp('*** Add USER and specify settings for loading/saving data ***')
end

AX21 = 'MAG'; % If not specified, 4D-MAG vill be placed in ax(2,1)
% or 'T1', 'FST1', 'dist', 'VENTS'

% ---------------------------------------------------------------------- %

currentFrame = 1; % Track current frame
maxFrame = 20;
scatterHandles = [];

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
SEG_RBH = RBH - 5;      % seg dropdown height
SEG_DY  = 2*RBH - 15;    % vertical step between seg dropdowns (was 2*RBH)

LOADMODE = 'MASKED';
mvBtn = uibutton(fig, ...
    'Text', 'Load Mode: Masked', ...
    'Position', [15, 835, 120, 30]);
mvBtn.ButtonPushedFcn = @(btn, evt) toggleLOADMODE();

uibutton(fig, 'Text','Load 4D data', ...
    'Position',[15 865 120 30], ...
    'ButtonPushedFcn', @(btn,event) loadData());

closeBtn = uibutton(fig, ...
    'Text', 'Close GUI', ...
    'Position', [15, 745, 120, 30], ...
    'ButtonPushedFcn', @(btn, evt) closeApp());

s_slice = uislider(g,'Limits',[1 10],'MajorTicks',[],'Orientation','vertical');
s_slice.Layout.Row = 2;
s_slice.Layout.Column = ncols+1;

% Link slider callbacks to update function
addlistener(s_slice,'ValueChanged',@(src,evt) updateFcn());

% Keyboard arrow control
fig.WindowKeyPressFcn = @keyControl;

% Initiate segmentation buttons and set y-positions
y1 = 750; 
y2 = y1 - SEG_DY;
y3 = y2 - SEG_DY;
y4 = y3 - SEG_DY;
y5 = y4 - SEG_DY;
szDropdown = uidropdown(fig, ...
    'Items', {'15', '25', '35', '45', '55', '65', '75', '85', '95', '105', '115', '125'}, ...
    'Value', '25', ...
    'Position', [RBX, y1, RBW, SEG_RBH], ...
    'Tooltip', 'Set Patch Width (voxels before interp.)');

thrDropdown = uidropdown(fig, ...
    'Items', {'10', '20', '30', '40', '50', '60', '70'}, ...
    'Value', '30', ... % edit to 30 % standard threshold 
    'Position', [RBX, y2, RBW, SEG_RBH], ...
    'Tooltip', 'Set Local Threshold (%)');

cscDropdown = uidropdown(fig, ...
    'Items', {'10', '20', '30', '40', '50', '60', '70', '80', '90', '100'}, ...
    'Value', '80', ...
    'Position', [RBX, y3, RBW, SEG_RBH], ...
    'Tooltip', 'Set Velocity Range (%)');

clipDropdown = uidropdown(fig, ...
    'Items', {'Off', '90', '95', '99', '100'}, ...
    'Value', '99', ...
    'Position', [RBX, y4, RBW, SEG_RBH], ...
    'Tooltip', 'Clip seg patch at percentile before threshold (Off = no clip)');

dilateDropdown = uidropdown(fig, ...
    'Items', {'0', '1', '2', '3'}, ...
    'Value', '0', ...
    'Position', [RBX, y5, RBW, SEG_RBH], ...
    'Tooltip', 'Segmentation mask dilation radius (0 = off)');

% Add labels for the dropdowns
szLabel = uilabel(fig, ...
    'Position', [RBX, y1+SEG_RBH, 55, 20], ...
    'Text', 'SegWidth');

thrLabel = uilabel(fig, ...
    'Position', [RBX, y2+SEG_RBH, 55, 20], ...
    'Text', 'Threshold');

cscLabel = uilabel(fig, ...
    'Position', [RBX, y3+SEG_RBH, 55, 20], ...
    'Text', 'VelScale');

clipLabel = uilabel(fig, ...
    'Position', [RBX, y4+SEG_RBH, 55, 20], ...
    'Text', 'ClipVal');

dilateLabel = uilabel(fig, ...
    'Position', [RBX, y5+SEG_RBH, 55, 20], ...
    'Text', 'Dilate');

thrDropdown.ValueChangedFcn = @(src, evt) updateFlowPlane();
szDropdown.ValueChangedFcn = @(src, evt) onPatchOrThreshChange();
thrDropdown.ValueChangedFcn = @(src, evt) onPatchOrThreshChange();
cscDropdown.ValueChangedFcn = @(src, evt) togglePlay();
clipDropdown.ValueChangedFcn = @(src, evt) onPatchOrThreshChange();
dilateDropdown.ValueChangedFcn = @(src, evt) onPatchOrThreshChange();

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

exportBtn = uibutton(fig, ...
    'Text', 'Export', ...
    'Position', [830, 805, 55, 30]);
exportBtn.ButtonPushedFcn = @(btn, evt) exportProj();

manuRoiBtn = uibutton(fig, ...
    'Text', 'ROI', ...
    'Position', [830, 775, 55, 30]);
manuRoiBtn.ButtonPushedFcn = @(btn, evt) runManuROI();

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
    'Text', 'Direction: PCA', ...
    'Position', [15, 805, 120, 30]);
dirBtn.ButtonPushedFcn = @(btn, evt) toggleDIRMODE();

WFSHOWMODE = 'FLOW';
wfsBtn = uibutton(fig, ...
    'Text', 'WF-Vis: Flow only', ...
    'Position', [15, 775, 120, 30]);
wfsBtn.ButtonPushedFcn = @(btn, evt) toggleWFSHOWMODE();

SHAPEMODE = 'ANY'; % ANY has no effect 
shapeBtn = uibutton(fig, ...
    'Text', 'Shape: Any', ...
    'Position', [140, 865, 120, 30]);
shapeBtn.ButtonPushedFcn = @(btn, evt) toggleSHAPEMODE();

LOCALCSMODE = 'FULL';
localCsBtn = uibutton(fig, ...
    'Text', 'Local CS: Full', ...
    'Position', [140, 835, 120, 30]);
localCsBtn.ButtonPushedFcn = @(btn, evt) toggleLOCALCSMODE();

SEGMODE = 'VxVSTD';
segBtn = uibutton(fig, ...
    'Text', 'SegVol: T2xVSTD', ...
    'Position', [140, 805, 120, 30]);
segBtn.ButtonPushedFcn = @(btn, evt) toggleSEGMODE();

CUBEMODE = 'CUBE';
cubeBtn = uibutton(fig, ...
    'Text', 'T2-w: CUBE', ...
    'Position', [140, 775, 120, 30], ...
    'Enable', 'off');
cubeBtn.ButtonPushedFcn = @(btn, evt) toggleCUBEAF();

autoReg3DBtn = uibutton(fig, ...
    'Text', 'AutoReg-3D', ...
    'Position', [140, 745, 120, 30]);
autoReg3DBtn.ButtonPushedFcn = @(btn, evt) runAutoReg3D();

manuRegBtn = uibutton(fig, ...
    'Text', 'ManuReg-3D', ...
    'Position', [140, 715, 120, 30]);
manuRegBtn.ButtonPushedFcn = @(btn, evt) runManuReg3D();

restoreT2wBtn = uibutton(fig, ...
    'Text', 'Restore T2-w', ...
    'Position', [140, 685, 120, 30]);
restoreT2wBtn.ButtonPushedFcn = @(btn, evt) restoreT2w();

showVolBtn = uibutton(fig, ...
    'Text', 'Show Volumes', ...
    'Position', [140, 655, 120, 30]);
showVolBtn.ButtonPushedFcn = @(btn, evt) showOverlay2D();

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

% TEMP edit: skip LV buttons for now; add 2x CA + pCSF instead  
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

subjectFolder = ''; % To store folder path from loadData
savefolder = '';
WF = [];

clickableAxes = [];

subjname = [];

    function t = anatTitle()
        if strcmp(CUBEMODE, 'AF')
            t = 'Anti-FLAIR';
        else
            t = 'T2 CUBE';
        end
    end

    function t = anatTitleWithShift()
        t = anatTitle();
        if isempty(fieldnames(data)) || ~isfield(data, 'shift')
            return;
        end
        if strcmp(CUBEMODE, 'AF') && isfield(data, 'hasAF') && data.hasAF
            sh = data.shift.AF;
        else
            sh = data.shift.T2CUBE;
        end
        t = sprintf('%s %s', t, formatRegShift(sh));
    end

    function t = vstdTitle()
        if strcmp(CUBEMODE, 'AF')
            t = 'AFxVSTD';
        else
            t = 'T2xVSTD';
        end
    end

    function f = segVolField()
        if strcmp(SEGMODE, 'ANAT')
            f = 'CUBE';
        else
            f = 'VxVSTD';
        end
    end

    function segLabel = segVolLabel()
        if strcmp(SEGMODE, 'ANAT')
            segLabel = anatTitle();
        else
            segLabel = vstdTitle();
        end
    end

    function updateSegBtnText()
        segBtn.Text = ['SegVol: ' segVolLabel()];
    end

    function refreshVolMaps()
        if ~isfield(data, 'T2CUBE')
            return
        end
        data.T2xVSTD = data.vstd .* data.T2CUBE;
        if isfield(data, 'hasAF') && data.hasAF
            data.AFxVSTD = data.vstd .* data.AF;
        end
        if strcmp(CUBEMODE, 'AF') && isfield(data, 'hasAF') && data.hasAF
            data.CUBE = data.AF;
            data.VxVSTD = data.AFxVSTD;
        else
            CUBEMODE = 'CUBE';
            data.CUBE = data.T2CUBE;
            data.VxVSTD = data.T2xVSTD;
        end
    end

    function updateCubeBtnState()
        if isfield(data, 'hasAF') && data.hasAF
            cubeBtn.Enable = 'on';
            if strcmp(CUBEMODE, 'AF')
                cubeBtn.Text = 'T2-w: Anti-FLAIR';
            else
                cubeBtn.Text = 'T2-w: CUBE';
            end
        else
            cubeBtn.Enable = 'off';
            CUBEMODE = 'CUBE';
            cubeBtn.Text = 'T2-w: CUBE';
        end
    end

    function tf = velocitiesLoaded()
        tf = isfield(data, 'hasVel') && data.hasVel;
    end

    function runAutoReg3D()
        if isempty(fieldnames(data))
            uialert(fig, 'Load data first.', 'AutoReg-3D');
            return;
        end

        slice = round(s_slice.Value);

        if isfield(data, 'hasAF') && data.hasAF
            data = coregSegmodeToMag3D(fig, data, ax(2,1), ax(2,2), slice, 'AF', 'T2CUBE');
            refreshVolMaps();
        end

        data = coregSegmodeToMag3D(fig, data, ax(2,1), ax(2,2), slice, 'T2CUBE', 'mag');
        refreshVolMaps();
        updateDisplays();
        if ~isempty(clickedX) && ~isempty(clickedY)
            updateWaveformsFromCoords(clickedX, clickedY, slice);
        end
    end

    function movingField = anatMovingField()
        if strcmp(CUBEMODE, 'AF') && isfield(data, 'hasAF') && data.hasAF
            movingField = 'AF';
        else
            movingField = 'T2CUBE';
        end
    end

    function runManuReg3D()
        if isempty(fieldnames(data))
            uialert(fig, 'Load data first.', 'ManuReg-3D');
            return;
        end

        movingField = anatMovingField();
        slice = round(s_slice.Value);
        center = [];
        if ~isempty(clickedX) && ~isempty(clickedY)
            center = [clickedX, clickedY, slice];
        end
        patch_width = str2double(szDropdown.Value);
        data = ManuReg3D(fig, data, ax(2,1), slice, movingField, segVolField(), ...
            center, direction, DIRMODE, patch_width);
        refreshVolMaps();
        updateCubeBtnState();
        updateSegBtnText();
        updateDisplays();
        if ~isempty(clickedX) && ~isempty(clickedY)
            updateWaveformsFromCoords(clickedX, clickedY, slice);
        end
    end

    function runManuROI()
        if isempty(fieldnames(data))
            uialert(fig, 'Load data first.', 'ROI');
            return;
        end
        if isempty(clickedX) || isempty(clickedY)
            uialert(fig, 'Click a point on the image first.', 'ROI');
            return;
        end

        x = clickedX;
        y = clickedY;
        z = round(s_slice.Value);

        patch_width = str2double(szDropdown.Value);
        local_thresh = str2double(thrDropdown.Value);
        if strcmpi(clipDropdown.Value, 'Off')
            clip_on = false;
            clip_val = 99;
        else
            clip_on = true;
            clip_val = str2double(clipDropdown.Value);
        end
        dilate_val = str2double(dilateDropdown.Value);

        [~, ~, patch_interp, ~, ~, ~] = extractThroughPlaneFlow_V3D( ...
            data, [x, y, z], direction, patch_width, segVolField(), local_thresh, ...
            DIRMODE, SHAPEMODE, clip_on, clip_val, dilate_val, [], LOCALCSMODE);

        bsegManual = ManuSegROI(patch_interp);
        if isempty(bsegManual) || ~any(bsegManual(:))
            return;
        end

        updateWaveformsFromCoords(x, y, z, bsegManual);
    end

    function restoreT2w()
        if isempty(fieldnames(data))
            uialert(fig, 'Load data first.', 'Restore T2-w');
            return;
        end
        if ~isfield(data, 'noreg') || ~isfield(data.noreg, 'T2')
            uialert(fig, 'No original T2-w volumes stored.', 'Restore T2-w');
            return;
        end

        data.T2CUBE = data.noreg.T2 + 0;
        if isfield(data, 'hasAF') && data.hasAF && isfield(data.noreg, 'AF')
            data.AF = data.noreg.AF + 0;
        end
        data = resetRegShifts(data);
        refreshVolMaps();
        updateCubeBtnState();
        updateSegBtnText();
        updateDisplays();
        if ~isempty(clickedX) && ~isempty(clickedY)
            updateWaveformsFromCoords(clickedX, clickedY, round(s_slice.Value));
        end
        disp('Restored T2CUBE and AF from data.noreg');
    end

% *** LOAD DATA ***
    function loadData()
        basefolder = uigetdir();
        if basefolder == 0, return; end

        subjectFolder = basefolder;  % Save basefolder for saving waveforms

        % Extract folders for RT2 and D.mat
        [foldername, subjname] = fileparts(subjectFolder);
        baseparent = fileparts(foldername);
        fig.Name = ['4D CSF Flow Viewer | ' subjname];
        disp(['Loading ' subjname]);

        % Load from server
        subjectFolder = fullfile(BASEPATH, subjname, 'processed');

        % Orig
        magfile = fullfile(subjectFolder, 'MAG.nii');
        cubefile = fullfile(subjectFolder, 'cr2mag', 'r4dT2.nii'); 
        antiflairfile = fullfile(subjectFolder, 'cr2mag', 'r4dAF.nii');
        dmatfile = fullfile(subjectFolder, 'matproc', 'D.mat'); 

        % TEMP WRAP CA analysis 
        savefolder = fullfile(subjectFolder, OUTFOLDER);
        if exist(savefolder, 'dir')
            disp(['Outfolder exists: ' savefolder])
        else
            disp(['Outfolder does not exist: ' savefolder])
            if exist(subjectFolder, 'dir')
                mkdir(savefolder);
                disp(['Outfolder created: ' savefolder])
            else
                disp('Subject folder missing')
                MAKEERROR;
            end
        end

        if exist(savefolder, 'dir')
            matFiles = dir(fullfile(savefolder, '*.mat'));
            if isempty(matFiles)
                disp('Outfolder has no .mat files')
            else
                disp('Outfolder .mat files:')
                for k = 1:numel(matFiles)
                    disp(['  ' matFiles(k).name])
                end
            end
        end

        if ~strcmp(LOADMODE, 'NONE')
            if LOCAL
                fvelsfolder = fullfile(LOCALVELS, subjname, 'fullvels');
                rvelsfolder = fullfile(LOCALVELS, subjname, 'brainvels'); % PRE MASKED TO REDUCE LOADING SIZE  
                disp('Loading from LOCAL brainvels folder')
            else
                fvelsfolder = fullfile(subjectFolder, 'fullvels');
                rvelsfolder = fullfile(subjectFolder, 'brainvels'); % PRE MASKED TO REDUCE LOADING SIZE 
                disp('Loading from REMOTE brainvels folder')
            end
        end

        disp('---Loading structural MRIs.---')
        data.mag = readNiiVol(magfile);
        disp('Loaded 4D-Mag (MAG.nii)');
        data.T2CUBE = readNiiVol(cubefile);
        disp('Loaded T2-CUBE (r4dT2.nii)');
        data.mag = imrotate(data.mag, -90);
        data.T2CUBE = imrotate(data.T2CUBE, -90);

        data.hasAF = false;
        if ~isempty(findNiiPath(antiflairfile))
            try
                data.AF = readNiiVol(antiflairfile);
                data.AF = imrotate(data.AF, -90);
                data.hasAF = true;
                disp('Loaded anti-FLAIR (r4dAF.nii)');
            catch ME
                disp(['anti-FLAIR load failed: ' ME.message]);
            end
        else
            disp('No anti-FLAIR found (r4dAF.nii)');
        end

        data.noreg.T2 = data.T2CUBE + 0;
        if data.hasAF
            data.noreg.AF = data.AF + 0;
        end

        data = initRegShifts(data);

        CUBEMODE = 'CUBE';
        SEGMODE = 'VxVSTD';

        disp('---Loading velocities.---')
        if strcmp(LOADMODE, 'MASKED')
            load([rvelsfolder, '/rx.mat'], 'rx');
            disp('rx velocity loaded')
            load([rvelsfolder, '/ry.mat'], 'ry');
            disp('ry velocity loaded')
            load([rvelsfolder, '/rz.mat'], 'rz');
            disp('rz velocity loaded')
            load([rvelsfolder, '/roi.mat'], 'roi');
            data.vx = rx;
            data.vy = ry;
            data.vz = rz;
            data.roi = roi;
            data.hasVel = true;
        elseif strcmp(LOADMODE, 'FULL') % Note these are still same format as masked velocities atm, just a mask of ones 
            load([fvelsfolder, '/rx.mat'], 'rx');
            disp('rx velocity loaded')
            load([fvelsfolder, '/ry.mat'], 'ry');
            disp('ry velocity loaded')
            load([fvelsfolder, '/rz.mat'], 'rz');
            disp('rz velocity loaded')
            data.vx = rx;
            data.vy = ry;
            data.vz = rz;
            data.roi = ones(size(data.mag));
            data.hasVel = true;
        else
            disp('Load mode NONE: skipping velocity load');
            data.hasVel = false;
        end

        % TEMP: adding T1-w and FS seg instead of distance for CA definition? 
        aafile = fullfile(subjectFolder, 'FSproc', 'aa_nn4d.nii.gz');
        t1file = fullfile(subjectFolder, 'FSproc', 'T1_nn4d.nii.gz');
        if strcmp(LOADMODE, 'NONE') % No need to load T1-w and FS output 
            aafile = ''; t1file = ''; 
        end
        try
            if ~isempty(t1file) && ~isempty(findNiiPath(t1file))
                data.T1 = readNiiVol(t1file);
                data.T1 = imrotate(data.T1, -90);
            else
                data.T1 = zeros(size(data.mag));
            end
            if ~isempty(aafile) && ~isempty(findNiiPath(aafile))
                data.aa = readNiiVol(aafile);
                data.aa = imrotate(data.aa, -90);
            else
                data.aa = zeros(size(data.mag));
            end
        catch ME
            disp(['T1/FS load failed: ' ME.message]);
            data.T1 = zeros(size(data.mag));
            data.aa = zeros(size(data.mag));
        end

        % If exists, load CSF geodesic distance (might need flip/rotation)
        try
            load(dmatfile, 'D')
            data.dist = D;
        catch
            data.dist = zeros(size(data.mag));
        end

        % Global ROI and within-ROI velocities
        if data.hasVel
            data.groi = data.roi > 0;
            data.ginds = find(data.groi(:));
            data.imap = zeros(size(data.groi));
            data.imap(data.ginds) = 1:numel(data.ginds);
            data.imap = flip(data.imap, 2);
        else
            direction.pca = [0; 0; 1];
            direction.man = [0; 0; 1];
            DIRMODE = 'man';
            dirBtn.Text = 'Direction: Manual';
        end

        vstdfile = fullfile(subjectFolder, 'VSTD.nii.gz');
        if isempty(findNiiPath(vstdfile))
            uialert(fig, ['VSTD not found in ' subjectFolder], 'Load error');
            return;
        end
        VSTD = readNiiVol(vstdfile);
        data.vstd = imrotate(VSTD, -90);

        refreshVolMaps();
        updateCubeBtnState();
        updateSegBtnText();

        % OPTIONAL DATA IN AX21
        data.VENTS = ismember(data.aa, [14 15 4 5 44 45]);
        disp(['AX21: ' AX21])
        if strcmp(AX21, 'T1')
            data.ax21 = data.T1; 
        elseif  strcmp(AX21, 'FST1')
            data.ax21 = data.T1 + ( data.VENTS * 0.5 * max(data.T1(:)) ); % TEMP 
        elseif strcmp(AX21, 'VENTS')
            data.ax21 = data.VENTS; 
        elseif strcmp(AX21, 'DIST')
            data.ax21 = data.dist; % CSF geodesic distance (spinal-to-cranial) 
        else
            data.ax21 = data.mag;
        end

        [~, ~, zres] = size(data.VxVSTD);
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

        toggleSaved();

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

        % TEMP: no zoom state preserve
        prevLim = axis(ax(2,1));  % Save zoom state
        cla(ax(2,1));
        imagesc(ax(2,1), squeeze(data.ax21(:,:,slice))');
        if ~isequal(prevLim, [0 1 0 1])  % If zoomed, restore view
            axis(ax(2,1), prevLim);
        else
            axis(ax(2,1), 'image');
        end
        title(ax(2,1), AX21); 

        prevLim = axis(ax(2,2));  % Save zoom state
        cla(ax(2,2));
        imagesc(ax(2,2), squeeze(data.CUBE(:,:,slice))');
        if ~isequal(prevLim, [0 1 0 1])  % If zoomed, restore view
            axis(ax(2,2), prevLim);
        else
            axis(ax(2,2), 'image');
        end
        title(ax(2,2), anatTitleWithShift()); 

        prevLim = axis(ax(2,3));  % Save zoom state
        cla(ax(2,3));
        vstdSl = squeeze(data.VxVSTD(:,:,slice))';
        vstdClip = prctile(data.VxVSTD(:), 99);
        if isempty(vstdClip) || ~isfinite(vstdClip) || vstdClip <= 0
            vstdClip = max(data.VxVSTD(:), [], 'omitnan');
        end
        if isempty(vstdClip) || ~isfinite(vstdClip) || vstdClip <= 0
            vstdClip = 1;
        end
        vstdSl(vstdSl > vstdClip) = vstdClip;
        imagesc(ax(2,3), vstdSl);
        if ~isequal(prevLim, [0 1 0 1])  % If zoomed, restore view
            axis(ax(2,3), prevLim);
        else
            axis(ax(2,3), 'image');
        end

        % Allow pointer to follow slice scroll / refresh patch at click
        if ~isempty(clickedX) && ~isempty(clickedY)
            updateWaveformsFromCoords(clickedX, clickedY, slice);
        end
        title(ax(2,3), vstdTitle());
        clim(ax(2,3), [0, 0.5 * vstdClip]);

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
        if velocitiesLoaded() && isfield(data, 'proj2D') && ~isempty(data.proj2D)
            frame = currentFrame;
            imagesc(ax(1,4), flip(data.proj2D(:,:,frame)', 1));
            colormap(ax(1,4), gray);
            d2d = data.proj2D .* data.bseg;
            CSCALE = 0.01*str2double(cscDropdown.Value);
            clim(ax(1,4), [CSCALE * min(d2d(:)), CSCALE * max(d2d(:))]);
            hold(ax(1,4), 'on');
            visboundaries(ax(1,4), flip(data.bseg', 1), 'Color', 'r', 'LineWidth', 1.5);
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

        imagesc(ax(1,4), flip(img', 1));
        axis(ax(1,4), 'image');
        colormap(ax(1,4), gray);
        title(ax(1,4), sprintf('Flow plane (%s), frame %d', toggle.Value, frame));
        hold(ax(1,4), 'on');
        visboundaries(ax(1,4), flip(data.bseg', 1), 'Color', 'r', 'LineWidth', 1.5);
        hold(ax(1,4), 'off');
        d2d = data.proj2D .* data.bseg;
        CSCALE = 0.01*str2double(cscDropdown.Value);
        clim(ax(1,4), [CSCALE * min(d2d(:)), CSCALE * max(d2d(:))]);
    end

    function saveWaveforms(filename)

        if strcmp(CUBEMODE, 'AF')
            filename = ['AF-' filename];
        end

        if ~velocitiesLoaded()
            uialert(fig, 'Load mode None: waveforms not available.', 'Save Error');
            return;
        end

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

        [xres, yres, zres, ~] = size(data.vstd);
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
        WF.size = str2num(szDropdown.Value); %#ok<*ST2NM>
        WF.thr = str2num(thrDropdown.Value);
        WF.dilate = str2double(dilateDropdown.Value);

        save(fullfile(savefolder, filename), 'WF');
        
        disp(['Saved waveforms to ', fullfile(savefolder, filename)]);
    end

    function updateWaveformsFromClick(event, slice)
        pos = round(event.IntersectionPoint(1:2));
        clickedX = pos(1);
        clickedY = pos(2);
        updateWaveformsFromCoords(clickedX, clickedY, slice);
    end

    function plotClickMarker(x, y)
        if ~isempty(marker)
            delete(marker(ishandle(marker)));
        end
        marker = gobjects(1, 3);
        for j = 1:3
            hold(ax(2,j), 'on');
            marker(j) = plot(ax(2,j), x, y, 'ro', 'MarkerSize', 4, ...
                'MarkerFaceColor', 'r', 'LineWidth', 1);
            hold(ax(2,j), 'off');
        end
    end

    function showPatchPanels(cube_patch, patch, bseg)
        imagesc(ax(1,2), flip(cube_patch', 1));
        axis(ax(1,2), 'image');
        title(ax(1,2), ['Flow plane: ' anatTitle()]);
        colormap(ax(1,2), gray);
        hold(ax(1,2), 'on');
        visboundaries(ax(1,2), flip(bseg', 1), 'Color', 'r', 'LineWidth', 1.5);
        hold(ax(1,2), 'off');
        ax(1,2).XTick = [];
        ax(1,2).YTick = [];

        imagesc(ax(1,3), flip(patch', 1));
        axis(ax(1,3), 'image');
        title(ax(1,3), ['Flow plane: ' segVolLabel()]);
        colormap(ax(1,3), gray);
        hold(ax(1,3), 'on');
        visboundaries(ax(1,3), flip(bseg', 1), 'Color', 'r', 'LineWidth', 1.5);
        hold(ax(1,3), 'off');
        ax(1,3).XTick = [];
        ax(1,3).YTick = [];
    end

    % This needs to be called for within ROI data
    function updateWaveformsFromCoords(x, y, z, bsegManual)

        if nargin < 4
            bsegManual = [];
        end

        [xres, yres, zres] = size(data.vstd);
        if any([x y z] < 1) || x > xres || y > yres || z > zres
            disp('Point outside volume');
            return;
        end

        hasVel = velocitiesLoaded();
        pc1 = [];
        PC1 = [];

        if hasVel
            ind = data.imap(x, y, z); % ind within ROI

            if ind == 0
                disp('Point outside ROI')
                return;
            end

            % Voxel (in current coordinate) waveforms 
            vx_t = squeeze(data.vx(ind, :))';
            vy_t = squeeze(data.vy(ind, :))';
            vz_t = squeeze(data.vz(ind, :))';

            % Volume PCA (shouldnt help alot when regularization is high)
            dilmap = zeros(size(data.imap));
            dilmap(x, y, z) = 1;
            dilmap = imdilate(dilmap, strel('sphere', 2));
            inds = data.imap(find(dilmap));
            inds(inds==0) = [];

            % Speherical volume (around voxel) waveforms
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
            V = [VX_t VY_t VZ_t]; % (frames x (voxels x 3))
            [~, score, ~] = pca(V, 'NumComponents', 1);
            PC1 = score(:,1); % Time series
        else
            cla(ax(3,1), 'reset');
        end

        patch_width = str2double(szDropdown.Value);
        local_thresh = str2double(thrDropdown.Value);
        if strcmpi(clipDropdown.Value, 'Off')
            clip_on = false;
            clip_val = 99;
        else
            clip_on = true;
            clip_val = str2double(clipDropdown.Value);
        end
        dilate_val = str2double(dilateDropdown.Value);

        % *** SPECIFY SEG/DIR/SHAPE MODES IN TOP OF GUI *** 
        [flow, cube_patch, patch, bseg, ~, proj2D] = ... 
            extractThroughPlaneFlow_V3D(data, [x, y, z], direction, patch_width, segVolField(), local_thresh, DIRMODE, SHAPEMODE, clip_on, clip_val, dilate_val, bsegManual, LOCALCSMODE);

        data.patch = patch;
        data.bseg = bseg;

        if hasVel
            data.proj2D = proj2D.proj;
            data.velx2D = proj2D.velx;
            data.vely2D = proj2D.vely;
            data.velz2D = proj2D.velz;
            data.flow = flow;
            data.pc1 = pc1;
            data.PC1 = PC1;
        else
            data.proj2D = [];
            data.velx2D = [];
            data.vely2D = [];
            data.velz2D = [];
            data.flow = [];
            cla(ax(1,4), 'reset');
            cla(ax(3,2), 'reset');
        end

        % Clear old arrow if it exists
        if isgraphics(flowArrow)
            delete(flowArrow);
            flowArrow = [];
        end
        if ~isempty(marker)
            delete(marker(ishandle(marker)));
        end
        marker = [];
        if hasVel
            dir_xy = direction.(DIRMODE)([2, 3]);
            arrowLength = 15;
            if norm(dir_xy) > 0
                dir_xy = dir_xy / norm(direction.(DIRMODE)) * arrowLength;
            else
                dir_xy = [0 arrowLength];
            end
            for j = 1:3, hold(ax(2,j), 'on'); end
            marker = [ ...
                quiver(ax(2,1), x, y, dir_xy(1), dir_xy(2), 'Color', 'r', 'LineWidth', 2.5, 'AutoScale', 'off', 'MaxHeadSize', 2.5);
                quiver(ax(2,2), x, y, dir_xy(1), dir_xy(2), 'Color', 'r', 'LineWidth', 2.5, 'AutoScale', 'off', 'MaxHeadSize', 2.5);
                quiver(ax(2,3), x, y, dir_xy(1), dir_xy(2), 'Color', 'r', 'LineWidth', 2.5, 'AutoScale', 'off', 'MaxHeadSize', 2.5)];
            for j = 1:3, hold(ax(2,j), 'off'); end
        else
            plotClickMarker(x, y);
        end

        showPatchPanels(cube_patch, patch, bseg);

        if hasVel
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
                legend(ax(3,2), 'Voxel PCA', 'Volume PCA', 'Through-plane flow')
            elseif strcmp(WFSHOWMODE, 'FLOW') % Just show the through-plane flow rate waveforms
                hold(ax(3,2), 'on'); % need to reapply hold on for each axis?
                plot(ax(3,2), flow, 'Color', 'r', 'LineWidth', 1.5);
                ylabel(ax(3,2), 'Flow rate (ml/s)');
                title(ax(3,2), sprintf('PC-1 dir: [%.2f  %.2f  %.2f]', direction.pca(1), direction.pca(2), direction.pca(3)));

                atpf = max(flow) - min(flow);
                legend(ax(3,2), ['TPF-amp: ' num2str(atpf)]);
            end

            % Plot coronal view in ax(2,4)
            cla(ax(2,4), 'reset');
            cor = squeeze(data.CUBE(clickedX, :, :)); % could use data.MIXED
            nz = size(cor, 2);
            ys = clickedY-nz/2; ye = clickedY+nz/2;
            if ye > size(data.CUBE, 2)
                ye = size(data.CUBE, 2);
            end
            cor = cor(ys:ye, :);
            imagesc(ax(2,4), cor);
            axis(ax(2,4), 'image');
            colormap(ax(2,4), gray);
            title(ax(2,4), [vstdTitle() ' (coronal)']);
            clim(ax(2,4), [0, 1.0 * max(data.CUBE(:))]);
            hold(ax(2,4), 'on');

            % Add velocity direction as quiver
            dir_xz = direction.(DIRMODE)([1, 3]);  % vx and vz components
            dir_xz = dir_xz / norm(direction.pca) * 10;  % scale for visibility

            FLIPDXZ = -1; % Based on FMo, it seems like this arrow needs to be flipped?
            quiver(ax(2,4), z, (nz/2) - 1, FLIPDXZ*dir_xz(1), dir_xz(2), ... 
                'Color', 'r', 'LineWidth', 1.5, 'AutoScale', 'off', 'MaxHeadSize', 1.5);

            hold(ax(2,4), 'off');
            ax(2,4).XTick = [];
            ax(2,4).YTick = [];

            updateFlowPlane();
        end

    end

    % TEMP: Voxel/Volume PCA or flow rate for propagation delay? 
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

    % TODO this clears correctly? not a toggle, just a show?  
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

            % TEMP: Flip sign if negatively correlated with SC? 
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

    % Feb 2026: simple export of current ROI data with user-defined name
    function showOverlay2D()
        if isempty(fieldnames(data))
            uialert(fig, 'Load data first.', 'Show Volumes');
            return;
        end
        ShowVolumes(fig, data, round(s_slice.Value));
    end

    function exportProj()
        proj2D = [];
        proj2D.tvel = data.proj2D;
        proj2D.patch = data.patch;
        proj2D.bseg = data.bseg;
        proj2D.flow = data.flow'; 
        proj2D.pc1 = data.pc1;
        if ~exist('EXPORT', 'dir')
            mkdir('EXPORT');
        end
        savename = input('Enter save name (e.g., CA): ', 's');
        fname = fullfile('EXPORT', sprintf('%s-%s.mat', savename, subjname));
        save(fname, 'proj2D');
    end

    function toggleLOADMODE()
        if strcmp(LOADMODE, 'MASKED')
            LOADMODE = 'FULL';
            mvBtn.Text = 'Load Mode: Full';
        elseif strcmp(LOADMODE, 'FULL')
            LOADMODE = 'NONE';
            mvBtn.Text = 'Load Mode: None';
        else
            LOADMODE = 'MASKED';
            mvBtn.Text = 'Load Mode: Masked';
        end
    end

    function toggleDIRMODE()
        if strcmp(DIRMODE, 'pca')
            DIRMODE = 'man';
            dirBtn.Text = 'Direction: Manual';
        else
            DIRMODE = 'pca';
            dirBtn.Text = 'Direction: PCA';
        end

        % Update flow-plane and waveforms if a point is selected
        if ~isempty(clickedX) && ~isempty(clickedY)
            updateWaveformsFromCoords(clickedX, clickedY, round(s_slice.Value));
            if velocitiesLoaded()
                updateFlowPlane();
            end
        end
    end

    function toggleSHAPEMODE()
        if strcmp(SHAPEMODE, 'CIRC')
            SHAPEMODE = 'RECT';
            shapeBtn.Text = 'Shape: Rect.';
        elseif strcmp(SHAPEMODE, 'RECT')
            SHAPEMODE = 'ANY';
            shapeBtn.Text = 'Shape: Any';
        elseif strcmp(SHAPEMODE, 'ANY')
            SHAPEMODE = 'CIRC';
            shapeBtn.Text = 'Shape: Circular';
        end

        % Update flow-plane and waveforms if a point is selected
        if ~isempty(clickedX) && ~isempty(clickedY)
            updateWaveformsFromCoords(clickedX, clickedY, round(s_slice.Value));
        end
    end

    function toggleLOCALCSMODE()
        if strcmp(LOCALCSMODE, 'FULL')
            LOCALCSMODE = 'CENTRAL';
            localCsBtn.Text = 'Local CS: Central';
        elseif strcmp(LOCALCSMODE, 'CENTRAL')
            LOCALCSMODE = 'LARGEST';
            localCsBtn.Text = 'Local CS: Largest';
        else
            LOCALCSMODE = 'FULL';
            localCsBtn.Text = 'Local CS: Full';
        end

        if ~isempty(clickedX) && ~isempty(clickedY)
            updateWaveformsFromCoords(clickedX, clickedY, round(s_slice.Value));
        end
    end

    function toggleSEGMODE()
        if strcmp(SEGMODE, 'VxVSTD')
            SEGMODE = 'ANAT';
        else
            SEGMODE = 'VxVSTD';
        end
        updateSegBtnText();
        if ~isempty(fieldnames(data))
            updateDisplays();
        end

        % Update flow-plane and waveforms if a point is selected
        if ~isempty(clickedX) && ~isempty(clickedY)
            updateWaveformsFromCoords(clickedX, clickedY, round(s_slice.Value));
        end
    end

    function toggleCUBEAF()
        if ~isfield(data, 'hasAF') || ~data.hasAF
            return
        end
        if strcmp(CUBEMODE, 'CUBE')
            CUBEMODE = 'AF';
        else
            CUBEMODE = 'CUBE';
        end
        refreshVolMaps();
        updateCubeBtnState();
        updateSegBtnText();
        updateDisplays();
        updateFlowPlane();
    end

    function toggleWFSHOWMODE()
        if strcmp(WFSHOWMODE, 'ALL')
            WFSHOWMODE = 'FLOW';
            wfsBtn.Text = 'WF-Vis: Flow only';
        elseif strcmp(WFSHOWMODE, 'FLOW')
            WFSHOWMODE = 'ALL';
            wfsBtn.Text = 'WF-Vis: Flow + PCA';
        end

        % Update flow-plane and waveforms if a point is selected
        if velocitiesLoaded() && ~isempty(clickedX) && ~isempty(clickedY)
            updateWaveformsFromCoords(clickedX, clickedY, round(s_slice.Value));
        end
    end

    function setRotation(idx, val)
        rotation(idx) = val;  % update rotation array

        % Normalize vector to length 1
        direction.man = rotation / norm(rotation);
        direction.man = direction.man(:);

        % If a point is selected, update flow-plane
        if ~isempty(clickedX) && ~isempty(clickedY)
            updateWaveformsFromCoords(clickedX, clickedY, round(s_slice.Value));
        end
    end

end
