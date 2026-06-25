% Quick coregistration viewer for 4D flow MAG, VSTD, T2 CUBE, and anti-FLAIR.

function checkCoregs

BASEPATH = '/Volumes/groups/CVMRIGroup/Users/txv016/WRAP2/niis/CURRENT/';
addpath('/Users/txv016/Documents/MATLAB/Freesurfer Tools');

delete(findall(0, 'Type', 'figure', 'Name', 'Check Coregs'));

data = struct();
subjname = '';
viewMode = 'MAG';
planeIdx = 1;
planeNames = {'Axial', 'Sagittal', 'Coronal'};

figW = 960;
figH = 520;
fig = uifigure('Name', 'Check Coregs', 'Position', [200 200 figW figH]);
g = uigridlayout(fig, [1, 2]);
g.ColumnWidth = {'1x', 120};
g.Padding = [4 4 4 4];

ax = uiaxes(g);
ax.Layout.Column = 1;
colormap(ax, gray);
ax.XTick = [];
ax.YTick = [];

ctrl = uigridlayout(g, [10, 1]);
ctrl.Layout.Column = 2;
ctrl.RowHeight = {30, 30, 30, 30, 30, 30, 22, 30, 30, '1x'};
ctrl.RowSpacing = 6;
ctrl.Padding = [0 0 0 0];

uibutton(ctrl, 'Text', 'Load subject', 'ButtonPushedFcn', @(~,~) loadSubject());
btnMag = uibutton(ctrl, 'Text', '4D-Mag', 'ButtonPushedFcn', @(~,~) setView('MAG'));
btnVstd = uibutton(ctrl, 'Text', '4D-VSTD', 'ButtonPushedFcn', @(~,~) setView('VSTD'));
btnCube = uibutton(ctrl, 'Text', 'T2-CUBE', 'ButtonPushedFcn', @(~,~) setView('CUBE'));
btnAf = uibutton(ctrl, 'Text', 'AntiFLAIR', 'ButtonPushedFcn', @(~,~) setView('AF'), 'Enable', 'off');
btnPermute = uibutton(ctrl, 'Text', 'Axial', 'ButtonPushedFcn', @(~,~) togglePlane());
uilabel(ctrl, 'Text', 'Slice', 'HorizontalAlignment', 'center');
s_slice = uislider(ctrl, 'Limits', [1 10], 'MajorTicks', [], ...
    'ValueChangedFcn', @(~,~) updateDisplay());
lblStatus = uilabel(ctrl, 'Text', 'Load a subject.', 'HorizontalAlignment', 'center', ...
    'FontSize', 10, 'WordWrap', 'on');

uibutton(ctrl, 'Text', 'Close GUI', 'ButtonPushedFcn', @(~,~) delete(fig));

fig.CloseRequestFcn = @(~,~) delete(fig);

    function loadSubject()
        basefolder = uigetdir();
        if basefolder == 0, return; end

        [~, subjname] = fileparts(basefolder);
        fig.Name = ['Check Coregs | ' subjname];
        disp(['checkCoregs: loading ' subjname]);

        subjectFolder = fullfile(BASEPATH, subjname, 'processed');
        magfile = fullfile(subjectFolder, 'MAG.nii');
        cubefile = fullfile(subjectFolder, 'cr2mag', 'r4dT2.nii');
        antiflairfile = fullfile(subjectFolder, 'cr2mag', 'r4dAF.nii');
        vstdfile = fullfile(subjectFolder, 'VSTD.nii.gz');

        if ~exist(magfile, 'file')
            uialert(fig, ['MAG not found: ' magfile], 'Load error');
            return;
        end
        if ~exist(cubefile, 'file')
            uialert(fig, ['T2 CUBE not found: ' cubefile], 'Load error');
            return;
        end
        if ~exist(vstdfile, 'file') && isempty(findNiiPath(vstdfile))
            uialert(fig, ['VSTD not found: ' vstdfile], 'Load error');
            return;
        end

        mag = readNiiVol(magfile);
        mag = imrotate(mag, -90);
        data.MAG = to3D(mag);

        data.T2CUBE = readNiiVol(cubefile);
        data.T2CUBE = imrotate(data.T2CUBE, -90);

        data.VSTD = readNiiVol(vstdfile);
        data.VSTD = imrotate(data.VSTD, -90);

        data.hasAF = false;
        if exist(antiflairfile, 'file')
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

        if data.hasAF
            btnAf.Enable = 'on';
        else
            btnAf.Enable = 'off';
        end

        viewMode = 'MAG';
        planeIdx = 1;
        btnPermute.Text = planeNames{planeIdx};
        resetSliceSlider(0.5);

        updateDisplay();
    end

    function setView(mode)
        if isempty(fieldnames(data))
            uialert(fig, 'Load a subject first.', 'Check Coregs');
            return;
        end
        if strcmp(mode, 'AF') && (~isfield(data, 'hasAF') || ~data.hasAF)
            uialert(fig, 'Anti-FLAIR not available for this subject.', 'Check Coregs');
            return;
        end
        viewMode = mode;
        updateDisplay();
    end

    function togglePlane()
        if isempty(fieldnames(data))
            uialert(fig, 'Load a subject first.', 'Check Coregs');
            return;
        end

        frac = (s_slice.Value - s_slice.Limits(1)) / max(1, diff(s_slice.Limits));
        planeIdx = mod(planeIdx, numel(planeNames)) + 1;
        btnPermute.Text = planeNames{planeIdx};
        resetSliceSlider(frac);
        updateDisplay();
    end

    function resetSliceSlider(frac)
        nSl = nSlicesForPlane(currentVolume(), planeIdx);
        s_slice.Limits = [1 nSl];
        s_slice.Value = max(1, min(nSl, round(frac * (nSl - 1) + 1)));
    end

    function updateDisplay()
        if isempty(fieldnames(data)), return; end

        slice = round(s_slice.Value);
        vol = currentVolume();
        im = planeSlice(vol, planeIdx, slice);
        nSl = nSlicesForPlane(vol, planeIdx);

        cla(ax);
        imagesc(ax, im);
        axis(ax, 'image');
        colormap(ax, gray);
        title(ax, sprintf('%s — %s', currentTitle(), planeNames{planeIdx}));

        lblStatus.Text = sprintf('%s | %s %s | slice %d/%d', ...
            subjname, currentTitle(), planeNames{planeIdx}, slice, nSl);
    end

    function vol = currentVolume()
        switch viewMode
            case 'MAG'
                vol = data.MAG;
            case 'VSTD'
                vol = data.VSTD;
            case 'CUBE'
                vol = data.T2CUBE;
            case 'AF'
                vol = data.AF;
            otherwise
                vol = data.MAG;
        end
    end

    function t = currentTitle()
        switch viewMode
            case 'MAG'
                t = '4D-Mag';
            case 'VSTD'
                t = '4D-VSTD';
            case 'CUBE'
                t = 'T2-CUBE';
            case 'AF'
                t = 'AntiFLAIR';
            otherwise
                t = '';
        end
    end

end

function vol3 = to3D(vol)
    if ndims(vol) > 3
        vol3 = mean(vol, 4);
    else
        vol3 = vol;
    end
end

function im = planeSlice(vol, planeIdx, slice)
    switch planeIdx
        case 1 % axial (dim 3)
            im = squeeze(vol(:,:,slice))';
        case 2 % sagittal (dim 1)
            im = squeeze(vol(slice,:,:))';
        case 3 % coronal (dim 2)
            im = squeeze(vol(:,slice,:))';
        otherwise
            im = squeeze(vol(:,:,slice))';
    end
end

function n = nSlicesForPlane(vol, planeIdx)
    switch planeIdx
        case 1
            n = size(vol, 3);
        case 2
            n = size(vol, 1);
        case 3
            n = size(vol, 2);
        otherwise
            n = size(vol, 3);
    end
end
