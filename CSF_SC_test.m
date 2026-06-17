% CSF_SC.m — Spinal canal (SI) waveform tool
% Loads z-velocity (brainvels), T2 CUBE, and anti-FLAIR; sagittal preview;
% define axial slice range (bottom–top, SI), segment per slice, plot vz waveforms.

function CSF_SC

% ---------------------------- SETTINGS -------------------------------- %
PREVIEW   = true;  % true: T2/AF/MAG only, no velocities or flow waveforms
LOCAL = true;
LOCALVELS = '/Users/txv016/Documents/BRAINVELS';
BASEPATH  = '/Volumes/groups/CVMRIGroup/Users/txv016/WRAP2/niis/CURRENT/';
CLIPVAL   = 99;   % global intensity clip percentile (fixed)
MAX_SEG_SLICES = 100;
addpath('/Users/txv016/Documents/MATLAB/Freesurfer Tools');

% ---------------------------------------------------------------------- %

delete(findall(0, 'Type', 'figure', 'Name', 'CSF Spinal Canal'));

data = struct();
subjname = '';
wfFig = [];
segFigs = gobjects(0);

figW = 1200;
figH = round(figW / 3);
scr = groot().ScreenSize;
figPos = [scr(1) + (scr(3) - figW) / 2, scr(2) + (scr(4) - figH) / 2, figW, figH];

fig = uifigure('Name', 'CSF Spinal Canal', 'Position', figPos);
g = uigridlayout(fig, [1, 4]);
g.ColumnWidth = {'1x', '1x', '1x', 100};
g.Padding = [4 4 4 4];
g.ColumnSpacing = 4;

ax = gobjects(1, 3);
panel3Title = 'MAG';
if ~PREVIEW
    panel3Title = '|vz| std'; %#ok<*UNRCH>
end
titles = {'T2 CUBE', 'Anti-FLAIR', panel3Title};
for j = 1:3
    ax(j) = uiaxes(g);
    ax(j).Layout.Column = j;
    colormap(ax(j), gray);
    ax(j).XTick = [];
    ax(j).YTick = [];
    title(ax(j), titles{j});
end
linkaxes(ax, 'xy');

ctrl = uigridlayout(g, [16, 1]);
ctrl.Layout.Column = 4;
ctrl.RowHeight = repmat({26}, 1, 16);
ctrl.RowSpacing = 4;
ctrl.Padding = [0 0 0 0];

uibutton(ctrl, 'Text', 'Load subject', 'ButtonPushedFcn', @(~,~) loadSubject());
uibutton(ctrl, 'Text', 'Segment range', 'ButtonPushedFcn', @(~,~) runSegmentation());
uibutton(ctrl, 'Text', 'Close', 'ButtonPushedFcn', @(~,~) close(fig));

uilabel(ctrl, 'Text', 'Preview sag', 'HorizontalAlignment', 'center');
s_preview = uislider(ctrl, 'Limits', [1 10], 'MajorTicks', [], ...
    'ValueChangedFcn', @(~,~) updateDisplay());

uilabel(ctrl, 'Text', 'Ax bottom', 'HorizontalAlignment', 'center');
s_axTop = uislider(ctrl, 'Limits', [1 10], 'MajorTicks', [], ...
    'ValueChangedFcn', @(~,~) onRangeChanged());

uilabel(ctrl, 'Text', 'Ax top', 'HorizontalAlignment', 'center');
s_axBot = uislider(ctrl, 'Limits', [1 10], 'MajorTicks', [], ...
    'ValueChangedFcn', @(~,~) onRangeChanged());

uilabel(ctrl, 'Text', 'Cor bottom', 'HorizontalAlignment', 'center');
s_corTop = uislider(ctrl, 'Limits', [1 10], 'MajorTicks', [], ...
    'ValueChangedFcn', @(~,~) onRangeChanged());

uilabel(ctrl, 'Text', 'Cor top', 'HorizontalAlignment', 'center');
s_corBot = uislider(ctrl, 'Limits', [1 10], 'MajorTicks', [], ...
    'ValueChangedFcn', @(~,~) onRangeChanged());

uilabel(ctrl, 'Text', 'Threshold %', 'HorizontalAlignment', 'center');
dd_thresh = uidropdown(ctrl, ...
    'Items', {'50', '60', '70', '80', '85', '90', '95'}, ...
    'Value', '80');

lbl_range = uilabel(ctrl, 'Text', 'Range: —', 'HorizontalAlignment', 'center', ...
    'FontSize', 10);
lbl_status = uilabel(ctrl, 'Text', 'Load a subject.', 'HorizontalAlignment', 'center', ...
    'FontSize', 10, 'WordWrap', 'on');

fig.CloseRequestFcn = @(~,~) closeApp();

    function loadSubject()
        basefolder = uigetdir();
        if basefolder == 0, return; end

        [~, subjname] = fileparts(basefolder);
        if PREVIEW
            fig.Name = ['CSF Spinal Canal (PREVIEW) | ' subjname];
        else
            fig.Name = ['CSF Spinal Canal | ' subjname];
        end
        disp(['CSF_SC: loading ' subjname]);

        subjectFolder = fullfile(BASEPATH, subjname, 'processed');
        magfile = fullfile(subjectFolder, 'MAG.nii');
        cubefile = fullfile(subjectFolder, 'cr2mag', 'r4dT2.nii');
        antiflairfile = fullfile(subjectFolder, 'cr2mag', 'r4dAF.nii');

        if LOCAL && ~PREVIEW
            rvelsfolder = fullfile(LOCALVELS, subjname, 'brainvels');
            disp('Loading velocities from LOCAL brainvels folder');
        elseif ~PREVIEW
            rvelsfolder = fullfile(subjectFolder, 'brainvels');
            disp('Loading velocities from REMOTE brainvels folder');
        end

        if ~exist(magfile, 'file')
            uialert(fig, ['MAG not found: ' magfile], 'Load error');
            return;
        end
        if ~exist(cubefile, 'file')
            uialert(fig, ['T2 CUBE not found: ' cubefile], 'Load error');
            return;
        end
        mag = MRIread(magfile).vol;
        data.T2CUBE = MRIread(cubefile).vol;
        disp('Loaded T2CUBE (r4dT2.nii)');
        mag = imrotate(mag, -90);
        data.T2CUBE = imrotate(data.T2CUBE, -90);
        data.MAG = mag;
        if ~isequal(size(mag), size(data.T2CUBE))
            uialert(fig, 'MAG and T2 CUBE size mismatch.', 'Load error');
            return;
        end

        data.hasAF = false;
        if exist(antiflairfile, 'file')
            try
                data.AF = MRIread(antiflairfile).vol;
                data.AF = imrotate(data.AF, -90);
                data.hasAF = true;
                disp('Loaded anti-FLAIR (r4dAF.nii)');
            catch ME
                disp(['anti-FLAIR load failed: ' ME.message]);
            end
        else
            disp('No anti-FLAIR found (r4dAF.nii)');
        end

        if PREVIEW
            disp('PREVIEW mode: skipping velocity load');
            data.hasVel = false;
        else
            if ~exist(fullfile(rvelsfolder, 'rz.mat'), 'file')
                uialert(fig, ['rz not found: ' rvelsfolder], 'Load error');
                return;
            end
            if ~exist(fullfile(rvelsfolder, 'roi.mat'), 'file')
                uialert(fig, ['roi not found: ' rvelsfolder], 'Load error');
                return;
            end
            S = load(fullfile(rvelsfolder, 'rz.mat'), 'rz');
            R = load(fullfile(rvelsfolder, 'roi.mat'), 'roi');
            data.vz = S.rz;
            data.roi = R.roi;
            if ~isequal(size(data.roi), size(data.T2CUBE))
                uialert(fig, 'roi.mat size does not match T2 CUBE.', 'Load error');
                return;
            end
            disp('rz velocity loaded');
            data = sc_buildImap(data);
            data.vzStd = sc_vzStdMap(data);
            data.hasVel = true;
        end

        data = sc_setVolumeAxes(data);

        sagLim = [1 data.nSag];
        axLim = [1 data.nAxial];
        corLim = [1 data.nCoronal];
        s_preview.Limits = sagLim;
        s_axTop.Limits = axLim;
        s_axBot.Limits = axLim;
        s_corTop.Limits = corLim;
        s_corBot.Limits = corLim;
        s_preview.Value = round(data.nSag / 2);
        s_axTop.Value = min(data.nAxial, round(data.nAxial * 0.75));
        s_axBot.Value = max(1, round(data.nAxial * 0.25));
        s_corTop.Value = min(data.nCoronal, round(data.nCoronal * 0.55));
        s_corBot.Value = max(1, round(data.nCoronal * 0.45));

        if PREVIEW
            lbl_status.Text = sprintf('%s | PREVIEW | sag %d ax %d cor %d', ...
                subjname, data.nSag, data.nAxial, data.nCoronal);
        else
            lbl_status.Text = sprintf('%s | %d fr | sag %d ax %d cor %d', ...
                subjname, size(data.vz, 2), data.nSag, data.nAxial, data.nCoronal);
        end
        updateDisplay();
        onRangeChanged();
        disp('CSF_SC: load complete.');
    end

    function [axLo, axHi, corLo, corHi] = getRanges()
        axLo = min(round(s_axTop.Value), round(s_axBot.Value));
        axHi = max(round(s_axTop.Value), round(s_axBot.Value));
        corLo = min(round(s_corTop.Value), round(s_corBot.Value));
        corHi = max(round(s_corTop.Value), round(s_corBot.Value));
    end

    function onRangeChanged()
        if ~isfield(data, 'T2CUBE'), return; end
        [axLo, axHi, corLo, corHi] = getRanges();
        lbl_range.Text = sprintf('Ax %d–%d | Cor %d–%d', axLo, axHi, corLo, corHi);
        updateDisplay();
    end

    function updateDisplay()
        if ~isfield(data, 'T2CUBE'), return; end
        sagIdx = round(s_preview.Value);
        [axLo, axHi, corLo, corHi] = getRanges();

        t2sag = sc_getPlaneSlice(data.T2CUBE, data.sagDim, sagIdx);
        sc_showSagittalWithRanges(ax(1), t2sag, axLo, axHi, corLo, corHi, data);
        title(ax(1), sprintf('T2 CUBE  sagittal %d/%d', sagIdx, data.nSag));

        if data.hasAF
            afsag = sc_getPlaneSlice(data.AF, data.sagDim, sagIdx);
            sc_showSagittalWithRanges(ax(2), afsag, axLo, axHi, corLo, corHi, data);
            title(ax(2), sprintf('Anti-FLAIR  sagittal %d/%d', sagIdx, data.nSag));
        else
            cla(ax(2));
            axis(ax(2), 'on');
            text(ax(2), 0.5, 0.5, 'No anti-FLAIR', 'HorizontalAlignment', 'center');
            axis(ax(2), 'off');
            title(ax(2), 'Anti-FLAIR');
        end

        if PREVIEW
            magsag = sc_getPlaneSlice(data.MAG, data.sagDim, sagIdx);
            sc_showSagittalWithRanges(ax(3), magsag, axLo, axHi, corLo, corHi, data);
            clim(ax(3), sc_clim(data.MAG));
            title(ax(3), sprintf('MAG  sagittal %d/%d', sagIdx, data.nSag));
        else
            vzsag = sc_getPlaneSlice(data.vzStd, data.sagDim, sagIdx);
            sc_showSagittalWithRanges(ax(3), vzsag, axLo, axHi, corLo, corHi, data);
            clim(ax(3), sc_clim(data.vzStd));
            title(ax(3), sprintf('|vz| std  sagittal %d/%d', sagIdx, data.nSag));
        end
    end

    function runSegmentation()
        if ~isfield(data, 'T2CUBE')
            uialert(fig, 'Load a subject first.', 'Segment');
            return;
        end

        [axLo, axHi, corLo, corHi] = getRanges();
        if axHi < axLo || corHi < corLo
            uialert(fig, 'Invalid slice range.', 'Segment');
            return;
        end

        nSlices = axHi - axLo + 1;
        if nSlices > MAX_SEG_SLICES
            msg = sprintf(['Total axial slice range is %d slices (%d–%d). ' ...
                'Maximum allowed is %d.'], nSlices, axLo, axHi, MAX_SEG_SLICES);
            fprintf('CSF_SC: %s\n', msg);
            uialert(fig, msg, 'Slice range');
            lbl_status.Text = msg;
            return;
        end

        threshPct = 10; % TEMP 
        lbl_status.Text = 'Segmenting...';
        drawnow;

        t2crop = data.T2CUBE(corLo:corHi, axLo:axHi, :);
        afcrop = data.AF(corLo:corHi, axLo:axHi, :);
        thr = [];
        thr.T2 = 0.01 * threshPct * max(t2crop(:));
        thr.AF = 0.01 * threshPct * max(afcrop(:));
        masks = [];
        masks.T2 = t2crop > thr.T2;
        masks.AF = afcrop > thr.AF;

        if ~PREVIEW
            waveforms = zeros(nSlices, size(data.vz, 2));
            dimOff = sc_dimOffsets(data.axialDim, axLo, data.coronalDim, corLo);
            for ia = 1:nSlices
                axIdx = axLo + ia - 1;
                waveforms(ia, :) = sc_sampleVzWaveform(data, masks(:, :, ia), ...
                    data.axialDim, axIdx, dimOff);
            end
            sc_plotSliceWaveforms(axLabels, waveforms, subjname, threshPct, CLIPVAL);
            if ishandle(wfFig)
                close(wfFig);
            end
            wfFig = gcf;
            wfFig.Name = ['SC waveforms | ' subjname];
        else
            disp('PREVIEW mode: segmentation only (no flow waveforms)');
        end

        if ~isempty(segFigs)
            close(segFigs(ishandle(segFigs)));
        end

        % TODO: plot 2D slices with imagesc(image), hold on, visboundaries(mask);

        updateDisplay();

        nPix = squeeze(sum(masks, [1 2]));
        lbl_status.Text = sprintf('Done: ax %d–%d cor %d–%d, thr=%g, %d–%d px/slice', ...
            axLo, axHi, corLo, corHi, globalThr, min(nPix), max(nPix));
        fprintf('CSF_SC: segmented ax %d–%d cor %d–%d, globalThr=%.4g (clip %d%%)\n', ...
            axLo, axHi, corLo, corHi, globalThr, CLIPVAL);
    end

    function closeApp()
        if ishandle(wfFig)
            close(wfFig);
        end
        if ~isempty(segFigs)
            close(segFigs(ishandle(segFigs)));
        end
        delete(fig);
    end

end

% ============================== LOCAL HELPERS ============================== %

function wf = sc_sampleVzWaveform(data, mask2d, axialDim, axIdx, dimOff)
[nx, ny, nz] = size(data.imap);
[rr, cc] = find(mask2d);
if isempty(rr)
    wf = nan(1, size(data.vz, 2));
    return;
end
nVox = numel(rr);
freeDims = setdiff(1:3, axialDim);
subs = cell(1, 3);
subs{axialDim} = repmat(axIdx, nVox, 1);
subs{freeDims(1)} = rr + dimOff(freeDims(1)) - 1;
subs{freeDims(2)} = cc + dimOff(freeDims(2)) - 1;
lin3d = sub2ind([nx, ny, nz], subs{1}, subs{2}, subs{3});
rows = data.imap(lin3d);
valid = rows > 0;
if ~any(valid)
    wf = nan(1, size(data.vz, 2));
    return;
end
wf = mean(double(data.vz(rows(valid), :)), 1, 'omitnan');
end

function data = sc_setVolumeAxes(data)
sz = size(data.T2CUBE);
data.volSize = sz;
data.sagDim = 3; 
data.axialDim = 2;
data.coronalDim = 1;
data.nSag = sz(data.sagDim);
data.nCoronal = sz(data.coronalDim);
data.nAxial = sz(data.axialDim);
fprintf('CSF_SC axes: sagittal dim %d (%d), coronal dim %d (%d), axial dim %d (%d)\n', ...
    data.sagDim, data.nSag, data.coronalDim, data.nCoronal, data.axialDim, data.nAxial);
end

function img = sc_getPlaneSlice(vol, alongDim, idx)
subs = {':', ':', ':'};
subs{alongDim} = idx;
img = squeeze(vol(subs{:}));
end

function lim = sc_clim(vol)
vmax = prctile(vol(isfinite(vol) & vol > 0), 99);
if isempty(vmax) || vmax <= 0
    lim = [0 1];
else
    lim = [0 vmax];
end
end

function sc_plotSliceWaveforms(axLabels, waveforms, subjname, threshPct, clipVal)
nSlices = numel(axLabels);
cmap = parula(max(nSlices, 2));
figure('Color', 'w', 'Position', [120 120 900 420]);
ax = axes;
hold(ax, 'on');
for i = 1:nSlices
    plot(ax, waveforms(i, :), 'Color', cmap(i, :), 'LineWidth', 1.2);
end
hold(ax, 'off');
xlabel(ax, 'Frame');
ylabel(ax, 'vz (SI)');
title(ax, sprintf('%s | SC ROI mean vz | axial %d–%d | thr %g%% clip %g%%', ...
    subjname, axLabels(1), axLabels(end), threshPct, clipVal));
if nSlices <= 20
    legStr = arrayfun(@(k) sprintf('ax=%d', k), axLabels, 'UniformOutput', false);
    legend(ax, legStr, 'Location', 'eastoutside', 'FontSize', 8);
else
    colormap(ax, cmap);
    cb = colorbar(ax);
    cb.Label.String = 'Axial slice';
    caxis(ax, [axLabels(1) axLabels(end)]); %#ok<*CAXIS>
end
grid(ax, 'on');
end
