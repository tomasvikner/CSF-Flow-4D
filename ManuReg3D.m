function data = ManuReg3D(parentFig, data, axRef, sliceIdx, movingField, displayField, center, direction, dirMode, patchWidth)
% Manual 3D registration: shift T2CUBE or AF toward MAG with arrow / slice controls.
% Panel 1: sagittal slice (zoomed to main GUI view) — vol(:,:,slice)'.
% Panel 2: local cross-section at center (current flow direction).
% Panel 3: axial at mid LR (dim 2) — vol(:,midLR,:)'.
% Nudge axes (volume AP, LR, SI = dims 1, 2, 3):
%   Up/Down -> LR, Left/Right -> AP, SI-/SI+ -> SI

    if nargin < 6 || isempty(displayField)
        displayField = movingField;
    end
    if nargin < 7
        center = [];
    end
    if nargin < 8
        direction = [];
    end
    if nargin < 9 || isempty(dirMode)
        dirMode = 'man';
    end
    if nargin < 10 || isempty(patchWidth)
        patchWidth = 21;
    end
    if ~isfield(data, 'mag') || isempty(data.mag)
        uialert(parentFig, 'MAG not loaded.', 'ManuReg-3D');
        return;
    end
    if ~isfield(data, movingField) || isempty(data.(movingField))
        uialert(parentFig, [movingField ' not loaded.'], 'ManuReg-3D');
        return;
    end

    magVol = volume3D(data.mag);
    workVol = data.(movingField);
    nz = min(size(magVol, 3), size(workVol, 3));
    sliceIdx = max(1, min(nz, round(sliceIdx)));

    magSl = sagittalDisplay(magVol, sliceIdx);
    [yIdx, xIdx] = zoomRoi(axRef, size(magSl));

    S = struct();
    S.magVol = magVol;
    S.workVol = workVol;
    S.vstd = data.vstd;
    S.sliceIdx = sliceIdx;
    S.nz = nz;
    S.yIdx = yIdx;
    S.xIdx = xIdx;
    S.viewMode = 'MAG';
    S.movingField = movingField;
    S.displayField = displayField;
    S.shiftBase = regShiftVec(data, movingField);
    S.shiftSession = [0 0 0];
    S.center = center;
    S.direction = direction;
    S.dirMode = dirMode;
    S.patchWidth = patchWidth;
    S.midLR = round(size(magVol, 2) / 2);

    regFig = uifigure('Name', sprintf('ManuReg-3D | %s', movingField), ...
        'Position', [60 150 1360 520], 'CloseRequestFcn', @(~,~) closeReg());
    g = uigridlayout(regFig, [1, 4]);
    g.ColumnWidth = {'1x', '1x', '1x', 130};
    g.Padding = [8 8 8 8];

    axSagittal = uiaxes(g);
    axSagittal.Layout.Column = 1;
    colormap(axSagittal, gray);

    axLocal = uiaxes(g);
    axLocal.Layout.Column = 2;
    colormap(axLocal, gray);

    axAxialMid = uiaxes(g);
    axAxialMid.Layout.Column = 3;
    colormap(axAxialMid, gray);

    ctrl = uigridlayout(g, [11, 1]);
    ctrl.Layout.Column = 4;
    ctrl.RowHeight = {30, 30, 30, 30, 30, 24, 24, 30, 30, 30, '1x'};
    ctrl.RowSpacing = 4;
    ctrl.Padding = [0 0 0 0];

    toggleBtn = uibutton(ctrl, 'Text', 'Show: MAG', 'ButtonPushedFcn', @(~,~) toggleView());
    uibutton(ctrl, 'Text', 'Up', 'ButtonPushedFcn', @(~,~) nudgeLR(-1));
    uibutton(ctrl, 'Text', 'Down', 'ButtonPushedFcn', @(~,~) nudgeLR(1));
    uibutton(ctrl, 'Text', 'Left', 'ButtonPushedFcn', @(~,~) nudgeAP(-1));
    uibutton(ctrl, 'Text', 'Right', 'ButtonPushedFcn', @(~,~) nudgeAP(1));
    uibutton(ctrl, 'Text', 'SI-', 'FontSize', 10, 'ButtonPushedFcn', @(~,~) nudgeSI(-1));
    uibutton(ctrl, 'Text', 'SI+', 'FontSize', 10, 'ButtonPushedFcn', @(~,~) nudgeSI(1));
    uilabel(ctrl, 'Text', sprintf('View slice %d/%d', sliceIdx, nz), ...
        'HorizontalAlignment', 'center');
    lblShift = uilabel(ctrl, 'Text', sprintf('Shift: %s', ...
        formatRegShift(S.shiftBase)), ...
        'HorizontalAlignment', 'center', 'FontSize', 10, 'WordWrap', 'on');
    uibutton(ctrl, 'Text', 'Done', 'ButtonPushedFcn', @(~,~) closeReg());

    S.toggleBtn = toggleBtn;
    S.lblShift = lblShift;
    regFig.UserData = S;
    redraw();
    uiwait(regFig);

    if isvalid(regFig)
        S = regFig.UserData;
        data.(movingField) = S.workVol;
        data = addRegShift(data, movingField, S.shiftSession);
        delete(regFig);
    end

    function toggleView()
        S = regFig.UserData;
        modes = {'MAG', 'SEG', 'OVERLAY'};
        idx = find(strcmp(S.viewMode, modes), 1);
        S.viewMode = modes{mod(idx, numel(modes)) + 1};
        regFig.UserData = S;
        redraw();
    end

    function nudgeAP(delta)
        nudge(delta, 0, 0);
    end

    function nudgeLR(delta)
        nudge(0, delta, 0);
    end

    function nudgeSI(delta)
        nudge(0, 0, delta);
    end

    function nudge(d1, d2, d3)
        S = regFig.UserData;
        S.workVol = circshift(S.workVol, [d1, d2, d3]);
        S.shiftSession = S.shiftSession + [d1, d2, d3];
        S.lblShift.Text = sprintf('Shift: %s', ...
            formatRegShift(S.shiftBase + S.shiftSession));
        regFig.UserData = S;
        redraw();
        disp(sprintf('ManuReg-3D %s shift total: %s', ...
            movingField, formatRegShift(S.shiftBase + S.shiftSession)));
    end

    function redraw()
        S = regFig.UserData;
        redrawSagittal(axSagittal, S);
        redrawLocal(axLocal, S);
        redrawAxial(axAxialMid, S);
        regFig.UserData = S;
    end

    function redrawSagittal(axh, S)
        cla(axh);
        switch S.viewMode
            case 'MAG'
                im = cropSlice(sagittalDisplay(S.magVol, S.sliceIdx), S.yIdx, S.xIdx);
                imagesc(axh, im);
                colormap(axh, gray);
                S.toggleBtn.Text = 'Show: MAG';
                title(axh, sprintf('Sagittal MAG | slice %d/%d', S.sliceIdx, S.nz));
            case 'SEG'
                im = cropSlice(sagittalDisplay(displayVol(S), S.sliceIdx), S.yIdx, S.xIdx);
                imagesc(axh, im);
                colormap(axh, gray);
                S.toggleBtn.Text = 'Show: Seg';
                title(axh, sprintf('Sagittal %s | slice %d/%d', S.displayField, S.sliceIdx, S.nz));
            case 'OVERLAY'
                magIm = cropSlice(sagittalDisplay(S.magVol, S.sliceIdx), S.yIdx, S.xIdx);
                segIm = cropSlice(sagittalDisplay(displayVol(S), S.sliceIdx), S.yIdx, S.xIdx);
                if ~isequal(size(magIm), size(segIm))
                    segIm = imresize(segIm, size(magIm));
                end
                image(axh, overlaySlice(magIm, segIm));
                S.toggleBtn.Text = 'Show: Overlay';
                title(axh, sprintf('Sagittal MAG + %s | slice %d/%d', ...
                    S.displayField, S.sliceIdx, S.nz));
        end
        axis(axh, 'image');
        axh.XTick = [];
        axh.YTick = [];
    end

    function redrawLocal(axh, S)
        cla(axh);
        if isempty(S.center) || numel(S.center) < 3 || isempty(S.direction) ...
                || ~isfield(S.direction, S.dirMode) || isempty(S.direction.(S.dirMode))
            title(axh, 'Local CS (click a point in main GUI)');
            axis(axh, 'image');
            axh.XTick = [];
            axh.YTick = [];
            return;
        end

        centerVec = S.center(:);
        try
            magPatch = extractLocalPatchInterp(S.magVol, centerVec, S.direction, S.dirMode, S.patchWidth);
            segPatch = extractLocalPatchInterp(displayVol(S), centerVec, S.direction, S.dirMode, S.patchWidth);
        catch ME
            title(axh, 'Local CS unavailable');
            axh.XTick = [];
            axh.YTick = [];
            warning('ManuReg3D:LocalCS', '%s', ME.message);
            return;
        end

        switch S.viewMode
            case 'MAG'
                imagesc(axh, patchDisplay(magPatch));
                colormap(axh, gray);
                title(axh, 'Local CS | MAG');
            case 'SEG'
                imagesc(axh, patchDisplay(segPatch));
                colormap(axh, gray);
                title(axh, sprintf('Local CS | %s', S.displayField));
            case 'OVERLAY'
                magIm = patchDisplay(magPatch);
                segIm = patchDisplay(segPatch);
                image(axh, overlaySlice(magIm, segIm));
                title(axh, sprintf('Local CS | MAG + %s', S.displayField));
        end
        axis(axh, 'image');
        axh.XTick = [];
        axh.YTick = [];
    end

    function redrawAxial(axh, S)
        cla(axh);
        midLR = S.midLR;
        nLR = size(S.magVol, 2);
        switch S.viewMode
            case 'MAG'
                im = cropAxialWide(axialDisplay(S.magVol, midLR));
                imagesc(axh, im);
                colormap(axh, gray);
                title(axh, sprintf('Axial MAG | LR %d/%d', midLR, nLR));
            case 'SEG'
                im = cropAxialWide(axialDisplay(displayVol(S), midLR));
                imagesc(axh, im);
                colormap(axh, gray);
                title(axh, sprintf('Axial %s | LR %d/%d', S.displayField, midLR, nLR));
            case 'OVERLAY'
                magIm = cropAxialWide(axialDisplay(S.magVol, midLR));
                segIm = cropAxialWide(axialDisplay(displayVol(S), midLR));
                if ~isequal(size(magIm), size(segIm))
                    segIm = imresize(segIm, size(magIm));
                end
                image(axh, overlaySlice(magIm, segIm));
                title(axh, sprintf('Axial MAG + %s | LR %d/%d', ...
                    S.displayField, midLR, nLR));
        end
        axis(axh, 'image');
        axh.XTick = [];
        axh.YTick = [];
    end

    function vol = displayVol(S)
        if strcmp(S.displayField, S.movingField)
            vol = S.workVol;
        elseif strcmp(S.displayField, 'VxVSTD')
            vol = S.vstd .* S.workVol;
        else
            vol = S.workVol;
        end
    end

    function closeReg()
        if isvalid(regFig)
            uiresume(regFig);
        end
    end

end

function im = patchDisplay(patch)
    im = flip(patch', 1);
end

function im = sagittalDisplay(vol, sliceIdx)
    % Main GUI ax(2,1–3) orientation — sagittal stack along dim 3.
    im = squeeze(vol(:,:,sliceIdx))';
end

function im = axialDisplay(vol, lrIdx)
    % Orthogonal axial cut at mid LR (dim 2).
    im = squeeze(vol(:, lrIdx, :))';
end

function im = cropAxialWide(im, trimEachSide)
    % Center-crop the wider in-plane axis (e.g. 512x80 -> 256x80).
    if nargin < 2
        trimEachSide = 128;
    end
    [h, w] = size(im);
    if w >= h && w > 2 * trimEachSide
        im = im(:, trimEachSide + 1:w - trimEachSide);
    elseif h > w && h > 2 * trimEachSide
        im = im(trimEachSide + 1:h - trimEachSide, :);
    end
end

function im = cropSlice(im, yIdx, xIdx)
    im = im(yIdx, xIdx);
end

function rgb = overlaySlice(imA, imB)
    a = normSlice(imA);
    b = normSlice(imB);
    rgb = cat(3, a, b, zeros(size(a)));
end

function im = normSlice(im)
    im = double(im);
    lo = prctile(im(:), 1);
    hi = prctile(im(:), 99);
    if hi > lo
        im = (im - lo) / (hi - lo);
    else
        im = zeros(size(im));
    end
    im = max(0, min(1, im));
end

function vol3 = volume3D(vol)
    if ndims(vol) > 3
        vol3 = mean(vol, 4);
    else
        vol3 = vol;
    end
end

function [yIdx, xIdx] = zoomRoi(axh, imSize)
    xlims = xlim(axh);
    ylims = ylim(axh);

    if isequal(xlims, [0 1]) && isequal(ylims, [0 1])
        yIdx = 1:imSize(1);
        xIdx = 1:imSize(2);
        return;
    end

    xIdx = round(max(1, xlims(1))):round(min(imSize(2), xlims(2)));
    yIdx = round(max(1, ylims(1))):round(min(imSize(1), ylims(2)));
end
