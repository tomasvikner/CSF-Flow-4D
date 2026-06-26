function ShowVolumes(parentFig, data, sliceIdx)
% Interactive MAG / T2 / Anti-FLAIR overlay viewer with plane and xlab toggles.

    if nargin < 3 || isempty(sliceIdx)
        sliceIdx = [];
    end

    if ~isfield(data, 'mag') || isempty(data.mag)
        uialert(parentFig, 'MAG not loaded.', 'Show Volumes');
        return;
    end
    if ~isfield(data, 'T2CUBE') || isempty(data.T2CUBE)
        uialert(parentFig, 'T2 CUBE not loaded.', 'Show Volumes');
        return;
    end

    oldFig = findall(0, 'Type', 'figure', 'Tag', 'ShowVolumes');
    delete(oldFig);

    magVol = volume3D(data.mag);
    cubeVol = data.T2CUBE;
    cubeXVol = data.T2xVSTD;
    hasAF = isfield(data, 'hasAF') && data.hasAF;
    if hasAF
        afVol = data.AF;
        afXVol = data.AFxVSTD;
    else
        afVol = [];
        afXVol = [];
    end

    S = struct();
    S.magVol = magVol;
    S.cubeVol = cubeVol;
    S.cubeXVol = cubeXVol;
    S.afVol = afVol;
    S.afXVol = afXVol;
    S.hasAF = hasAF;
    S.plane = 'SAGITTAL';
    S.useXlab = false;
    if isempty(sliceIdx)
        sliceIdx = defaultSliceForPlane(magVol, S.plane);
    end
    S.sliceIdx = sliceIdx;

    viewerFig = uifigure('Name', 'Show Volumes', ...
        'Tag', 'ShowVolumes', ...
        'Position', [120 120 1180 520]);
    g = uigridlayout(viewerFig, [1, 4]);
    g.ColumnWidth = {'1x', '1x', '1x', 130};
    g.Padding = [8 8 8 8];

    ax1 = uiaxes(g);
    ax1.Layout.Column = 1;
    ax2 = uiaxes(g);
    ax2.Layout.Column = 2;
    ax3 = uiaxes(g);
    ax3.Layout.Column = 3;

    ctrl = uigridlayout(g, [6, 1]);
    ctrl.Layout.Column = 4;
    ctrl.RowHeight = {30, 30, 30, 30, '1x', 30};
    ctrl.RowSpacing = 6;
    ctrl.Padding = [0 0 0 0];

    xlabBtn = uibutton(ctrl, 'Text', 'Anat: CUBE / AF', ...
        'ButtonPushedFcn', @(~,~) toggleXlab());
    planeBtn = uibutton(ctrl, 'Text', 'Plane: Sagittal', ...
        'ButtonPushedFcn', @(~,~) togglePlane());
    lblSlice = uilabel(ctrl, 'Text', '', 'HorizontalAlignment', 'center');
    sliceSlider = uislider(ctrl, 'Limits', [1 10], 'MajorTicks', [], ...
        'ValueChangedFcn', @(~,~) onSliceChange());
    uilabel(ctrl, 'Text', '');
    uibutton(ctrl, 'Text', 'Close', 'ButtonPushedFcn', @(~,~) delete(viewerFig));

    S.ax = [ax1, ax2, ax3];
    S.xlabBtn = xlabBtn;
    S.planeBtn = planeBtn;
    S.sliceSlider = sliceSlider;
    S.lblSlice = lblSlice;
    viewerFig.UserData = S;

    updateSliceControls();
    redraw();

    function toggleXlab()
        S = viewerFig.UserData;
        S.useXlab = ~S.useXlab;
        if S.useXlab
            S.xlabBtn.Text = 'Anat: xVSTD';
        else
            S.xlabBtn.Text = 'Anat: CUBE / AF';
        end
        viewerFig.UserData = S;
        redraw();
    end

    function togglePlane()
        S = viewerFig.UserData;
        switch S.plane
            case 'SAGITTAL'
                S.plane = 'CORONAL';
                S.planeBtn.Text = 'Plane: Coronal';
            case 'CORONAL'
                S.plane = 'AXIAL';
                S.planeBtn.Text = 'Plane: Axial';
            otherwise
                S.plane = 'SAGITTAL';
                S.planeBtn.Text = 'Plane: Sagittal';
        end
        S.sliceIdx = defaultSliceForPlane(S.magVol, S.plane);
        viewerFig.UserData = S;
        updateSliceControls();
        redraw();
    end

    function onSliceChange()
        S = viewerFig.UserData;
        S.sliceIdx = round(S.sliceSlider.Value);
        viewerFig.UserData = S;
        updateSliceLabel();
        redraw();
    end

    function updateSliceControls()
        S = viewerFig.UserData;
        nz = nPlaneSlices(S.magVol, S.plane);
        S.sliceIdx = max(1, min(nz, round(S.sliceIdx)));
        S.sliceSlider.Limits = [1, nz];
        S.sliceSlider.Value = S.sliceIdx;
        viewerFig.UserData = S;
        updateSliceLabel();
    end

    function updateSliceLabel()
        S = viewerFig.UserData;
        nz = nPlaneSlices(S.magVol, S.plane);
        S.lblSlice.Text = sprintf('%s %d / %d', planeLabel(S.plane), S.sliceIdx, nz);
    end

    function redraw()
        S = viewerFig.UserData;
        if S.useXlab
            cubeShow = S.cubeXVol;
            afShow = S.afXVol;
            cubeTag = 'T2xVSTD';
            afTag = 'AFxVSTD';
        else
            cubeShow = S.cubeVol;
            afShow = S.afVol;
            cubeTag = 'T2 CUBE';
            afTag = 'Anti-FLAIR';
        end

        redrawOverlay(S.ax(1), S.magVol, cubeShow, sprintf('MAG / %s', cubeTag));
        if S.hasAF
            redrawOverlay(S.ax(2), S.magVol, afShow, sprintf('MAG / %s', afTag));
            redrawOverlay(S.ax(3), cubeShow, afShow, sprintf('%s / %s', cubeTag, afTag));
        else
            redrawPlaceholder(S.ax(2), 'Anti-FLAIR not loaded');
            redrawPlaceholder(S.ax(3), 'Anti-FLAIR not loaded');
        end
        viewerFig.UserData = S;
    end

    function redrawOverlay(axh, volA, volB, titleStr)
        St = viewerFig.UserData;
        prevLim = captureAxisLim(axh);
        cla(axh);
        if isempty(volB)
            redrawPlaceholder(axh, 'Volume not available');
            return;
        end
        imA = volumePlaneSlice(volA, St.plane, St.sliceIdx);
        imB = volumePlaneSlice(volB, St.plane, St.sliceIdx);
        if ~isequal(size(imA), size(imB))
            imB = imresize(imB, size(imA));
        end
        image(axh, overlayRgb(imA, imB));
        title(axh, titleStr);
        restoreAxisLim(axh, prevLim);
    end

    function redrawPlaceholder(axh, msg)
        prevLim = captureAxisLim(axh);
        cla(axh);
        text(axh, 0.5, 0.5, msg, 'HorizontalAlignment', 'center', ...
            'Units', 'normalized', 'Color', [0.4 0.4 0.4]);
        title(axh, msg);
        restoreAxisLim(axh, prevLim);
    end

end

function sliceIdx = defaultSliceForPlane(vol, plane)
    sliceIdx = round(nPlaneSlices(vol, plane) / 2);
end

function n = nPlaneSlices(vol, plane)
    switch plane
        case 'SAGITTAL'
            n = size(vol, 3);
        case 'CORONAL'
            n = size(vol, 1);
        case 'AXIAL'
            n = size(vol, 2);
        otherwise
            n = size(vol, 3);
    end
end

function label = planeLabel(plane)
    switch plane
        case 'SAGITTAL'
            label = 'Sagittal';
        case 'CORONAL'
            label = 'Coronal';
        case 'AXIAL'
            label = 'Axial';
        otherwise
            label = plane;
    end
end

function im = volumePlaneSlice(vol, plane, sliceIdx)
    switch plane
        case 'SAGITTAL'
            im = squeeze(vol(:,:,sliceIdx))';
        case 'CORONAL'
            im = squeeze(vol(sliceIdx, :, :));
        case 'AXIAL'
            im = squeeze(vol(:, sliceIdx, :))';
        otherwise
            im = squeeze(vol(:,:,sliceIdx))';
    end
end

function vol3 = volume3D(vol)
    if ndims(vol) > 3
        vol3 = mean(vol, 4);
    else
        vol3 = vol;
    end
end

function lims = captureAxisLim(axh)
    if ~isvalid(axh) || isempty(axh.Children)
        lims = [];
        return;
    end
    xl = xlim(axh);
    yl = ylim(axh);
    if isequal(xl, [0 1]) && isequal(yl, [0 1])
        lims = [];
    else
        lims = [xl, yl];
    end
end

function restoreAxisLim(axh, lims)
    if isempty(lims)
        axis(axh, 'image');
    else
        xlim(axh, lims(1:2));
        ylim(axh, lims(3:4));
    end
    axh.XTick = [];
    axh.YTick = [];
end

function rgb = overlayRgb(imA, imB)
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
