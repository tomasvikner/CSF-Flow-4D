function data = coregSegmodeToMag3D(fig, data, axFixed, axMoving, sliceIdx, movingField, fixedField, slabHalf)
% 3D translation of movingField toward fixedField using zoomed in-plane ROI
% and +/- slabHalf slices around sliceIdx (default 7).

    if nargin < 7 || isempty(fixedField)
        fixedField = 'mag';
    end
    if nargin < 8 || isempty(slabHalf)
        slabHalf = 7;
    end

    if ~isequal(xlim(axFixed), xlim(axMoving)) || ~isequal(ylim(axFixed), ylim(axMoving))
        disp('AutoReg-3D: ax(2,1) and ax(2,2) zoom differ; using ax(2,1) ROI.');
    end
    if ~isfield(data, fixedField) || isempty(data.(fixedField))
        uialert(fig, [fixedField ' not loaded.'], 'AutoReg-3D');
        return;
    end
    if ~isfield(data, movingField) || isempty(data.(movingField))
        uialert(fig, [movingField ' not loaded.'], 'AutoReg-3D');
        return;
    end

    disp(sprintf('coregSegmodeToMag3D(): %s -> %s', movingField, fixedField));

    fixedVol = single(volume3D(data.(fixedField)));
    movingVol = single(volume3D(data.(movingField)));

    nz = min(size(fixedVol, 3), size(movingVol, 3));
    if sliceIdx < 1 || sliceIdx > nz
        uialert(fig, 'Slice index out of range.', 'AutoReg-3D');
        return;
    end

    fixedSl = squeeze(fixedVol(:,:,sliceIdx))';
    [yIdx, xIdx] = zoomRoi(axFixed, size(fixedSl));
    if isempty(yIdx) || isempty(xIdx)
        uialert(fig, 'Zoom in on ax(2,1) / ax(2,2) before running AutoReg-3D.', 'AutoReg-3D');
        return;
    end

    apIdx = xIdx;
    lrIdx = yIdx;
    siLo = max(1, sliceIdx - slabHalf);
    siHi = min(nz, sliceIdx + slabHalf);
    siIdx = siLo:siHi;

    fixedROI3 = fixedVol(apIdx, lrIdx, siIdx);
    movingROI3 = movingVol(apIdx, lrIdx, siIdx);

    if numel(fixedROI3) < 100
        uialert(fig, 'Zoomed region is too small for registration.', 'AutoReg-3D');
        return;
    end

    [optimizer, metric] = imregconfig('multimodal');
    optimizer.MaximumIterations = 300;
    optimizer.GrowthFactor = 1.025;
    optimizer.InitialRadius = 2.5e-3;
    optimizer.Epsilon = 1.5e-6;

    % imregtform needs each dim >= 2^(PyramidLevels+1); +/-7 slab can be 15 on SI.
    pyramidLevels = pyramidLevelsForRoi(fixedROI3);
    disp(sprintf('AutoReg-3D ROI size [%s], PyramidLevels=%d', ...
        num2str(size(fixedROI3)), pyramidLevels));

    tform3d = imregtform(movingROI3, fixedROI3, 'translation', optimizer, metric, ...
        'PyramidLevels', pyramidLevels);
    tr = tform3d.Translation;
    disp(sprintf('AutoReg-3D %s -> %s (slice %d, SI %d:%d): tx=%.3f, ty=%.3f, tz=%.3f', ...
        movingField, fixedField, sliceIdx, siLo, siHi, tr(1), tr(2), tr(3)));

    Rout = imref3d(size(fixedVol));
    Rin = imref3d(size(movingVol));
    registeredVol = imwarp(movingVol, Rin, tform3d, 'OutputView', Rout);

    data.(movingField) = cast(registeredVol, class(data.(movingField)));
    data = addRegShift(data, movingField, tform3dToVolShift(tform3d));
end

function levels = pyramidLevelsForRoi(vol)
    minDim = min(size(vol));
    levels = 1;
    for k = 3:-1:1
        if minDim >= 2^(k + 1)
            levels = k;
            return;
        end
    end
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
