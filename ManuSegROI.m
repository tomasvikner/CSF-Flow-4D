function bseg = ManuSegROI(patch_interp)
% Popup to draw a manual ROI on the current interpolated patch (ax(1,3) orientation).

    bseg = [];

    regFig = uifigure('Name', 'Manu-ROI', ...
        'Position', [350 200 520 520]);
    ax = uiaxes(regFig, 'Position', [30 60 460 420]);
    colormap(ax, gray);

    imagesc(ax, flip(patch_interp', 1));
    axis(ax, 'image');
    title(ax, 'Draw ROI (double-click to finish)');
    drawnow;

    hold(ax, 'on');
    h = drawfreehand(ax);
    if isempty(h)
        hold(ax, 'off');
        delete(regFig);
        return;
    end
    wait(h);
    if ~isvalid(h)
        hold(ax, 'off');
        delete(regFig);
        return;
    end
    maskShow = createMask(h);
    delete(h);
    hold(ax, 'off');

    bseg = logical(flip(maskShow, 1)');
    if ~isequal(size(bseg), size(patch_interp))
        bseg = imresize(bseg, size(patch_interp), 'nearest');
    end
    if ~any(bseg(:))
        warning('ManuSegROI:EmptyROI', 'Manual ROI is empty.');
        bseg = [];
    end

    delete(regFig);
end
