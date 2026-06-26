function segMask = extractSegment(bwMask, mode)
% extractSegment  Select connected component(s) from a binary patch mask.
%   mode 'FULL'    - keep all components (no change)
%   mode 'CENTRAL' - keep component with centroid closest to patch center
%   mode 'LARGEST' - keep largest component by area

    bwMask = logical(bwMask);
    if nargin < 2 || isempty(mode)
        mode = 'CENTRAL';
    end

    if strcmp(mode, 'FULL')
        segMask = bwMask;
        return;
    end

    CC = bwconncomp(bwMask);
    if CC.NumObjects <= 1
        segMask = bwMask;
        return;
    end

    sizes = cellfun(@numel, CC.PixelIdxList);
    if strcmp(mode, 'LARGEST')
        [~, pickIdx] = max(sizes);
    elseif strcmp(mode, 'CENTRAL')
        centerX = size(bwMask, 2) / 2;
        centerY = size(bwMask, 1) / 2;
        stats = regionprops(CC, 'Centroid');
        centroids = cat(1, stats.Centroid);
        dists = sqrt((centroids(:, 1) - centerX).^2 + (centroids(:, 2) - centerY).^2);
        [~, pickIdx] = min(dists);
    else
        error('extractSegment:BadMode', 'Unknown mode: %s', mode);
    end

    segMask = false(size(bwMask));
    segMask(CC.PixelIdxList{pickIdx}) = true;
end
