function centralMask = extractCentral(bwMask)
% extractCentral  Keep one connected component from a binary patch mask.
%
% If the largest component is at least 2x the second-largest, keep largest.
% Otherwise keep the component whose centroid is closest to the patch center.

centerX = size(bwMask, 2) / 2;
centerY = size(bwMask, 1) / 2;

CC = bwconncomp(bwMask);

if CC.NumObjects <= 1
    centralMask = bwMask;
    return;
end

sizes = cellfun(@numel, CC.PixelIdxList);
[sortedSizes, sortIdx] = sort(sizes, 'descend');

if sortedSizes(1) >= 2 * sortedSizes(2)
    pickIdx = sortIdx(1);
else
    stats = regionprops(CC, 'Centroid');
    centroids = cat(1, stats.Centroid);
    dists = sqrt((centroids(:,1) - centerX).^2 + (centroids(:,2) - centerY).^2);
    [~, pickIdx] = min(dists);
end

centralMask = false(size(bwMask));
centralMask(CC.PixelIdxList{pickIdx}) = true;
end
