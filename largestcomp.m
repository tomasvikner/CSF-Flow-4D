function vol = largestcomp(vol)
    if sum(vol) == 0
        return;
    end
    CC = bwconncomp(vol);
    nl = zeros(CC.NumObjects, 1);
    for k = 1:CC.NumObjects
        nl(k) = numel(CC.PixelIdxList{k});
    end
    [~, mi] = max(nl);
    vol = zeros(size(vol));
    vol(CC.PixelIdxList{mi}) = 1;
end