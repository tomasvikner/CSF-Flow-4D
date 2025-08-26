function centralMask = extractCentral(bwMask)
    %GETCENTRAL Extract the connected component closest to the image center.
    %   centralMask = GETCENTRAL(bwMask) returns a binary mask containing only the
    %   connected component whose centroid is closest to the center of the input mask.
    %
    %   Input:
    %       bwMask - logical binary mask (2D)
    %
    %   Output:
    %       centralMask - logical binary mask with only the central connected component

    % Get image center
    centerX = size(bwMask, 2) / 2;
    centerY = size(bwMask, 1) / 2;

    % Label connected components
    CC = bwconncomp(bwMask);

    if CC.NumObjects <= 1
        centralMask = bwMask; % only one component or none, return as is
        return;
    end

    % Compute centroids of all connected components
    stats = regionprops(CC, 'Centroid');
    centroids = cat(1, stats.Centroid);

    % Calculate distance from each centroid to image center
    dists = sqrt((centroids(:,1) - centerX).^2 + (centroids(:,2) - centerY).^2);

    % Find the component closest to the center
    [~, idxClosest] = min(dists);

    % Create output mask with only the closest component
    centralMask = false(size(bwMask));
    centralMask(CC.PixelIdxList{idxClosest}) = true;
end