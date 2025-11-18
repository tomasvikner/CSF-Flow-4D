function bseg_rect = extractRectangular(bseg_interp, patch_interp, thresh)

% --- Find centroid (same as circular version) ---
stats = regionprops(bseg_interp, 'Centroid');
center = stats.Centroid;   % [x, y]

% --- Build horizontal-only distance map ---
[X, Y] = ndgrid(1:size(bseg_interp,1), 1:size(bseg_interp,2));

% Horizontal distance from vertical midline (centroid.x)
% NOTE: center(1) = x-coordinate in MATLAB’s regionprops
horizontal_distance = abs(Y - center(1));

% Normalize
horizontal_norm = horizontal_distance ./ max(horizontal_distance(:));

% --- Normalize intensity patch ---
intensity_norm = patch_interp / max(patch_interp(:));

% --- Weighted combined score ---
DISTWEIGHT = 0.4;   % TEMP for ACA A2 SAS
combined_score = intensity_norm - DISTWEIGHT * horizontal_norm;

% --- Apply threshold ---
bseg_rect = bseg_interp & (combined_score > (thresh / 100));

% Optional: removed voxels count
% rnvox = sum(bseg_interp(:)) - sum(bseg_rect(:));
% disp(['Removed due to rectangular constraint: ' num2str(rnvox)]);
end
