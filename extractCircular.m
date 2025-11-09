
function bseg_circular = extractCircular(bseg_interp, patch_interp, thresh)

stats = regionprops(bseg_interp, 'Centroid');
center = stats.Centroid;  % Get the center (Centroid) of the region

[X, Y] = ndgrid(1:size(bseg_interp,1), 1:size(bseg_interp,2));
distance_map = sqrt((X - center(2)).^2 + (Y - center(1)).^2);

intensity_norm = patch_interp / max(patch_interp(:));
distance_norm = distance_map ./ max(distance_map(:));

DISTWEIGHT = 0.2; % TEMP for the CA 
combined_score = intensity_norm - DISTWEIGHT*distance_norm;

bseg_circular = bseg_interp & (combined_score > (thresh / 100));

rnvox = sum(bseg_interp(:)) - sum(bseg_circular(:));
% if rnvox > 0
%     disp(['Voxels removed due to circularity constraint: ' num2str(rnvox) '/' num2str(sum(bseg_interp(:)))])
% end

end
