% Not even used atm 
function PC1 = extractPCAFromRMS(data, x, y, z)
PC1 = [];

% Define local ROI as 5 mm radius sphere
SE = strel('sphere', 2);  % ~5mm if voxel size ~1.2–1.5 mm
cvol = zeros(size(data.rms));
cvol(x, y, z) = 1;
cvol = imdilate(cvol, SE);
global_inds = find(cvol);

% Extract local RMS values
local_rms = data.rms(global_inds);
[nx, ny] = size(local_rms);
reduced_rms = zeros(nx, ny);
reduced_rms(nx/2+1:end-nx/2, nx/2+1:end-nx/2) = 1;
reduced_rms = reduced_rms .* local_rms;
rms_thresh = 0.5 * max(reduced_rms);

% Keep only voxels with RMS above threshold
selected_global_inds = global_inds(local_rms > rms_thresh);
if isempty(selected_global_inds), return; end

[xx, yy, zz] = ind2sub(size(data.rms), selected_global_inds);
n_vox = numel(xx);
n_frames = size(data.vx, 4);
if n_vox < 2, return; end

% Extract velocity time series
cvx_t = zeros(n_vox, n_frames);
cvy_t = zeros(n_vox, n_frames);
cvz_t = zeros(n_vox, n_frames);
for i = 1:n_vox
    cvx_t(i,:) = squeeze(data.vx(xx(i), yy(i), zz(i), :));
    cvy_t(i,:) = squeeze(data.vy(xx(i), yy(i), zz(i), :));
    cvz_t(i,:) = squeeze(data.vz(xx(i), yy(i), zz(i), :));
end

% PCA on [frames x (3 * n_vox)] matrix
CV = [cvx_t', cvy_t', cvz_t'];
[~, SCORE, ~] = pca(CV, 'NumComponents', 1);
PC1 = SCORE(:,1);
end