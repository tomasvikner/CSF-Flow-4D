% TODO, apply the ndils and nrodes
function [PC1, coords] = extractPCAFromVols(data, x, y, z, ndils, nrodes)
PC1 = [];

SE = strel('sphere', ndils);
cvol = zeros(size(data.mag));
cvol(x, y, z) = 1;
bvol = data.dist > 0;
cvol = imdilate(cvol, SE);
cvol = cvol .* bvol;
global_inds = find(cvol);

% Probably complicating things 
KMEANS = false; 
if KMEANS
    icv.cube = data.cube(global_inds); %#ok<*UNRCH> 
    icv.mag  = data.mag(global_inds);
    icv.rms  = data.rms(global_inds); 
    X = [icv.cube(:), icv.mag(:), icv.rms(:)];
    X = zscore(X);
    X(:,2) = -X(:,2);  % High cube, low mag, high rms
    [IDX, C] = kmeans(X, 2);
    [~, bestCluster] = max(sum(C, 2));
    selected_local_inds = find(IDX == bestCluster);
    selected_global_inds = global_inds(selected_local_inds); %#ok<FNDSB>
    [xx, yy, zz] = ind2sub(size(data.mag), selected_global_inds);
else % Simply just sample within the distance mask
    for i = 1:nrodes
        evol = imerode(cvol, strel('sphere', 1));
        if sum(evol(:)) > 1
            cvol = evol; 
        end
    end
    [xx, yy, zz] = ind2sub(size(data.mag), global_inds);
end

n_vox = numel(xx);
n_frames = size(data.vx, 4);

% if n_vox < 2
%     return;
% end

cvx_t = zeros(n_vox, n_frames);
cvy_t = zeros(n_vox, n_frames);
cvz_t = zeros(n_vox, n_frames);

for i = 1:n_vox
    cvx_t(i,:) = squeeze(data.vx(xx(i), yy(i), zz(i), :));
    cvy_t(i,:) = squeeze(data.vy(xx(i), yy(i), zz(i), :));
    cvz_t(i,:) = squeeze(data.vz(xx(i), yy(i), zz(i), :));
end

CV = [cvx_t', cvy_t', cvz_t'];  % frames x (n_vox * 3)
[~, SCORE, ~] = pca(CV, 'NumComponents', 1);
PC1 = SCORE(:,1);
coords = [];
coords.xx = xx;
coords.yy = yy;
coords.zz = zz;
coords.seg2D = cvol(:, :, z)';
coords.cvol = cvol;
end