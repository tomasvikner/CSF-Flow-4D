% 
function [flow, cube_patch, patch_interp, bseg_interp, projV_interp, proj2D] ... 
    = extractThroughPlaneFlow_V3D(data, center, direction, patch_width, SEGMODE, thresh, DIRMODE, SHAPE, CLIPON, CLIPVAL, DILATE, bsegManual, LOCALCS)
    % Extracts through-plane flow using 2x interpolated patch and velocities
    % Fully vectorized version for velocity extraction and projection

    if nargin < 13 || isempty(LOCALCS)
        LOCALCS = 'FULL';
    end
    if nargin < 12
        bsegManual = [];
    end
    if nargin < 11 || isempty(DILATE)
        DILATE = 0;
    end
    if nargin < 10 || isempty(CLIPVAL)
        CLIPVAL = 99;
    end
    if nargin < 9 || isempty(CLIPON)
        CLIPON = true;
    end

    % Ensure column vectors
    center = center(:); % center in AP, SI, LR 
    direction = direction.(DIRMODE)(:);

    % Image dimensions
    [nx, ny, nz] = size(data.vstd);
    if isfield(data, 'hasVel') && data.hasVel
        n_frames = size(data.vx, 2);
    else
        n_frames = 0;
    end

    % *** From arrow direction to direction in data ***
    % --- GUI-arrow direction is LR, AP, SI
    % --- Data direction is AP, SI, LR
    n = zeros(3, 1);
    n(1) = direction(2); % GUI AP to data AP 
    n(2) = direction(3); % GUI SI to data SI 
    n(3) = direction(1); % GUI LR to data LR 

    % Now normal should be in data-direction, not arrow-direction
    n = n/norm(n);
    a = [1; 0; 0];
    u = (a - dot(a,n)*n);  u = u/norm(u);
    v = cross(n,u);        v = v/norm(v);
    % u, v should span the local 2D-plane

    % Create regular grid in-plane
    pad_width = (patch_width - 1)/2;
    [U, V] = meshgrid(-pad_width:pad_width, -pad_width:pad_width);
    spoints = center + u * U(:)' + v * V(:)';  % [3 x N]
    spoints = round(spoints);
    % Sample points (spoints) should contain the x, y, z coordinates in the
    % data.(SEGMODE) coordinates (nx, ny, nz)

    % Clip to valid indices
    spoints(1,:) = min(max(spoints(1,:), 1), nx);
    spoints(2,:) = min(max(spoints(2,:), 1), ny);
    spoints(3,:) = min(max(spoints(3,:), 1), nz);

    % Linear index for data patch extraction
    lin_idx = sub2ind([nx ny nz], spoints(1,:), spoints(2,:), spoints(3,:));
    patch_vals = data.(SEGMODE)(lin_idx);
    patch = reshape(patch_vals, patch_width, patch_width);

    % Interpolate patch to 2× resolution
    [y0, x0] = size(patch);
    xi = linspace(1, x0, 2*x0);
    yi = linspace(1, y0, 2*y0);
    [XI, YI] = meshgrid(xi, yi);
    patch_interp = interp2(patch, XI, YI, 'cubic');

    cube_patch_vals = data.CUBE(lin_idx);
    cube_patch = reshape(cube_patch_vals, patch_width, patch_width);
    cube_patch = interp2(cube_patch, XI, YI, 'cubic');
    
    if CLIPON
        clip_val = prctile(patch_interp(:), CLIPVAL);
        patch_interp(patch_interp > clip_val) = clip_val;
    end

    if ~isempty(bsegManual)
        bseg_interp = logical(bsegManual);
        if ~isequal(size(bseg_interp), size(patch_interp))
            error('extractThroughPlaneFlow_V3D:BadManualMask', ...
                'Manual mask size must match interpolated patch.');
        end
    else
        % Segment interpolated patch; threshold max from central disk only
        patch_core = patchValsInCenterDisk(patch_interp, patch_width);
        mval = max(patch_core(:), [], 'omitnan');
        if isempty(mval) || ~isfinite(mval)
            mval = max(patch_interp(:), [], 'omitnan');
        end
        bseg_interp = patch_interp > (thresh/100) * mval;
        if patch_width < 50 % avoid applying this for SC 
            bseg_interp = imfill(bseg_interp, 'holes'); 
            bseg_interp = bwareaopen(bseg_interp, 10);
        end
        bseg_interp = applyLocalCsSeg(logical(bseg_interp), LOCALCS);

        if strcmp(SHAPE, 'CIRC')
            bseg_circular = extractCircular(bseg_interp, patch_interp, thresh);
            bseg_interp = applyLocalCsSeg(bseg_circular, LOCALCS);
        elseif strcmp(SHAPE, 'RECT')
            bseg_rectangular = extractRectangular(bseg_interp, patch_interp, thresh);
            bseg_interp = applyLocalCsSeg(bseg_rectangular, LOCALCS);
        end

        if DILATE > 0
            SE1 = strel('sphere', DILATE);
            bseg_interp = imdilate(bseg_interp, SE1);
        end
    end

    hasVel = isfield(data, 'hasVel') && data.hasVel;
    if ~hasVel
        flow = [];
        proj2D = struct('proj', [], 'velx', [], 'vely', [], 'velz', []);
        return;
    end

    % Vectorized velocity extraction
    N = patch_width^2;

    % Extract velocities for all frames at once (N x n_frames)
    vx_vals = zeros(N, n_frames);
    vy_vals = zeros(N, n_frames);
    vz_vals = zeros(N, n_frames);

    % New sample points mapping from global-roi-patch index 
    idx_g = sub2ind([nx ny nz], spoints(1,:), spoints(2,:), spoints(3,:));
    idx_r = data.imap(idx_g); % This maps global inds to local inds if MASKED loading
    deli = idx_r==0; % remove these 
    idx_p = find(idx_r ~= 0);
    idx_r(deli) = [];
    vx_vals(idx_p, :) = data.vx(idx_r, :);
    vy_vals(idx_p, :) = data.vy(idx_r, :);
    vz_vals(idx_p, :) = data.vz(idx_r, :);

    % Reshape to patch_width x patch_width x n_frames
    vx_frame = reshape(vx_vals, patch_width, patch_width, n_frames);
    vy_frame = reshape(vy_vals, patch_width, patch_width, n_frames);
    vz_frame = reshape(vz_vals, patch_width, patch_width, n_frames);

    % Interpolate velocity fields for each frame
    [ny2, nx2] = size(patch_interp);
    vx_plane = zeros(ny2, nx2, n_frames);
    vy_plane = zeros(ny2, nx2, n_frames);
    vz_plane = zeros(ny2, nx2, n_frames);
    for t = 1:n_frames
        vx_plane(:,:,t) = interp2(vx_frame(:,:,t), XI, YI, 'cubic');
        vy_plane(:,:,t) = interp2(vy_frame(:,:,t), XI, YI, 'cubic');
        vz_plane(:,:,t) = interp2(vz_frame(:,:,t), XI, YI, 'cubic');
    end

    % Vectorized velocity projection onto normal direction
    vx_vec = reshape(vx_plane, [], n_frames);
    vy_vec = reshape(vy_plane, [], n_frames);
    vz_vec = reshape(vz_plane, [], n_frames);
    normal = direction / norm(direction);
    projV_vec = normal(1)*vx_vec + normal(2)*vy_vec + normal(3)*vz_vec;
    projV_interp = reshape(projV_vec, ny2, nx2, n_frames);

    % Compute final flow by summing over mask pixels for each frame
    [yy, xx] = find(bseg_interp);
    lin_mask_idx = sub2ind([ny2, nx2], yy, xx);
    projV_masked = projV_vec(lin_mask_idx, :);
    flow = sum(projV_masked, 1, 'omitnan') * 0.0625 * 1e-03; % mm3/s to ml/s 

    % Return struct of velocity fields
    proj2D = struct();
    proj2D.proj = projV_interp;
    proj2D.velx = vx_plane;
    proj2D.vely = vy_plane;
    proj2D.velz = vz_plane;

end

function mask = applyLocalCsSeg(mask, localCsMode)
    if strcmp(localCsMode, 'FULL')
        mask = logical(mask);
    else
        mask = extractSegment(logical(mask), localCsMode);
    end
end

function vals = patchValsInCenterDisk(patch, patchWidth)
    % Central disk for local threshold: radius = 40% of patch width
    % (diameter ~80% of patch), centered on the interpolated plane.
    [ny, nx] = size(patch);
    cx = (nx + 1) / 2;
    cy = (ny + 1) / 2;
    radius = 0.4 * patchWidth * (nx / patchWidth);
    [X, Y] = meshgrid(1:nx, 1:ny);
    disk = (X - cx).^2 + (Y - cy).^2 <= radius^2;
    vals = patch(disk);
end
