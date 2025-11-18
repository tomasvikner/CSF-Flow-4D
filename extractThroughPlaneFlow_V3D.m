% Todo need to update this sampling a bit 
function [flow, cube_patch, patch_interp, bseg_interp, projV_interp, proj2D] ... 
    = extractThroughPlaneFlow_V3D(data, center, direction, patch_width, SEGMODE, thresh, DIRMODE, SHAPE)
    % Extracts through-plane flow using 2x interpolated patch and velocities
    % Fully vectorized version for velocity extraction and projection

    % Ensure column vectors
    center = center(:);

    direction = direction.(DIRMODE)(:);
    % disp(DIRMODE)

    % Image dimensions
    [nx, ny, nz] = size(data.rms);
    n_frames = size(data.vx, 2);

    % Define local plane basis (in x-z plane)
    dir_sag = [direction(1); direction(3)];
    dir_sag(2) = -dir_sag(2);  % flip z axis for image convention
    dir_sag = dir_sag / norm(dir_sag);

    u = [dir_sag(1); 0; dir_sag(2)];
    v = [-dir_sag(2); 0; dir_sag(1)];

    pad_width = (patch_width - 1)/2;

    % Create regular grid in-plane
    [U, V] = meshgrid(-pad_width:pad_width, -pad_width:pad_width);
    sample_points = center + u * U(:)' + v * V(:)';  % [3 x N]
    sample_points = round(sample_points);

    % Clip to valid indices
    sample_points(1,:) = min(max(sample_points(1,:), 1), nx);
    sample_points(2,:) = min(max(sample_points(2,:), 1), ny);
    sample_points(3,:) = min(max(sample_points(3,:), 1), nz);

    % Linear index for data patch extraction
    lin_idx = sub2ind([nx ny nz], sample_points(1,:), sample_points(2,:), sample_points(3,:));
    patch_vals = data.(SEGMODE)(lin_idx);
    patch = reshape(patch_vals, patch_width, patch_width);

    % Interpolate patch to 2× resolution
    [y0, x0] = size(patch);
    [x, y] = meshgrid(1:x0, 1:y0);
    xi = linspace(1, x0, 2*x0);
    yi = linspace(1, y0, 2*y0);
    [XI, YI] = meshgrid(xi, yi);
    patch_interp = interp2(patch, XI, YI, 'cubic');

    cube_patch_vals = data.cube(lin_idx);
    cube_patch = reshape(cube_patch_vals, patch_width, patch_width);
    cube_patch = interp2(cube_patch, XI, YI, 'cubic');
    
    % Clip at 99
    clip_val = prctile(patch_interp(:), 99);
    patch_interp(patch_interp > clip_val) = clip_val;

    % Segment interpolated patch
    mval = max(patch_interp(:), [], 'omitnan');
    bseg_interp = patch_interp > (thresh/100) * mval;
    if patch_width < 50 % avoid applying this for SC 
        bseg_interp = imfill(bseg_interp, 'holes'); 
        bseg_interp = bwareaopen(bseg_interp, 10);
    end
    bseg_interp = extractCentral(logical(bseg_interp));

    % TEMP: using circularity constraint
    if strcmp(SHAPE, 'CIRC')
        bseg_circular = extractCircular(bseg_interp, patch_interp, thresh);
        bseg_interp = extractCentral(bseg_circular);
    elseif strcmp(SHAPE, 'RECT')
        bseg_rectangular = extractRectangular(bseg_interp, patch_interp, thresh);
        bseg_interp = extractCentral(bseg_rectangular);
    end

    % TEMP: dilate 2 voxels for CA
    SE1 = strel('sphere', 2);
    bseg_interp = imdilate(bseg_interp, SE1);

    % Vectorized velocity extraction
    N = patch_width^2;
    idx = sub2ind([nx ny nz], sample_points(1,:), sample_points(2,:), sample_points(3,:)); % 1 x N

    % Extract velocities for all frames at once (N x n_frames)
    vx_vals = zeros(N, n_frames);
    vy_vals = zeros(N, n_frames);
    vz_vals = zeros(N, n_frames);

    % New sample points mapping from global-roi-patch index 
    idx_g = sub2ind([nx ny nz], sample_points(1,:), sample_points(2,:), sample_points(3,:));
    idx_r = data.imap(idx_g);
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
