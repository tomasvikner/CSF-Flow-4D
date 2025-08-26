function [flow, patch, patch_interp, bseg_interp, projV_interp, proj2D] = extractThroughPlaneFlow_interp2(data, center, direction, patch_width, SEGMODE, thresh)
    % Extracts through-plane flow using 2x interpolated patch and velocities
    % Fully vectorized version for velocity extraction and projection

    % Ensure column vectors
    center = center(:);
    direction = direction(:);

    % Image dimensions
    [nx, ny, nz] = size(data.rms);
    n_frames = size(data.vx, 4);

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

    % Vectorized velocity extraction
    N = patch_width^2;
    idx = sub2ind([nx ny nz], sample_points(1,:), sample_points(2,:), sample_points(3,:)); % 1 x N

    % Extract velocities for all frames at once (N x n_frames)
    vx_vals = zeros(N, n_frames);
    vy_vals = zeros(N, n_frames);
    vz_vals = zeros(N, n_frames);
    
    for t = 1:n_frames
        idx_t = sub2ind([nx ny nz], sample_points(1,:), sample_points(2,:), sample_points(3,:));
        vx_vals(:, t) = data.vx(idx_t + (t-1)*nx*ny*nz); % Linear index offset for frame t
        vy_vals(:, t) = data.vy(idx_t + (t-1)*nx*ny*nz);
        vz_vals(:, t) = data.vz(idx_t + (t-1)*nx*ny*nz);
    end

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
