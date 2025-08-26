function [flow, patch, patch_interp, bseg_interp, projV_interp, proj2D] = extractThroughPlaneFlow_interp(data, center, direction, patch_width, MODE, thresh)
    % Extracts through-plane flow using 2x interpolated patch and velocities

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

    % Linear index for data
    lin_idx = sub2ind([nx ny nz], sample_points(1,:), sample_points(2,:), sample_points(3,:));
    patch_vals = data.(MODE)(lin_idx);
    patch = reshape(patch_vals, patch_width, patch_width);

    % Interpolate patch to 2× resolution
    [y0, x0] = size(patch);
    [x, y] = meshgrid(1:x0, 1:y0);
    xi = linspace(1, x0, 2*x0);
    yi = linspace(1, y0, 2*y0);
    [XI, YI] = meshgrid(xi, yi);
    patch_interp = interp2(patch, XI, YI, 'cubic');

    % Segment interpolated patch
    mval = max(patch_interp(:), [], 'omitnan');
    bseg_interp = patch_interp > (thresh/100) * mval;
    bseg_interp = imfill(bseg_interp, 'holes');
    bseg_interp = bwareaopen(bseg_interp, 10);
    bseg_interp = extractCentral(logical(bseg_interp));

    % Interpolate velocity fields
    [ny2, nx2] = size(patch_interp);
    vx_plane = zeros(ny2, nx2, n_frames);
    vy_plane = zeros(ny2, nx2, n_frames);
    vz_plane = zeros(ny2, nx2, n_frames);

    sample_points_3D = reshape(sample_points.', patch_width, patch_width, 3);
    for t = 1:n_frames
        vx_frame = zeros(patch_width, patch_width);
        vy_frame = zeros(patch_width, patch_width);
        vz_frame = zeros(patch_width, patch_width);
        for i = 1:patch_width
            for j = 1:patch_width
                xi = sample_points_3D(i,j,1);
                yi = sample_points_3D(i,j,2);
                zi = sample_points_3D(i,j,3);
                vx_frame(i,j) = data.vx(xi,yi,zi,t);
                vy_frame(i,j) = data.vy(xi,yi,zi,t);
                vz_frame(i,j) = data.vz(xi,yi,zi,t);
            end
        end
        vx_plane(:,:,t) = interp2(vx_frame, XI, YI, 'cubic');
        vy_plane(:,:,t) = interp2(vy_frame, XI, YI, 'cubic');
        vz_plane(:,:,t) = interp2(vz_frame, XI, YI, 'cubic');
    end

    % Project velocity at every pixel (not just in the mask)
    projV_interp = zeros(ny2, nx2, n_frames);
    normal = direction(:) / norm(direction); 
    for t = 1:n_frames
        for i = 1:ny2
            for j = 1:nx2
                v3d = [vx_plane(i,j,t), vy_plane(i,j,t), vz_plane(i,j,t)];
                projV_interp(i,j,t) = dot(v3d, normal);
            end
        end
    end

    % Final flow curve (summed projected velocities)
    [yy, xx] = find(bseg_interp);
    projV_vec = reshape(projV_interp, nx2*ny2, n_frames);
    flow = sum(projV_vec(sub2ind([ny2, nx2], yy, xx), :), 1, 'omitnan');

    % Return full velocity fields
    proj2D = struct();
    proj2D.proj = projV_interp;
    proj2D.velx = vx_plane;
    proj2D.vely = vy_plane;
    proj2D.velz = vz_plane;

end
