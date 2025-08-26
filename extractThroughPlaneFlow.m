function [flow, patch, patch_interp, bseg, bseg_interp, projV, proj2D] = extractThroughPlaneFlow(data, center, direction, sz, MODE)
        % Assumes sample_points are valid and data.rms is same size as velocity frame

        % Ensure column vectors
        center = center(:);
        direction = direction(:);
    
        % Extract size info
        [nx, ny, nz] = size(data.rms);
        n_frames = size(data.vx, 4);
        full_sz = 2*sz + 1;
    
        % Build plane basis (in x-z plane)
        dir_sag = [direction(1); direction(3)];
        dir_sag(2) = -dir_sag(2);  % Flip z for image Y
        dir_sag = dir_sag / norm(dir_sag);
    
        u = [dir_sag(1); 0; dir_sag(2)];  % in-plane X direction (x-z)
        v = [-dir_sag(2); 0; dir_sag(1)]; % orthogonal in-plane Y
    
        % Build 2D grid
        [U, V] = meshgrid(-sz:sz, -sz:sz);
        sample_points = center + u * U(:)' + v * V(:)';  % [3 x N]
        sample_points = round(sample_points);
    
        % Clip to valid range
        sample_points(1,:) = min(max(sample_points(1,:), 1), nx);
        sample_points(2,:) = min(max(sample_points(2,:), 1), ny);
        sample_points(3,:) = min(max(sample_points(3,:), 1), nz);

        % Linear indices
        lin_idx = sub2ind([nx ny nz], sample_points(1,:), sample_points(2,:), sample_points(3,:));

        % Extract patch
        patch_vals = data.(MODE)(lin_idx);
        patch = reshape(patch_vals, full_sz, full_sz);

        [ny, nx] = size(patch);
        [x, y] = meshgrid(1:nx, 1:ny);
        xi = linspace(1, nx, 4*nx);
        yi = linspace(1, ny, 4*ny);
        [XI, YI] = meshgrid(xi, yi);
        patch_interp = interp2(patch, XI, YI, 'cubic');
        mval = max(patch_interp(:), [], 'omitnan');
        bseg = patch_interp > 0.5 * mval;
        bseg = imfill(bseg, 'holes');         % fill holes to get a proper region
        bseg = bwareaopen(bseg, 10);          % remove small noisy regions (adjust 10 as needed)
        bseg = logical(bseg);                  % ensure logical
        bseg_interp = extractCentral(bseg); % currently for show

        % Segment region of interest
        mval = max(patch(:), [], 'omitnan');
        bseg = patch > 0.5 * mval;
        bseg = imfill(bseg, 'holes'); 
        bseg = bwareaopen(bseg, 10);
        bseg = extractCentral(bseg); 
        seg_points_3D = sample_points(:, bseg);
        n_pts = size(seg_points_3D, 2);

        % Project flow
        projV = zeros(n_pts, n_frames);
        normal = direction / norm(direction);

        % if we just the interpolated bseg, gotta interp velocities too
        for i = 1:n_pts
            xi = seg_points_3D(1, i);
            yi = seg_points_3D(2, i);
            zi = seg_points_3D(3, i);
            vx_t = squeeze(data.vx(xi, yi, zi, :));
            vy_t = squeeze(data.vy(xi, yi, zi, :));
            vz_t = squeeze(data.vz(xi, yi, zi, :));
            v3d = [vx_t, vy_t, vz_t];  % [frames x 3]
            projV(i, :) = v3d * normal(:);
        end

        flow = sum(projV, 1, 'omitnan');

        % Get 2D plane 
        all_inds_2D = 1:numel(bseg);
        all_points_3D = sample_points(:, all_inds_2D);
        na_pts = size(all_points_3D, 2);
        projP = zeros(na_pts, n_frames);
        vx = zeros(na_pts, n_frames);
        vy = zeros(na_pts, n_frames);
        vz = zeros(na_pts, n_frames);
        for i = 1:na_pts
            xi = all_points_3D(1, i);
            yi = all_points_3D(2, i);
            zi = all_points_3D(3, i);
            vx_t = squeeze(data.vx(xi, yi, zi, :));
            vy_t = squeeze(data.vy(xi, yi, zi, :));
            vz_t = squeeze(data.vz(xi, yi, zi, :));
            v3d = [vx_t, vy_t, vz_t];  % [frames x 3]
            projP(i, :) = v3d * normal(:);
            vx(i, :) = vx_t;
            vy(i, :) = vy_t;
            vz(i, :) = vz_t;
        end
        nx = sqrt(na_pts); 
        proj2D = [];
        proj2D.proj = reshape(projP, nx, nx, n_frames);
        proj2D.velx = reshape(vx, nx, nx, n_frames);
        proj2D.vely = reshape(vy, nx, nx, n_frames);
        proj2D.velz = reshape(vz, nx, nx, n_frames);

    end