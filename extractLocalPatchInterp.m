function patch_interp = extractLocalPatchInterp(vol, center, direction, dirMode, patch_width)
% Interpolated 2x local cross-section patch (same plane as extractThroughPlaneFlow_V3D).

    center = center(:);
    direction = direction.(dirMode)(:);

    [nx, ny, nz] = size(vol);

    % GUI-arrow direction is LR, AP, SI; data direction is AP, SI, LR
    n = zeros(3, 1);
    n(1) = direction(2);
    n(2) = direction(3);
    n(3) = direction(1);
    n = n / norm(n);

    a = [1; 0; 0];
    u = (a - dot(a, n) * n);
    u = u / norm(u);
    v = cross(n, u);
    v = v / norm(v);

    pad_width = (patch_width - 1) / 2;
    [U, V] = meshgrid(-pad_width:pad_width, -pad_width:pad_width);
    spoints = center + u * U(:)' + v * V(:)';
    spoints = round(spoints);

    spoints(1, :) = min(max(spoints(1, :), 1), nx);
    spoints(2, :) = min(max(spoints(2, :), 1), ny);
    spoints(3, :) = min(max(spoints(3, :), 1), nz);

    lin_idx = sub2ind([nx, ny, nz], spoints(1, :), spoints(2, :), spoints(3, :));
    patch_vals = vol(lin_idx);
    patch = reshape(patch_vals, patch_width, patch_width);

    [y0, x0] = size(patch);
    xi = linspace(1, x0, 2 * x0);
    yi = linspace(1, y0, 2 * y0);
    [XI, YI] = meshgrid(xi, yi);
    patch_interp = interp2(patch, XI, YI, 'cubic');
end
