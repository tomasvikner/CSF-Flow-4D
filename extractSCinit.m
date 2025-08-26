function inds = extractSCinit(D)
    zs = 11; ze = 70; % edge exclusion (temp, or do 11 to 70)
    ymask = zeros(size(D));
    zmask = zeros(size(D));
    zmask(:, :, zs:ze) = 1;
    SCinit = (D > 0) .* zmask;
    CSFsum = squeeze(sum(sum(SCinit, 1), 3));  % sum across x and z
    SCy = find(CSFsum > 0);
    ymask(:, SCy(end-50:end), :) = 1;
    SCinit = SCinit .* ymask;
    inds = find(SCinit);
end