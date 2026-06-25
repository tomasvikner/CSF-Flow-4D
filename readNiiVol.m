function vol = readNiiVol(fpath)
%READNIIVOL Load NIfTI volume (.nii or .nii.gz).
% Uses niftiread for gzip files (MRIread appends .gz twice on .nii.gz paths).

    fpath = findNiiPath(fpath);
    if isempty(fpath)
        error('readNiiVol:FileNotFound', 'NIfTI file not found.');
    end

    if endsWith(fpath, '.gz', 'IgnoreCase', true)
        vol = niftiread(fpath);
        return;
    end

    mri = MRIread(fpath);
    if isstruct(mri) && isfield(mri, 'vol')
        vol = mri.vol;
    else
        error('readNiiVol:BadMRIread', 'MRIread failed for %s', fpath);
    end
end
