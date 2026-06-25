function fpath = findNiiPath(fpath)
%FINDNIIPATH Return path to .nii or .nii.gz if either exists, else ''.

    if isempty(fpath)
        return;
    end
    if exist(fpath, 'file')
        return;
    end
    if endsWith(fpath, '.gz', 'IgnoreCase', true)
        alt = fpath(1:end-3);
        if exist(alt, 'file')
            fpath = alt;
        else
            fpath = '';
        end
    else
        gzpath = [fpath '.gz'];
        if exist(gzpath, 'file')
            fpath = gzpath;
        else
            fpath = '';
        end
    end
end
