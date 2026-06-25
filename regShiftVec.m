function v = regShiftVec(data, movingField)
    if ~isfield(data, 'shift')
        v = [0, 0, 0];
        return;
    end
    if strcmp(movingField, 'T2CUBE')
        v = data.shift.T2CUBE;
    elseif strcmp(movingField, 'AF')
        v = data.shift.AF;
    else
        v = [0, 0, 0];
    end
    v = double(v(:))';
end
