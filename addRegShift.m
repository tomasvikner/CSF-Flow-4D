function data = addRegShift(data, movingField, delta)
    if ~isfield(data, 'shift')
        data = initRegShifts(data);
    end
    delta = round(double(delta(:))');
    if strcmp(movingField, 'T2CUBE')
        data.shift.T2CUBE = data.shift.T2CUBE + delta;
    elseif strcmp(movingField, 'AF')
        data.shift.AF = data.shift.AF + delta;
    end
end
