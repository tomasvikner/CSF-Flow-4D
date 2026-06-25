function d = tform2dToVolShift(tform2d)
    % Display-plane tx, ty -> volume [AP, LR, SI] on vol(:,:,slice)'.
    d = [round(tform2d.T(3, 1)), round(tform2d.T(3, 2)), 0];
end
