% Prob no need, completed ITK coreg may be better anyway 
function coregCubeToMag(fig, data)
    if isempty(data.cube) || isempty(data.mag)
        uialert(fig, 'Cube or Mag not loaded.', 'Coregistration Error');
        return;
    end

    disp('coregCubeToMag()')

    % Prepare fixed and moving images (convert to single/double)
    fixed = data.mag;
    moving = data.cube;

    % Define the optimizer and metric for rigid registration
    [optimizer, metric] = imregconfig('multimodal');
    optimizer.MaximumIterations = 300;
    optimizer.GrowthFactor = 1.025;
    optimizer.InitialRadius = 2.5e-03;
    optimizer.Epsilon = 1.5e-6;

    % Perform rigid 3D registration
    tform = imregtform(moving, fixed, 'rigid', optimizer, metric);

    % Apply the transformation to the moving image (cube)
    Rfixed = imref3d(size(fixed));
    registeredCube = imwarp(moving, tform, 'OutputView', Rfixed);

    % Update data.cube with registered volume
    data.cube = registeredCube;

    updateDisplays();
end