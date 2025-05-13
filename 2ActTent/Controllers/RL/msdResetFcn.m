function in = msdResetFcn(in)
    % Random initial position and velocity
    theta0 = 0;     % Initial angle
    thetad0 = 0;    % Initial angular velocity
    x0 = 0.5;       % Initial position x
    y0 = 0;         % Initial position y

    % Assign them to Simulink model variables
    in = setVariable(in, 'theta0', theta0);
    in = setVariable(in, 'thetad0', thetad0);
    in = setVariable(in, 'x0', x0);
    in = setVariable(in, 'y0', y0);
end