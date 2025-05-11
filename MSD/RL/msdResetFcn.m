function in = msdResetFcn(in)
    % Random initial position and velocity
    x0 = 0;   % Random initial position
    v0 = 0;   % Random initial velocity

    % Assign them to Simulink model variables
    in = setVariable(in, 'x0', x0);
    in = setVariable(in, 'v0', v0);
end