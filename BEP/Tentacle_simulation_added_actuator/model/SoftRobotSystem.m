classdef SoftRobotSystem < matlab.System
    % SoftRobotSystem: System block to simulate one time step of a soft robot model
    % using an object loaded from matlab.mat (must contain 'T1')
    % Manually change nr of Dof of model (determine ndof by typing in command window (Linkagename).ndof)
    % Manually change nr of nsig of model (determine nsig by typing in
    %   command window (linkagename).nsig)

    properties (Access = private)
        Robot               % SoRoSim model object (T1)
        ModelLoaded = false % Makes sure the model is initialized for the simulation
        dt = 0.01           % Fixed timestep
    end

    methods (Access = protected)
        function setupImpl(obj)
            % Only load the model once
            if ~obj.ModelLoaded
                s = load('model\Tentacle1_model.mat');  % must contain object T1
                obj.Robot = s.Tentacle1;                % Store Soft body under obj.Robot
                obj.ModelLoaded = true;                 % Model is initialized and this piece of code will be skipped further on in the simulation
            end
        end

        function [Angle, AngularVelocity, sigcoordinates] = stepImpl(obj, u1, u2, ang, avel)
            % u: control input (column vector)
            % ang: current angles
            % avel: current angluar velocity

            % Combine input state
            state = [ang(:); avel(:)];

            % Make sure contol inputs only cause tentsion


            % Simulate from current state using one timestep
            [~, qqd] = obj.Robot.dynamics(state, ...
                @(t) deal([u1;u2], [], []), ...
                dt = obj.dt, t_start = 0, t_end = 2 * obj.dt, ...
                Integrator = 'ode45');

            % Return middle step result 
            state_out = (qqd(1, :) + qqd(2,:) + qqd(3,:))/3;    % Average of 3 output values
            Angle = state_out(1:length(ang))';                  % Angles of current state
            AngularVelocity = state_out(length(ang)+1:end)';    % Angular Velocities of current state

            % Return x, y, z coordinates of Tentacle
            point_co = obj.Robot.FwdKinematics(Angle);          % Optain Forward kinematics for angle (4 * 4) * nsig
            Trajec = zeros(3*8,1);                              % 3 coordinates, 8 significant points
            for j = 1:(length(point_co)/4)
                point_co(4*j-3:4*j-1,4);                        % Tranformation matrix [[R p]; [0, 1]] -> R = Rotation matrix, p = position vector
                Trajec(3*j-2:3*j,1) = point_co(4*j-3:4*j-1,4);  % Take position vector and store in column [x1;y1;z1;x2;y2;z2;enz] -> 3*8x1 vector
            end
            sigcoordinates = Trajec;
        end

        function resetImpl(obj)
            % No reset needed
        end

        % Function below are referenced during initialization to let
        % simulink know what Simulink can expect in terms of data
        % structures and they thus do not know any of the model parameters

        function flag = isInputSizeMutableImpl(~, ~)
            flag = false;
        end

        function [posOut, velOut, cooOut] = getOutputSizeImpl(obj)
            % You must match these to your model's DoFs
            N = 4;%obj.ndof;%<- change this to your number of position variables
            nsig = 8;
            posOut = [N, 1];
            velOut = [N, 1];
            cooOut = [nsig*3, 1];
        end

        function [posOut, velOut, cooOut] = getOutputDataTypeImpl(~)
            posOut = 'double';
            velOut = 'double';
            cooOut = 'double';
        end

        function [posOut, velOut, cooOut] = isOutputComplexImpl(~)
            posOut = false;
            velOut = false;
            cooOut = false;
        end

        function [posOut, velOut, cooOut] = isOutputFixedSizeImpl(~)
            posOut = true;
            velOut = true;
            cooOut = true;
        end
    end
end
