classdef SoftRobotSystem < matlab.System
    % SoftRobotSystem: System block to simulate one time step of a soft robot model
    % using an object loaded from matlab.mat (must contain 'T1')

    properties (Access = private)
        Robot          % SoRoSim model object (T1)
        ModelLoaded = false
        dt = 0.01      % Fixed timestep
    end

    methods (Access = protected)
        function setupImpl(obj)
            % Only load the model once
            if ~obj.ModelLoaded
                s = load('matlab.mat');  % must contain object T1
                obj.Robot = s.T1;
                obj.ModelLoaded = true;
            end
        end

        function [position, velocity] = stepImpl(obj, u, pos, vel)
            % u: control input (column vector)
            % pos: current position
            % vel: current velocity

            % Combine input state
            state = [pos(:); vel(:)];

            % Simulate from current state using one timestep
            [~, qqd] = obj.Robot.dynamics(state', ...
                @(t) deal(50 * u(:), [], []), ...
                dt = obj.dt, t_start = 0, t_end = obj.dt, ...
                Integrator = 'ode45');

            % Return middle step result (or final)
            state_out = qqd(end, :);
            position = state_out(1:length(pos))';
            velocity = state_out(length(pos)+1:end)';
        end

        function resetImpl(obj)
            % No reset needed
        end

        function flag = isInputSizeMutableImpl(~, ~)
            flag = false;
        end

        function [posOut, velOut] = getOutputSizeImpl(obj)
            % You must match these to your model's DoFs
            N = 4; % <- change this to your number of position variables
            posOut = [N, 1];
            velOut = [N, 1];
        end

        function [posOut, velOut] = getOutputDataTypeImpl(~)
            posOut = 'double';
            velOut = 'double';
        end

        function [posOut, velOut] = isOutputComplexImpl(~)
            posOut = false;
            velOut = false;
        end

        function [posOut, velOut] = isOutputFixedSizeImpl(~)
            posOut = true;
            velOut = true;
        end
    end
end
