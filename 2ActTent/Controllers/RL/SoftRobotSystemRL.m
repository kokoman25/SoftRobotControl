classdef SoftRobotSystemRL < matlab.System
    % Fast Restart–compliant soft robot simulation block

    properties (Access = private)
        Robot   % SoRoSim model object
        dt = 0.01
    end

    methods (Access = protected)
        function setupImpl(obj)
            % Initialize model only once
            obj.Robot = getTentacleModel();
        end

        function [Angle, AngularVelocity, sigcoordinates] = stepImpl(obj, u1, u2, ang, avel)
            % Combine current state
            state = [ang(:); avel(:)];

            % Simulate dynamics using preloaded model
            [~, qqd] = obj.Robot.dynamics(state, ...
                @(t) deal([u1; u2], [], []), ...
                dt = obj.dt, t_start = 0, t_end = 2 * obj.dt, ...
                Integrator = 'ode45');

            % Average of 3 steps
            state_out = mean(qqd(1:3,:), 1);
            Angle = state_out(1:length(ang))';
            AngularVelocity = state_out(length(ang)+1:end)';

            % Forward kinematics
            point_co = obj.Robot.FwdKinematics(Angle);
            nsig = 8;  % Hardcoded; match your robot
            Trajec = zeros(3 * nsig, 1);

            for j = 1:nsig
                Trajec(3*j-2:3*j) = point_co(4*j-3:4*j-1, 4);
            end

            sigcoordinates = Trajec;
        end

        function resetImpl(obj)
            % Fast Restart–safe: no reinitialization
        end

        function s = saveObjectImpl(obj)
            s = saveObjectImpl@matlab.System(obj);
            s.Robot = obj.Robot;
        end

        function obj = loadObjectImpl(obj, s, wasLocked)
            obj = loadObjectImpl@matlab.System(obj, s, wasLocked);
            obj.Robot = s.Robot;
        end

        function flag = isInputSizeMutableImpl(~,~)
            flag = false;
        end

        function [posOut, velOut, cooOut] = getOutputSizeImpl(~)
            posOut = [4, 1];
            velOut = [4, 1];
            cooOut = [3*8, 1];
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

function model = getTentacleModel()
    % Persistent cache of tentacle model
    persistent cachedModel
    if isempty(cachedModel)
        data = load('TentacleFiles.mat');
        cachedModel = data.Tentacle;
    end
    model = cachedModel;
end
