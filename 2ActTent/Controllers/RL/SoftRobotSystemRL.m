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
            % Fast Restart–safe: no reinitialization
        end

        function s = saveObjectImpl(obj)
            s = saveObjectImpl@SoftRobotSystemRL(obj);
        end

        function obj = loadObjectImpl(obj, s, wasLocked)
            obj = loadObjectImpl@SoftRobotSystemRL(obj, s, wasLocked);
            obj.Robot = getTentacleModel();
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
