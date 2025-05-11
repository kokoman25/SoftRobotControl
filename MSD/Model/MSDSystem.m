classdef MSDSystem < matlab.System
    % MassSpringDamper simulates a mass-spring-damper system
    %   This system object models the dynamics of a mass-spring-damper
    %   system with parameters m (mass), c (damping), and k (spring
    %   constant). It accepts an external force input and outputs the
    %   position and velocity of the mass.

    properties
        m = 1;      % Mass (kg)
        c = 0.5;    % Damping coefficient (N*s/m)
        k = 10;     % Spring constant (N/m)
    end

    properties(Access = private)
        Ts = 0.01;  %Sample time (s)
    end

    methods (Access = protected)
        %% Common functions
        function [x_out, v_out] = stepImpl(obj, u)
            % Compute acceleration
            a = (u - obj.c * obj.v - obj.k * obj.x) / obj.m;

            %update velocity and position using Euler integration
            obj.v = obj.v + obj.Ts * a;
            obj.x = obj.x + obj.Ts * obj.v;

            % Output current position and velocity
            x_out = obj.x;
            v_out = obj.v;
        end

        function resetImpl(obj)
            % Reset states to initial conditions
            obj.x = randn();
            obj.v = randn();
        end

        function sts = getSampleTimeImpl(obj)
            % Define sample time
            sts = createSampleTime(obj, 'Type', 'Discrete', 'SampleTime', obj.Ts);
        end

        function [out1, out2] = getOutputSizeImpl(~)
            % Define output sizes
            out1 = [1, 1];
            out2 = [1, 1];
        end

        function [out1, out2] = getOutputDataTypeImpl(~)
            % Define output data types
            out1 = 'double';
            out2 = 'double';
        end

        function [out1, out2] = isOutputComplexImpl(~)
            % Define output complexity
            out1 = false;
            out2 = false;
        end

        function [out1, out2] = isOutputFixedSizeImpl(~)
            % Define output fixed size
            out1 = true;
            out2 = true;
        end

        function icon = getIconImpl(~)
            % Define icon for system block
            icon = 'Mass-Spring-Damper';
        end

        function names = getInputNamesImpl(~)
            % Define input port names
            names = {'u'};
        end

        function names = getOutputNamesImpl(~)
            % Define output port names
            names = {'Position', 'Velocity'};
        end
    end
end
