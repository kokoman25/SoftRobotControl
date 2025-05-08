function [pos, vel] = softRobotTestFunc(u)
% Discrete-time mass-spring-damper system
% u: input force (N)
% pos: output position (m)
% vel: output velocity (m/s)

% Load SoRoSim Model
% T1 = import('C:\Users\Robin\Documents\Robin\Delft TU\BEP\TestModelControl\matlab.mat').T1; %Load Link and Linkage 
% 
% Ts = 0.01;  % Sample time (s)
% 
% % Persistent state variables
% persistent q_prev qd_prev qqd_prev
% if isempty(q_prev)
%     q_prev = [0, 0, 0, 0]; % Initial position
%     qd_prev = [0, 0, 0, 0]; % Initial velocity
%     qqd_prev = [q_prev, qd_prev];
% end
% 
% % Compute acceleration
% actuation = @(t) deal([500 * sin(2*pi*u)], zeros(0,1), zeros(0,1));
% 
% % Update states using dynamics simulation step
% [t, qqd] = T1.dynamics(qqd_prev(end, :), actuation, dt = Ts, ...
%         t_start=0, t_end=2*Ts, Integrator='ode45');
% 
% % Average 3 of 3 time step output from simulation
% q = (qqd(1,1:4) + qqd(2,1:4) + qqd(3,1:4))/3;
% qd = (qqd(1,5:8) + qqd(1,5:8) + qqd(3,1:4))/3;
% 
% % Update persistent variables
% q_prev = q;
% qd_prev = qd;
% qqd_prev = [q, qd];
% 
% Output
in = load('matlab.mat');
Ts = 0.01;
u
pos = u;        %q;
vel = u + Ts*u; %qd;
end