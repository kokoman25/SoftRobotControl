function [pos, vel] = massSpringDamperDynamics(u)
% Discrete-time mass-spring-damper system
% u: input force (N)
% pos: output position (m)
% vel: output velocity (m/s)

% System parameters
m = 1.0;    % Mass (kg)
c = 0.5;    % Damping coefficient (Ns/m)
k = 2.0;    % Spring constant (N/m)
Ts = 0.01;  % Sample time (s)

% Persistent state variables
persistent x_prev v_prev
if isempty(x_prev)
    x_prev = 0; % Initial position
    v_prev = 0; % Initial velocity
end

% Compute acceleration
a = (u - c * v_prev - k * x_prev) / m;

% Update states using Euler integration
v = v_prev + Ts * a;
x = x_prev + Ts * v;

% Update persistent variables
v_prev = v;
x_prev = x;

% Output
pos = x;
vel = v;
end