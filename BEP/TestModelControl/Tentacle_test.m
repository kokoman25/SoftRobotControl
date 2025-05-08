clear;
%startup;
load('matlab.mat');
%open('DynamcicsSolution.mat')

% Definition of Time simulation
simTim = 0.1;
delt = 0.01;
tim = 0:delt:simTim;

% Initial conditions
x0 = [0, 0, 0, 0];  % initial condition position
xd0 = [0, 0, 0, 0]; % initial condition velocity
qqd = [x0, xd0];     % Initial condition

% Initialization of actuation function
actuation = @(t) deal([500 * sin(2*pi*t)], zeros(0,1), zeros(0,1));

% Simulation data storage
t_list = [0];   % time list for simulation output 
q = [x0];       % position list for simulation output
qd  = [xd0];    % Velocity list for simulation output

for i = 1:(length(tim)-2)
    [t, qqd] = T1.dynamics(qqd(end, :), actuation, dt = delt, ...
        t_start=tim(i), t_end=tim(i+2), Integrator='ode45');
    t = tim(i+1);
    t_list = [t_list; t];
    q = [q; (qqd(1,1:4) + qqd(2,1:4) + qqd(3,1:4))/3];
    qd = [qd; (qqd(1,5:8)+qqd(2,5:8)+qqd(3,5:8))/3];
end

figure;
subplot(2,2,1);
plot(t_list, q)
legend
title('Position')

subplot(2,2,2);
plot(t_list, qd)
legend
title('Velocity')
xlabel('Time (s)')

subplot(2,2,3);
plot(q(:,1), q(:,2));

subplot(2,2,4);
plot(q(:,3), q(:,4));
% subplot(2,2,3);
% T1.plotq(q(1,:))
% 
% subplot(2,2,4);
% T1.plotq(q(end,:))
