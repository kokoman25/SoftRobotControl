clear;                              % Clear Workspace to make sure that only data from this script is used 
%startup;                           % Initialize SoRoSim Toolbox
load('model\TentacleFiles.mat');    % Load model paramters into the workspace
%open('DynamcicsSolution.mat')

% Definition of Time simulation
T = 10;              % Total duration of simulation (s)
dt = 0.01;          % Time step (s)
Ts = 0:dt:T;        % List of time steps until end of simulation

% Initial conditions
q0 = [0, 0, 0, 0];      % initial condition position
qd0 = [0, 0, 0, 0];     % initial condition velocity
qqd = [q0, qd0];        % Initial condition

% Initialization of actuation function
actuation = @(t) deal([5 * sin(2*pi*t); 5*cos(4*pi*t)], zeros(0,1), zeros(0,1)); 

% Simulation data storage
t_list = zeros(1,length(Ts));                           % Time list for simulation output 
q = zeros(length(Ts),length(q0));      q(1,:) = q0;     % Angle list for simulation output
qd = zeros(length(Ts), length(qd0));   qd(1,:) = qd0;   % Angluar velocity list for simulation output
Bco = zeros(3*Tentacle.nsig, length(Ts));               % Body coordinate list

% Simulate Dynamics of Soft Body
for i = 1:(length(Ts)-2)

    % Calculate Dynamic state vector
    [t, qqd] = Tentacle.dynamics(qqd(end, :), actuation, dt = dt, ...
        t_start=Ts(i), t_end=Ts(i+2), Integrator='ode45');      % Dynamic simulation of 1 timestep, output is 3 dynamic solutions for 3 timesteps
    t_list(i) = t(2);                                           % add current time to time list

    % Get angle and angluar velocity and add to state space vector
    q(i,:) = (qqd(1,1:4)+qqd(2,1:4)+qqd(3,1:4))/3;  % Take average of angle simulation output and add to list
    qd(i,:) = (qqd(1,5:8)+qqd(2,5:8)+qqd(3,5:8))/3; % Take average of angluar velocity 

    % Get coordinates of significant points along soft body
    point_co = Tentacle.FwdKinematics(q(i,:));
    for j = 1:(length(point_co)/4)
        point_co(4*j-3:4*j-1,4);
        Bco(3*j-2:3*j,i) = point_co(4*j-3:4*j-1,4);
    end

    % Compute the Jacobian matrix at all significant points
    J = Tentacle.Jacobian(q(i,:));
    
    % Determine the number of significant points
    n_sig = Tentacle.nsig;
    
    % Determine the number of degrees of freedom
    n_dof = Tentacle.ndof;
    
    % Extract the Jacobian for the last significant point
    % Each point's Jacobian occupies 6 rows in the matrix
    J_tip = J((n_sig-1)*6 + 1 : n_sig*6, :);

    t = Ts(i+1); % progress timestep
end


% hold On
% for n = 1:T1.nsig
%     plot(Trajec(3*n-2,1:end-2), Trajec(3*n-1,1:end-2))
% end
% 
% size(Ts(1:end))
% size([q, qd])

Tentacle.plotqt(Ts(1:end), [q, qd])
 
% % figure;
% % subplot(2,2,1);
% % plot(t_list, q)
% % legend
% % title('Position')
% % 
% % subplot(2,2,2);
% % plot(t_list, qd)
% % legend
% % title('Velocity')
% % xlabel('Time (s)')
% % 
% % subplot(2,2,3);
% % plot(q(:,1), q(:,2));
% % 
% % subplot(2,2,4);
% % plot(q(:,3), q(:,4));
% % subplot(2,2,3);
% % T1.plotq(q(1,:))
% % 
% % subplot(2,2,4);
% % T1.plotq(q(end,:))
% %
% % %lists for Dynamic solver
% % y_list = [];
% % C_list = [];
% % B_action_list = [];
% % action_list = [];
% %
% % See what Dynamic solver does
% % [y,C,B_action,action] = T1.dynamicsSolver(tim(1), qqd(1,:)', 10);
% % y(1:T1.ndof)
% % y_list = [y_list, y];
% % C_list = [C_list, C];
% % B_action_list = [B_action_list, B_action];
% % action_list = [action_list, action];