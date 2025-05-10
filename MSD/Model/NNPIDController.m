% Run this file only after running MSD Simulink model
% This file expects data in the workspace which is outputted by the
% Simulink model.
%% Access Data
y = squeeze(out.simout);    % Data vector
t = out.tout;               % Time vector

%% Assign data to variables
u = y(1,:)';     % Control signal vector
q = y(2,:)';     % Position vector
qd = y(3,:)';    % Velocity vector
r = y(4,:)';     % Reference signal vector

%% Plot data vs time
figure;
plot(t, u, 'DisplayName','Control'); hold on
plot(t, q, 'DisplayName','Position');
plot(t, qd, 'DisplayName','Velocity');
plot(t, r, 'DisplayName','Reference');
legend

%% Prepare data for Training
input = [r, q, qd];
output = u;

%% Train Neural Network
net = feedforwardnet([10 10 10]);
net.layers{1}.transferFcn = 'logsig';
net.layers{2}.transferFcn = 'purelin';
net.layers{3}.transferFcn = 'tansig';
net.trainParam.showWindow = true;
net = train(net, input', output');