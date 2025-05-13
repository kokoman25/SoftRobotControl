% Setup
n_steps = 200;  % number of time setps
Ts = 0.01;      % sample time

%initial states 
position = zeros(1, n_steps);
velocity = zeros(1, n_steps);
control = zeros(1, n_steps);

position(1) = 0;
velocity(1) = 0;
control(1) = 0;

% Target
target_position = 1.0;

% Controller gains (simple proportional controller)
Kp = 5; % Proportional gain

% Control loop
for k = 1:n_steps-1
    % Error to target
    error = target_position - position(k);

    % Simple proportional control
    control(k) = Kp * error;

    % Build input for neural network
    nn_input= [control(k); velocity(k); position(k)];

    % Predict next state
    nn_output = net(nn_input);

    % Update state
    control(k+1) = nn_output(1);
    velocity(k+1) = nn_output(2);
    position(k+1) = nn_output(3);
end

% Time vector
time = (0:n_steps-1)*Ts;

% Plot results
figure;
subplot(3,1,1);
plot(time, position, 'b', 'LineWidth', 1.5);
yline(target_position, 'r--');
xlabel('Time (s)');
ylabel('Position (m)');
title('Position Tracking');
legend('Predicted Position', 'Target');

subplot(3,1,2);
plot(time, velocity, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Velocity (m/s)');
title('Velocity');

subplot(3,1,3);
plot(time, control, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Control Force (N)');
title('Control Input');