% Setup
n_steps = 200; 
Ts = 0.01;     

% Initial states
position = zeros(1, n_steps);
velocity = zeros(1, n_steps);
control  = zeros(1, n_steps);

position(1) = 0; % initial
velocity(1) = 0; % initial
control(1)  = 0; % initial

% Target
target_position = 1.0;

% Controller gains (simple proportional-derivative controller)
Kp = 10;   % Proportional gain
Kd = 2;    % Derivative gain

for k = 1:n_steps-1
    % Real feedback: current position and velocity
    current_position = position(k);
    current_velocity = velocity(k);
    
    % Error to target
    position_error = target_position - current_position;
    velocity_error = -current_velocity; % we want velocity -> 0
    
    % Simple PD control
    control(k) = Kp * position_error + Kd * velocity_error;
    
    % Build input for neural network
    nn_input = [control(k); current_velocity; current_position];
    
    % Predict next step
    nn_output = net(nn_input);
    
    % Update state
    control(k+1)  = nn_output(1); % predicted control for consistency (not used here)
    velocity(k+1) = nn_output(2);
    position(k+1) = nn_output(3);
end

% Time vector
time = (0:n_steps-1)*Ts;

% Plot results
figure;
subplot(3,1,1);
plot(time, position, 'b', 'LineWidth', 1.5);
yline(target_position, 'r--', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Position (m)');
title('Closed-Loop Position Control');
%legend('Position', 'Target');

subplot(3,1,2);
plot(time, velocity, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Velocity (m/s)');
title('Velocity');

subplot(3,1,3);
plot(time, control(1:end), 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Control Force (N)');
title('Control Input');
