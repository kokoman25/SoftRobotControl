% System parameters
m = 1.0;   % Mass (kg)
c = 0.5;   % Damping coefficient (N·s/m)
k = 2.0;   % Spring constant (N/m)

% % Make state vector list
% x = zeros(2, 1); % State vector [velocity; position]
% x0 = [0; 0]; % Initial state [velocity; position]
% x(:,1) = x0; % Initial state is first position of the state list
% 
% % Simulation parameters
% tspan = [0 10]; % Simulation time (s)
% dt = 0.01;      % Time step (s)
% t = tspan(1):dt:tspan(end); % Time vector
% N = length(t);

% % Access data
% simData = out.simData;
% t = simData.time;
% y = simData.Data;
% 
% % Position and Velocity
% position = y(:,1);
% velocity = y(:,2);
% input = y(:,3);
% 
% % Plot results
% figure;
% subplot(2,1,1);
% plot(t, position, 'LineWidth', 1.5);
% title('Mass-Spring-Damper Position');
% xlabel('Time (s)'); ylabel('Position (m)');
% 
% subplot(2,1,2);
% plot(t, velocity, 'LineWidth', 1.5);
% title('Mass-Spring-Damper Velocity');
% xlabel('Time (s)'); ylabel('Velocity (m/s)');


% % PID gains
% Kp = 10;    % Proportional gain
% Ki = 2;     % Integral gain
% Kd = 0.5;   % Derivative gain
% 
% % Initialize PID variables
% integral_error = 0;
% prev_error = 0;
% F_control = zeros(1, N); % Control force array
% ref = 1.0 * ones(size(t)); % Reference signal (step input)
% 
% sensor_noise_std = 0.01; % Sensor noise (standard deviation)
% F_max = 10;             % Actuator saturation limit (N)
% 
%  Preallocate state and measurement arrays
% x = zeros(2, N); % True states [position; velocity]
% x_measured = zeros(2, N); % Noisy measurements
% x(:,1) = x0;
% 
% % Simulation loop
% for i = 1:N-1
%     %---------------------------------------
%     % 1. Sensor measurements (with noise)
%     %---------------------------------------
%     x_measured(:,i) = x(:,i) + sensor_noise_std * randn(2,1);
% 
%     %---------------------------------------
%     % 2. Compute control force (PID)
%     %---------------------------------------
%     error = ref(i) - x_measured(1,i); % Position error
%     integral_error = integral_error + error * dt;
%     derivative_error = (error - prev_error) / dt;
% 
%     F = Kp * error + Ki * integral_error + Kd * derivative_error;
% 
%     % Actuator saturation
%     F = max(min(F, F_max), -F_max);
%     F_control(i) = F;
%     prev_error = error;
% 
%     %---------------------------------------
%     % 3. Update system dynamics (Euler integration)
%     %---------------------------------------
%     x(:,i+1) = x(:,i) + dt * [x(2,i); (F - c*x(2,i) - k*x(1,i))/m];
% end
% 
% figure;
% subplot(3,1,1);
% plot(t, x(1,:), 'b', t, ref, 'r--', 'LineWidth', 1.5);
% title('Position vs. Reference');
% xlabel('Time (s)'); ylabel('Position (m)');
% legend('Actual', 'Desired');
% 
% subplot(3,1,2);
% plot(t, F_control, 'LineWidth', 1.5);
% title('Control Force');
% xlabel('Time (s)'); ylabel('Force (N)');
% 
% subplot(3,1,3);
% plot(t, x(2,:), 'LineWidth', 1.5);
% title('Velocity');
% xlabel('Time (s)'); ylabel('Velocity (m/s)');