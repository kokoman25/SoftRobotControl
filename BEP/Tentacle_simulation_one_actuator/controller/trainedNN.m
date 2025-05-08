new_input = input';

current_state = new_input(:,1);

num_steps = 2000;
dt = 0.01;
time = 1:dt:((num_steps*dt)-dt+1);
predicted_rollout = zeros(size(current_state, 1), num_steps);

for k = 1:num_steps
    predicted_next = net(current_state);
    predicted_rollout(:,k) = predicted_next;

    %update current state
    current_state = predicted_next;
end

predicted_control = predicted_rollout(1,:);
predicted_velocity = predicted_rollout(2,:);
predicted_position = predicted_rollout(3,:);

figure;
subplot(3,1,1);
plot(time, predicted_control(1,:))
title('control prediction')

subplot(3,1,2);
plot(time, predicted_velocity(1,:))
title('velocity prediction')

subplot(3,1,3)
plot(time, predicted_position(1,:))
title('position prediction')
xlabel('Time (s)')


% true_control = output(:, 1)';
% true_velocity = output(:, 2)';
% true_position = output(:, 3)';
% 
% time = t(2:end);
% 
% % Plot results
% figure;
% subplot(3,1,1);
% plot(time, true_control, 'b', time, predicted_control, 'r--', 'LineWidth', 1.5);
% legend('True Control', 'NN Predicted Control');
% title('Control Prediction');
% 
% subplot(3,1,2);
% plot(time, true_velocity, 'b', time, predicted_velocity, 'r--', 'LineWidth', 1.5);
% legend('True Velocity', 'NN Predicted Velocity');
% title('Velocity Prediction');
% 
% subplot(3,1,3);
% plot(time, true_position, 'b', time, predicted_position, 'r--', 'LineWidth', 1.5);
% legend('True Position', 'NN Predicted Position');
% title('Position Prediction');
% xlabel('Time (s)');

