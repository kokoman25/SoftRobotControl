% Access data
simData = out.simData;
t = simData.time;
y = simData.Data;

% Position and Velocity
control = y(:,1);
velocity = y(:,2);
position = y(:,3);
reference = y(:,4);
% 
% input = [reference, position, velocity];
% output = control;
% % 
% net = feedforwardnet([10 10 10]);
% net.layers{1}.transferFcn = 'logsig';
% net.layers{2}.transferFcn = 'purelin';
% net.layers{3}.transferFcn = 'tansig';
% net.trainParam.showWindow = true;
% net = train(net, input', output');