mdl = 'SRT2ARL';
agentBlk = [mdl, '/RL Agent'];

% Define observation and action specifications
obsInfo = rlNumericSpec([12 1]);
obsInfo.Name = 'error_ref_ypos_xpos_avel_ang';

actInfo = rlNumericSpec([2 1], 'LowerLimit', [-10; -10], 'UpperLimit', [0; 0]);
actInfo.Name = 'ac1_ac2';

% Create the environment
env = rlSimulinkEnv(mdl, agentBlk, obsInfo, actInfo);
env.ResetFcn = @(in) msdResetFcn(in);
env.UseFastRestart = 'off';

%% RL Agent Design
% Define the observation and action dimensions
numObs = obsInfo.Dimension(1);
numAct = actInfo.Dimension(1);

% Create the actor network
actorLayerSizes = [10, 10];
actorNetwork = [
    featureInputLayer(numObs, 'Normalization', 'none', 'Name', 'state')
    fullyConnectedLayer(actorLayerSizes(1), 'Name', 'fc1')
    reluLayer('Name', 'relu1')
    fullyConnectedLayer(actorLayerSizes(2), 'Name', 'fc2')
    reluLayer('Name', 'relu2')
    fullyConnectedLayer(numAct, 'Name', 'fc3')
    tanhLayer('Name', 'tanh1')
    scalingLayer('Name', 'ActorScaling', 'Scale', max(abs(actInfo.UpperLimit)))
    ];

actorOptions = rlRepresentationOptions('LearnRate',1e-4,'GradientThreshold',1, 'UseDevice', 'gpu');
actor = rlDeterministicActorRepresentation(actorNetwork, obsInfo, actInfo, ...
    'Observation', {'state'}, 'Action', {'ActorScaling'}, actorOptions);

% Create the critic network
statePath = [
    featureInputLayer(numObs, 'Normalization', 'none', 'Name', 'state')
    fullyConnectedLayer(24, 'Name', 'fc1')
    reluLayer('Name', 'relu1')];

actionPath = [
    featureInputLayer(numAct, 'Normalization', 'none', 'Name', 'action')
    fullyConnectedLayer(24, 'Name', 'fc2')
    ];

commonPath = [
    additionLayer(2, 'Name', 'add')
    reluLayer('Name', 'relu2')
    fullyConnectedLayer(1, 'Name', 'fc3')];

criticNetwork = layerGraph(statePath);
criticNetwork = addLayers(criticNetwork, actionPath);
criticNetwork = addLayers(criticNetwork, commonPath);
criticNetwork = connectLayers(criticNetwork, 'relu1', 'add/in1');
criticNetwork = connectLayers(criticNetwork, 'fc2', 'add/in2');

criticOptions = rlRepresentationOptions('LearnRate',1e-3,'GradientThreshold',1, 'UseDevice', 'gpu');
critic = rlQValueRepresentation(criticNetwork, obsInfo, actInfo, ...
    'Observation', {'state'}, 'Action', {'action'}, criticOptions);

% Create the DDPG agent
agentOptions = rlDDPGAgentOptions(...
    'SampleTime', 0.01, ...
    'TargetSmoothFactor', 1e-3, ...
    'ExperienceBufferLength', 1e6, ...
    'MiniBatchSize', 64, ...
    'DiscountFactor', 0.99);

agent = rlDDPGAgent(actor, critic, agentOptions);

%% Training
trainOpts = rlTrainingOptions(...
    'MaxEpisodes', 10, ...
    'MaxStepsPerEpisode', 200, ...
    'ScoreAveragingWindowLength', 5, ...
    'Verbose', false, ...
    'Plots', 'training-progress', ...
    'StopTrainingCriteria', 'AverageReward', ...
    'StopTrainingValue', 500);

% trainingStats = train(agent, env, trainOpts);

% in = Simulink.SimulationInput(mdl)';
% in = msdResetFcn(in);
% sim(in);
sim('SRT2ARL')