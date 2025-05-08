trajec = [];
tipTrajectory = [];

for i = 1:(length(q_pos))
    tipTrajectory = [];
    tip_co = T1.FwdKinematics(q_pos(i,:));
    for j = 1:(length(tip_co)/4)
        p = p + 1;
        tipTrajectory = [tipTrajectory; tip_co(4*j-3:4*j-1,4)];
    end
    if isempty(trajec)
        trajec = [tipTrajectory];
    else
        trajec = [trajec, tipTrajectory];
    end
    size(trajec)
end
T1.dynamicsSolver()
hold on
for n = 1:height(trajec)/3
    plot(trajec(3*n-2,:), trajec(3*n-1,:))
end
figure

% plot3(trajec(19,:), trajec(20,:), trajec(21,:))
% plot3(trajec(22,:), trajec(23,:), trajec(24,:))
grid
xlabel('X (m)')
ylabel('Y (m)')
hold off
%xlim([-0.5, 0.5])
%ylim([-0.5, 0.5])