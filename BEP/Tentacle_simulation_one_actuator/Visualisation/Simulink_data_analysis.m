load('model\matlab.mat')

% Access data
simData = out.simout;

% Open Data under names
q = simData.q.Data;
qd = simData.qdot.Data;
coordinates = simData.coordinates.Data;
ts = out.tout;

q_pos_perm = permute(q, [3,1,2]);
q = reshape(permute(q, [3,1,2]), length(q), height(q));

qd_vel_perm = permute(qd, [3,1,2]);
qd = reshape(qd_vel_perm, length(qd), height(qd));

coo = zeros(3*8, length(ts));
for n = 1:width(coo)
    n
    coo(:,n) = coordinates(:,1,n);
end

% hold On
% for n = 1:T1.nsig
%     plot(coo(3*n-2,1:end-2), coo(3*n-1,1:end-2))
% end

T1.plotqt(ts, [q, qd]); 

% subplot(4,2,1)
% plot(ts, q_pos(:,1))
% title('Position')
% 
% subplot(4,2,3)
% plot(ts, q_pos(:,2))
% 
% subplot(4,2,5)
% plot(ts, q_pos(:,3))
% 
% subplot(4,2,7)
% plot(ts, q_pos(:,4))
% xlabel('Time (s)')
% 
% subplot(4,2,2)
% plot(ts, qd_vel(:,1))
% title('velocity')
% 
% subplot(4,2,4)
% plot(ts, qd_vel(:,2))
% 
% subplot(4,2,6)
% plot(ts, qd_vel(:,3))
% 
% subplot(4,2,8)
% plot(ts, qd_vel(:,4))
% xlabel('Time (s)')