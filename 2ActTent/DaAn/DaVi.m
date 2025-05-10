signals = squeeze(out.simout);
time = out.tout;

q = signals(1:4,:)';
qd = signals(5:8,:)';

Tentacle.plotqt(time, [q, qd])