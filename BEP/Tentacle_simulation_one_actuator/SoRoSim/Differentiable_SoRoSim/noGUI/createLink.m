%This is for those who hate GUI
C3 = SorosimLink('empty'); %create an empty class element

%To create a rigid link
C3.jointtype='N';
C3.linktype='s';
C3.npie=1;

CS = 'C'; %'C' for circular, 'R' for circular, 'E' for circular
L = 2;
Rho = 1000;
r = @(X1)0.03; %function of normalzied X. X1 in [0 1]

[M,cx] = RigidBodyProperties(CS,L,Rho,r); %for circular cross section
% [M,cx] = RigidBodyProperties(CS,L,Rho,h,w); %for rectangular cross section
% [M,cx] = RigidBodyProperties(CS,L,Rho,a,b); %for elliptical cross section

C3.ld{1}=2;
C3.L= L;
C3.CS= 'C';
C3.r{1}= @(X1)0.03;
C3.h= [];
C3.w= [];
C3.a= [];
C3.b= [];
C3.cx= cx;
C3.gi = eye(4);
C3.gf = eye(4);
C3.M= M;

C3.Rho= Rho;
C3.Kj= [];
C3.Dj= [];

C3.n_l= 25;
C3.n_r= 18; %should be 5 for rectanglular cross section
C3.color= [0.9572 0.4854 0.8003];
C3.alpha= 1;
C3.CPF= false;
C3.PlotFn= @(g)CustomShapePlot(g);
C3.Lscale= 0.0947;