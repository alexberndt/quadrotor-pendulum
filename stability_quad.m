% Simple quad for MPC stability

clc
clear

%% 

I = 3e-3; % inertia
h = 0.1; % 0.01 s sampling time
M = 0.3; % kg

% define linearized system

A = [1 0 0 0 0 0;
     h 1 0 0 0 0;
     0 0 1 0 0 0;
     0 0 h 1 0 0;
     0 0 0 0 1 0;
     0 0 0 0 h 1];
 
B = [h/I -h/I;
     0     0;
     h/M  h/M;
     0     0;
     0     0;
     0     0];
 
C = ctrb(A,B)
rank(C)     
     