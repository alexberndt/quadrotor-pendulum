%% QUADROTOR BALANCING PENDULUM MODEL PREDICTIVE CONTROL SIMULATION
%
% MATLAB simulation of the paper A Flying Inverted Pendulum by Markus Hehn 
% and Raffaello D'Andrea using a Model Predictive Controller

%% INIT
clc
clear
addpath('functions/');

%% DEFINE CONSTANTS
g = 9.81;       % m/s^2
m = 0.5;        % kg
L = 0.565;      % meters (Length of pendulum to center of mass)
l = 0.17;       % meters (Quadrotor center to rotor center)
I_yy = 3.2e-3;  % kg m^2 (Quadrotor inertia around y-axis)
I_xx = I_yy;    
I_zz = 5.5e-3;  % kg m^2 (Quadrotor inertia around z-axis)

%% DEFINE STATE SPACE SYSTEM
sysc = init_system_dynamics(g,m,L,l,I_xx,I_yy,I_zz);
check_controllability(sysc);

%% DISCRETIZE SYSTEM
h = 0.1;
sysd = c2d(sysc,h);

A = sysd.A;
B = sysd.B;
C = sysd.C;

%% MODEL PREDICTIVE CONTROL

% simulation time
T = 3.5;
x0 = [0.02 0 0.1 0 0 0  0.05 0 0.4 0 0 0  0.2 0  0.3 0];

% reference sequence
r = [ 0*ones(1,((T/h)+1));
      0*ones(1,((T/h)+1));
      0.5*ones(1,((T/h)+1));
      1.0*ones(1,((T/h)+1))];
  
Q = eye(size(A));
R = 0.1*eye(length(B(1,:)));

[K,S,e] = dlqr(A,B,Q,R,[]); 

B_ref = zeros(16,4);
B_ref(3,1) = 1;
B_ref(9,2) = 1;
B_ref(13,3) = 1;
B_ref(15,4) = 1;

% define closed loop system with LQR control law
sysd_cl_unnormalized = ss(A-B*K,B_ref,C,[],h);

% normalize closed-loop reference tracking gains 
dcgain_cl = dcgain(sysd_cl_unnormalized);
B_ref(3,1) = 1/dcgain_cl(3,1);
B_ref(9,2) = 1/dcgain_cl(9,2);
B_ref(13,3) = 1/dcgain_cl(13,3);
B_ref(15,4) = 1/dcgain_cl(15,4);


N = T/h;
x = zeros(length(A(:,1)),N);
u = zeros(length(B(1,:)),N);
y = zeros(length(C(:,1)),N);
t = zeros(1,N);

x(:,1) = x0';

for k = 1:1:N
    t(k) = (k-1)*h;
    
    % compute control action
    u(:,k) = -K*x(:,k);   
    
    % apply control action
    x(:,k+1) = A*x(:,k) + B*u(:,k) + B_ref*r(:,k);
    y(:,k) = C*x(:,k);
end

% states_trajectory: Nx16 matrix of trajectory of 16 states
states_trajectory = y';

%% PLOT RESULTS
% plot 2D results
plot_2D_plots(t, states_trajectory);

% show 3D simulation
X = states_trajectory(:,[3 9 13 11 5 15 1 7]);
visualize_quadrotor_trajectory(states_trajectory(:,[3 9 13 11 5 15 1 7]));
