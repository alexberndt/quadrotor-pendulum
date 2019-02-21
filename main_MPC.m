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
h = 0.2;
sysd = c2d(sysc,h);

A = sysd.A;
B = sysd.B;
C = sysd.C;

%% MODEL PREDICTIVE CONTROL

% simulation time
t_final = 8;
x0 = [0.02 0 0.1 0 0 0  0.05 0 0.4 0 0 0  0.2 0  0.3 0];

% reference sequence
r = [ 0*ones(1,((t_final/h)+1));
      0*ones(1,((t_final/h)+1));
      0.5*ones(1,((t_final/h)+1));
      1.0*ones(1,((t_final/h)+1))];
  
Q_lqr = eye(size(A));
R_lqr = 0.1*eye(length(B(1,:)));

[K,~,~] = dlqr(A,B,Q_lqr,R_lqr,[]); 

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


T = t_final/h;
x = zeros(length(A(:,1)),T);
u = zeros(length(B(1,:)),T);
y = zeros(length(C(:,1)),T);
t = zeros(1,T);

x0 = x0';
x(:,1) = x0';

% Define MPC Control Problem

% MPC cost function
%          N-1
% V(u_N) = Sum 1/2[ x(k)'Qx(k) + u(k)'Ru(k) ] + x(N)'Sx(N) 
%          k = 0

% tuning weights
Q = 10*eye(size(A));
R = 0.1*eye(length(B(1,:)));
S = 10*eye(size(A));

% prediction horizon
N = 8; 

Qbar = kron(Q,eye(N));
Rbar = kron(R,eye(N));
Sbar = S;

P = [eye(size(A,1)); A; A^2; A^3; A^4; A^5; A^6; A^7];
z = zeros(16,4);
Z = [  z     z     z     z     z     z    z   z ;
       B     z     z     z     z     z    z   z ;
      A*B    B     z     z     z     z    z   z ;
     A^2*B  A*B    B     z     z     z    z   z ;
     A^3*B A^2*B  A*B    B     z     z    z   z ;
     A^4*B A^3*B A^2*B  A*B    B     z    z   z ;
     A^5*B A^4*B A^3*B A^2*B  A*B    B    z   z ; 
     A^6*B A^5*B A^4*B A^3*B A^2*B  A*B   B   z ];
W = [A^7*B A^6*B A^5*B A^4*B A^3*B A^2*B  A*B   B];
              
H = (Z'*Qbar*Z + Rbar + 2*W'*Sbar*W);
d = (x0'*P'*Qbar*Z + 2*x0'*(A^N)'*Sbar*W)';
 
%%

for k = 1:1:T
    t(k) = (k-1)*h;
    
    % compute control action
    x0 = x(:,k);
    d = (x0'*P'*Qbar*Z + 2*x0'*(A^N)'*Sbar*W)';
    
    cvx_begin quiet
        variable u_N(4*N)
        minimize ( (1/2)*quad_form(u_N,H) + d'*u_N )
        u_N >= -1000;
        u_N <=  1000;
    cvx_end
    
    u(:,k) = u_N(1:4); % MPC control action
    
    % u(:,k) = -K*x(:,k); % LQR control action
    
    % apply control action
    x(:,k+1) = A*x(:,k) + B*u(:,k); % + B_ref*r(:,k);
    y(:,k) = C*x(:,k);
end

% states_trajectory: Nx16 matrix of trajectory of 16 states
states_trajectory = y';

%% PLOT RESULTS
% plot 2D results
% plot_2D_plots(t, states_trajectory);

% show 3D simulation
X = states_trajectory(:,[3 9 13 11 5 15 1 7]);
visualize_quadrotor_trajectory(states_trajectory(:,[3 9 13 11 5 15 1 7]),0.1);
