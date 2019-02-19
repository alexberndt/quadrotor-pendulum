%% QUADROTOR BALANCING PENDULUM MODEL PREDICTIVE CONTROL SIMULATION
%
% MATLAB simulation of the paper A Flying Inverted Pendulum by Markus Hehn 
% and Raffaello D'Andrea 

%% INITIALIZATION 
clc
clear

% Constants
g = 9.81;       % m/s^2
L = 0.565;      % meters (Length of pendulum to center of mass)
l = 0.17;       % meters (Quadrotor center to rotor center)
I_yy = 3.2e-3;  % kg m^2 (Quadrotor inertia around y-axis)
I_xx = I_yy;
I_zz = 5.5e-3;  % kg m^2 (Quadrotor inertia around z-axis) 

% continuous-time state space matrices
Ac = [0   1 0 0  0 0 ;
     g/L  0 0 0 -g 0 ;
     0    0 0 1  0 0 ;
     0    0 0 0  g 0 ;
     0    0 0 0  0 1 ;
     0    0 0 0  0 0];
 
Bc = [0         0;
      0         0;
      0         0;  
      0         0;  
      0         0;
      l/I_yy -l/I_yy];

Cc = eye(size(Ac));
% outputs we are interested are
% - r1:     the displacement of 
% - r2:     velocity of pendulum relative to quadrotor
% - x1:     x-direction displacement of quadrotor relative to inertial
%           coordinate frame O
% - x2:     x-direction velocity of quadrotor relative to inertial
%           coordinate frame O
% - beta:   pitch angle of quad (rotation angle around y-axis)

% continuous-time state-space model
sysc = ss(Ac,Bc,Cc,[]);
 
%% CHECK CONTROLLABILITY OF CONTINUOUS TIME SYSTEM
 
Ctrb_rank = rank(ctrb(Ac,Bc));
disp('Rank of controllability matrix');
disp(Ctrb_rank);

% not controllable
% use Hautus test 

eigvals = eig(Ac);

for idx = 1:numel(eigvals)
 lambda = eigvals(idx);
 disp('Eigenvalue: ');
 disp(lambda)
 rk = rank([(eye(size(Ac))*lambda-Ac) Bc]);
 disp('rank: ');
 disp(rk);
 disp('----------------');
end

%% DISCRETIZE SYSTEM

Ts = 0.05;  % 10 ms sampling time
sysd = c2d(sysc,Ts);

A = sysd.A;
B = sysd.B;
C = sysd.C;

%% LINEAR QUADRATIC REGULATOR

% we have three uncontrollable modes lambda = 0,0,0

Q = eye(size(A));
R = 0.1*eye(length(B(1,:)));

[K,S,e] = dlqr(A,B,Q,R,[]);

%% SIMUALTION

sysd_closedloop = ss(A-B*K,B,C,[],Ts);

% [states_trajectory, time] = step(sysd_closedloop);

T = 10;
dt = Ts;
u = [zeros(2,T/dt) [0;0]];
t = 0:dt:T;
x0 = [0.05;0.05;0;0;0;0];

y = lsim(sysd_closedloop,u,t,x0);

states_trajectory = y;
time = t;

%% PLOT SIMULATION RESULTS
% 
% Simple 2D plots to quickly see the performance characteristics of each
% decoupled controller

% Show 5 States of x-direction control

show_x_horizontal_performance_plots = true;
if show_x_horizontal_performance_plots
    figure(1);
    clf;
    subplot 511;
    stairs(time, states_trajectory(:,1), 'm-'); grid();
    ylabel('$r_1$ [m]','interpreter','latex')
    
    subplot 512;
    stairs(time, states_trajectory(:,2), 'm-');  grid();
    ylabel('$r_2$ [m/s]','interpreter','latex')
    
    subplot 513;
    stairs(time, states_trajectory(:,3), 'b-');  grid();
    ylabel('$x_1$ [m]','interpreter','latex')
    
    subplot 514;
    stairs(time, states_trajectory(:,4), 'b-');  grid();
    ylabel('$x_2$ [m/s]','interpreter','latex')
    
    subplot 515;
    stairs(time, states_trajectory(:,5), 'k-');  grid();
    ylabel('$\beta$ [rad]','interpreter','latex')
    
    xlabel('Time [s]');
end

%% PLOT SIMULATION IN 3D
%
% A full 3D visualization of the quadrotor given previous simulations

show_3D_visualization = true;
if show_3D_visualization 
    x =     states_trajectory(:,3);
    y =     states_trajectory(:,3);
    z =     2*ones(length(time),1);

    roll =  states_trajectory(:,5);
    pitch = states_trajectory(:,5);
    yaw =   zeros(length(time),1);
    
    r =     states_trajectory(:,1);
    s =     states_trajectory(:,1);

    X =     [x,y,z,roll,pitch,yaw,r,s];
    visualize_quadrotor_trajectory(X);
end















































% .