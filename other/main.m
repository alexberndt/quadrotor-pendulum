%% QUADROTOR BALANCING PENDULUM MODEL PREDICTIVE CONTROL SIMULATION
%
% MATLAB simulation of the paper A Flying Inverted Pendulum by Markus Hehn 
% and Raffaello D'Andrea 

%% INITIALIZATION AND DEFINING SYSTEM DYNAMICS
clc
clear

% DEFINE CONSTANTS
g = 9.81;       % m/s^2
m = 0.5;        % kg
L = 0.565;      % meters (Length of pendulum to center of mass)
l = 0.17;       % meters (Quadrotor center to rotor center)
I_yy = 3.2e-3;  % kg m^2 (Quadrotor inertia around y-axis)
I_xx = I_yy;    
I_zz = 5.5e-3;  % kg m^2 (Quadrotor inertia around z-axis) 

% SUBSYSTEM 1 - pitch angle dynamics (around y-axis)
Ac1 = [  0    1 0 0  0 0 ;      % - r1:     the displacement of
         g/L  0 0 0 -g 0 ;      % - r2:     velocity of pendulum relative to quadrotor
         0    0 0 1  0 0 ;      % - x1:     x-direction displacement of quadrotor relative to inertial coordinate frame O
         0    0 0 0  g 0 ;      % - x2:     x-direction velocity of quadrotor relative to inertial coordinate frame O
         0    0 0 0  0 1 ;      % - beta1:  pitch angle of quad (rotation angle around y-axis)
         0    0 0 0  0 0 ];     % - beta2:  pitch angle rate of quad (rotation angle around y-axis)
 
Bc1 = [0         0;
       0         0;
       0         0;  
       0         0;  
       0         0;
      -l/I_yy  l/I_yy];

Cc1 = eye(size(Ac1));

% SUBSYSTEM 2 - pitch angle dynamics (around x-axis)
Ac2 = [   0   1 0 0  0  0 ;     % - s1:     the displacement of
         g/L  0 0 0 -g  0 ;     % - s2:     velocity of pendulum relative to quadrotor
         0    0 0 1  0  0 ;     % - y1:     x-direction displacement of quadrotor relative to inertial coordinate frame O
         0    0 0 0 -g  0 ;     % - y2:     x-direction velocity of quadrotor relative to inertial coordinate frame O
         0    0 0 0  0  1 ;     % - gamma1: pitch angle of quad (rotation angle around x-axis)
         0    0 0 0  0  0 ];    % - gamma2: pitch angle rate of quad (rotation angle around x-axis)
     
Bc2 = [0         0;
       0         0;
       0         0;  
       0         0;  
       0         0;
      -l/I_xx  l/I_xx];
  
Cc2 = -eye(size(Ac2)); % negative to have positive y-direction

% SUBSYSTEM 3 - verticle translational dynamics
Ac3 = [0 1;
       0 0];

Bc3 = [0   0   0   0   ;
       1/m 1/m 1/m 1/m];
   
Cc3 = eye(size(Ac3));

% SUBSYSTEM 4 - yaw angle of quadrotor
Ac4 = [0 1;
       0 0];
Bc4 = [0        0         0       0   ;
       -1/I_zz  -1/I_zz  1/I_zz 1/I_zz];
Cc4 = eye(size(Ac4));

% FULL SYSTEM CONCATENATION
Ac = blkdiag(Ac1,Ac2,Ac3,Ac4);

Bc = [Bc1        zeros(6,2);
      zeros(6,2) Bc2;
      Bc3 ;
      Bc4 ];
  
Cc = blkdiag(Cc1,Cc2,Cc3,Cc4);

% continuous-time state-space model
sysc = ss(Ac,Bc,Cc,[]);
 
%% CHECK CONTROLLABILITY OF CONTINUOUS TIME SYSTEM
 
Ctrb_rank = rank(ctrb(Ac,Bc));
disp('Number of states');
disp(size(Ac));
disp('Rank of controllability matrix');
disp(Ctrb_rank);

% Hautus test
% hautus_test(Ac,Bc);

%% DISCRETIZE SYSTEM

Ts = 0.02;  
sysd = c2d(sysc,Ts);

A = sysd.A;
B = sysd.B;
C = sysd.C;

%% LINEAR QUADRATIC REGULATOR

% we have three uncontrollable modes lambda = 0,0,0

Q = eye(size(A));
R = 0.1*eye(length(B(1,:)));

[K,S,e] = dlqr(A,B,Q,R,[]);

%% DEFINE THE CLOSED LOOP SYSTEM WITH REFERENCE TRACKING INPUT

% reference tracking input control
B_ref = [0 0 0 0;  % r1
         0 0 0 0;  % r2
         1 0 0 0;  % x1 
         0 0 0 0;  % x2
         0 0 0 0;  % beta1
         0 0 0 0;  % beta2
         0 0 0 0;  % s1 
         0 0 0 0;  % s2
         0 1 0 0;  % y1
         0 0 0 0;  % y2
         0 0 0 0;  % gamma1
         0 0 0 0;  % gamma2
         0 0 1 0;  % z1
         0 0 0 0;  % z2
         0 0 0 1;  % 
         0 0 0 0]; % 

% define closed loop system with LQR control law
sysd_cl_unnormalized = ss(A-B*K,B_ref,C,[],Ts);

% normalize closed-loop reference tracking gains 
dcgain_cl = dcgain(sysd_cl_unnormalized);
B_ref(3,1) = 1/dcgain_cl(3,1);
B_ref(9,2) = 1/dcgain_cl(9,2);
B_ref(13,3) = 1/dcgain_cl(13,3);
B_ref(15,4) = 1/dcgain_cl(15,4);

% define closed-loop system with normalized input gains
sysd_cl = ss(A-B*K,B_ref,C,[],Ts);
dcgain(sysd_cl);

%% RUN SIMULATION WITH CLOSED LOOP SYSTEM

T = 5;
dt = Ts;

input_sequence = [0*ones(1,((T/dt)+1));
                  0*ones(1,((T/dt)+1));
                  0*ones(1,((T/dt)+1));
                  0*ones(1,((T/dt)+1))];
time = 0:dt:T;

% initial conditions
x0 = [0.02 0 0.1 0 0 0  0.05 0 0.4 0 0 0  0.2 0  0.3 0];

states_trajectory = lsim(sysd_cl,input_sequence,time,x0);

%% PLOT SIMULATION RESULTS
% 
% 2D plots to quickly see the performance characteristics of each
% decoupled controller

% Show 6 States of x-direction control
show_x_horizontal_performance_plots = true;
if show_x_horizontal_performance_plots
    figure(1);
    clf;
    sgtitle('x-direction horizontal motion and pitch angles');
    subplot 511;
    stairs(time, states_trajectory(:,1), 'm-');  grid();
    ylabel('$r$ [m]','interpreter','latex');
    
    subplot 512;
    stairs(time, states_trajectory(:,2), 'm-');  grid();
    ylabel('$\dot{r}$ [m/s]','interpreter','latex');
    
    subplot 513;
    stairs(time, states_trajectory(:,3), 'b-');  grid();
    ylabel('$x$ [m]','interpreter','latex');
    
    subplot 514;
    stairs(time, states_trajectory(:,4), 'b-');  grid();
    ylabel('$\dot{x}$ [m/s]','interpreter','latex');
    
    subplot 515;
    stairs(time, states_trajectory(:,5), 'k-');  grid();
    ylabel('$\beta$ [rad]','interpreter','latex');
    
    xlabel('Time [s]');
end

% Show 6 States of y-direction control
show_y_horizontal_performance_plots = true;
if show_y_horizontal_performance_plots
    figure(2);
    clf;
    sgtitle('y-direction horizontal motion and roll angles');
    subplot 511;
    stairs(time, states_trajectory(:,7), 'm-');  grid();
    ylabel('$s$ [m]','interpreter','latex');
    
    subplot 512;
    stairs(time, states_trajectory(:,8), 'm-');  grid();
    ylabel('$\dot{s}$ [m/s]','interpreter','latex');
    
    subplot 513;
    stairs(time, states_trajectory(:,9), 'b-');  grid();
    ylabel('$y$ [m]','interpreter','latex');
    
    subplot 514;
    stairs(time, states_trajectory(:,10), 'b-');  grid();
    ylabel('$\dot{y}$ [m/s]','interpreter','latex');
    
    subplot 515;
    stairs(time, states_trajectory(:,11), 'k-');  grid();
    ylabel('$\gamma$ [rad]','interpreter','latex');
    
    xlabel('Time [s]');
end

%% PLOT SIMULATION IN 3D
%
% A full 3D visualization of the quadrotor given previous simulations

show_3D_visualization = true;
if show_3D_visualization 
    x =     states_trajectory(:,3);
    y =     states_trajectory(:,9);
    z =     states_trajectory(:,13);

    pitch = states_trajectory(:,5);
    roll =  states_trajectory(:,11);
    yaw =   states_trajectory(:,15);
    
    r =     states_trajectory(:,1);
    s =     states_trajectory(:,7);

    X =     [x,y,z,roll,pitch,yaw,r,s];
    visualize_quadrotor_trajectory(X);
end















































% .