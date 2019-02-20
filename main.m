%% QUADROTOR BALANCING PENDULUM MODEL PREDICTIVE CONTROL SIMULATION
%
% MATLAB simulation of the paper A Flying Inverted Pendulum by Markus Hehn 
% and Raffaello D'Andrea 

%% INITIALIZATION 
clc
clear

% Constants
g = 9.81;       % m/s^2
m = 0.5;        % kg
L = 0.565;      % meters (Length of pendulum to center of mass)
l = 0.17;       % meters (Quadrotor center to rotor center)
I_yy = 3.2e-3;  % kg m^2 (Quadrotor inertia around y-axis)
I_xx = I_yy;
I_zz = 5.5e-3;  % kg m^2 (Quadrotor inertia around z-axis) 

% SUBSYSTEM 1 - pitch angle dynamics (around y-axis)
Ac1 = [  0    1 0 0  0 0 ;
         g/L  0 0 0 -g 0 ;
         0    0 0 1  0 0 ;
         0    0 0 0  g 0 ;
         0    0 0 0  0 1 ;
         0    0 0 0  0 0 ];
 
Bc1 = [0         0;
       0         0;
       0         0;  
       0         0;  
       0         0;
      -l/I_yy  l/I_yy];

Cc1 = eye(size(Ac1));
% outputs we are interested are
% - r1:     the displacement of 
% - r2:     velocity of pendulum relative to quadrotor
% - x1:     x-direction displacement of quadrotor relative to inertial
%           coordinate frame O
% - x2:     x-direction velocity of quadrotor relative to inertial
%           coordinate frame O
% - beta1:  pitch angle of quad (rotation angle around y-axis)
% - beta2:  pitch angle rate of quad (rotation angle around y-axis)

% SUBSYSTEM 2 - pitch angle dynamics (around x-axis)
Ac2 = [   0   1 0 0  0  0 ;
         g/L  0 0 0 -g  0 ;
         0    0 0 1  0  0 ;
         0    0 0 0 -g  0 ;
         0    0 0 0  0  1 ;
         0    0 0 0  0  0 ];
 
Bc2 = [0         0;
       0         0;
       0         0;  
       0         0;  
       0         0;
      -l/I_xx  l/I_xx];

Cc2 = -eye(size(Ac2)); % negative to have positive y-direction
% outputs we are interested are
% - s1:     the displacement of 
% - s2:     velocity of pendulum relative to quadrotor
% - y1:     x-direction displacement of quadrotor relative to inertial
%           coordinate frame O
% - y2:     x-direction velocity of quadrotor relative to inertial
%           coordinate frame O
% - gamma1: pitch angle of quad (rotation angle around x-axis)
% - gamma2: pitch angle rate of quad (rotation angle around x-axis)

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

% % IF not controllable
% % use Hautus test 
% 
% eigvals = eig(Ac);
% 
% for idx = 1:numel(eigvals)
%  lambda = eigvals(idx);
%  disp('Eigenvalue: ');
%  disp(lambda)
%  rk = rank([(eye(size(Ac))*lambda-Ac) Bc]);
%  disp('rank: ');
%  disp(rk);
%  disp('----------------');
% end

%% DISCRETIZE SYSTEM

Ts = 0.02;  % 10 ms sampling time
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

B_ref = [0 0 0 0;
         0 0 0 0;
         1 0 0 0;
         0 0 0 0;
         0 0 0 0;
         0 0 0 0;
         0 0 0 0;
         0 0 0 0;
         0 1 0 0;
         0 0 0 0;
         0 0 0 0;
         0 0 0 0;
         0 0 1 0;
         0 0 0 0;
         0 0 0 1;
         0 0 0 0];

sysd_closedloop = ss(A-B*K,B_ref,C,[],Ts);

dcgain_cl = dcgain(sysd_closedloop);

B_ref(3,1) = 1/dcgain_cl(3,1);
B_ref(9,2) = 1/dcgain_cl(9,2);
B_ref(13,3) = 1/dcgain_cl(13,3);
B_ref(15,4) = 1/dcgain_cl(15,4);

sysd_closedloop_adjusted = ss(A-B*K,B_ref,C,[],Ts);

dcgain(sysd_closedloop_adjusted)

% [states_trajectory, time] = step(sysd_closedloop);

%%

T = 5;
dt = Ts;
u = 1*ones(3,((T/dt)+1));

u = [0*ones(1,((T/dt)+1));
     0*ones(1,((T/dt)+1));
     0*ones(1,((T/dt)+1));
     0*ones(1,((T/dt)+1))];

t = 0:dt:T;
% initial pendulum offset, quadrotor offset and quadrotor yaw angle
x0 = [0.02 0 0.1 0 0 0  0.05 0 0.4 0 0 0  0.2 0  0.3 0];

y = lsim(sysd_closedloop_adjusted,u,t,x0);

states_trajectory = y;
time = t;

%% PLOT SIMULATION RESULTS
% 
% Simple 2D plots to quickly see the performance characteristics of each
% decoupled controller

% Show 6 States of x-direction control
show_x_horizontal_performance_plots = true;
if show_x_horizontal_performance_plots
    figure(1);
    clf;
    sgtitle('x-direction horizontal motion and pitch angles');
    subplot 511;
    stairs(time, states_trajectory(:,1), 'm-'); grid();
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

% Show 6 States of x-direction control
show_y_horizontal_performance_plots = true;
if show_y_horizontal_performance_plots
    figure(2);
    clf;
    sgtitle('y-direction horizontal motion and roll angles');
    subplot 511;
    stairs(time, states_trajectory(:,7), 'm-'); grid();
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