%% QUADROTOR BALANCING PENDULUM MODEL PREDICTIVE CONTROL SIMULATION
%
% MATLAB simulation of the paper A Flying Inverted Pendulum by Markus Hehn 
% and Raffaello D'Andrea 

%% INIT
clc
clear

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
h = 0.02;
sysd = c2d(sysc,h);

A = sysd.A;
B = sysd.B;
C = sysd.C;

%% MODEL PREDICTIVE CONTROL

% simulation time
T = 10;

% states_trajectory: Nx16 matrix of trajectory of 16 states

%% PLOT RESULTS
% plot 2D results
% plot_2D_plots(time, states_trajectory);
% 
% % show 3D simulation
% X = states_trajectory(:,[3 9 13 11 5 15 1 7]);
% visualize_quadrotor_trajectory(states_trajectory(:,[3 9 13 11 5 15 1 7]));
