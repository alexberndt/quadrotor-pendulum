%% Disturbance Estimation

clc
clear
addpath('../functions/');

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

%% RUN SIMULATION

simTime = 10;
h = 0.02;

T = simTime/h;

t = zeros(1,T);
x = zeros(16,T);
u = zeros(4,T);
y = zeros(8,T);

sysd = c2d(sysc,h);

Ad = sysd.A;
Bd = sysd.B;
Cd = zeros(8,16);
Cd(1,1) = 1; Cd(2,3) = 1;                 % pendulum position
Cd(3,5) = 1; Cd(4,7) = 1; Cd(5,9) = 1;     % quad position
Cd(6,11) = 1; Cd(7,13) = 1; Cd(8,15) = 1;  % quad rotation angle

x0 = [0.01 0 -0.02 0 0 0 0 0 0 0 0 0 0 0 0 0];
x(:,1) = x0;

for k = 1:T
    % TIME 
    t(k) = (k-1)*h;    
    
    % LQR CONTROL ACTION
    Q_lqr = eye(16);
    R_lqr = eye(4);
    [~,L,K_lqr] = dare(Ad,Bd,Q_lqr,R_lqr);
    
    u(:,k) = -K_lqr*x(:,k);
    
    x(:,k+1) = Ad*x(:,k) + Bd*u(:,k);
    y(:,k) = Cd*x(:,k);
    
end

x = x(:,1:end-1);

%% PLOT RESULTS

figure(1);
clf;
stairs(t,x(1,:));
hold on
stairs(t,x(7,:));
legend('r','s');
grid on;

figure(2);
clf;
stairs(t,x(,:));
hold on
stairs(t,x(7,:));

stairs(t,x(9,:));

legend('x','y','z');
grid on;


































%.