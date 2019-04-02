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

% simulation time in seconds
simTime = 8;
h = 0.04;

sysd = c2d(sysc,h);
T = simTime/h;

A = sysd.A;
B = sysd.B;
Cplot = sysd.C;

C = zeros(8,16);
C(1,1) = 1; C(2,3) = 1;                 % pendulum position
C(3,5) = 1; C(4,7) = 1; C(5,9) = 1;     % quad position
C(6,11) = 1; C(7,13) = 1; C(8,15) = 1;  % quad rotation
% C(9,12) = 1;
% C(10,14) = 1;
% C(11,15) = 1;

%% MODEL PREDICTIVE CONTROL

% initial state
x0 = [0.05 0 0.1 0 0 0  0.05 0 0.4 0 0 0  0.2 0  0.3 0]';
xhat0 = zeros(1,16);

% desired reference (x,y,z,yaw)
r = [zeros(1,T/2) ones(1,T/2);      % x reference
     zeros(1,T/2) ones(1,T/2);      % y reference
     zeros(1,T);                    % z reference
     zeros(1,T)];                   % yaw reference

% B_ref relates reference to states x_ref = B_ref*r
B_ref = zeros(16,4);
B_ref(3,1) = 1;
B_ref(9,2) = -1;
B_ref(13,3) = 1;
B_ref(15,4) = 1;

x = zeros(length(A(:,1)),T);    % state trajectory
yplot = zeros(length(A(:,1)),T);% output to plot
xhat = zeros(length(A(:,1)),T); % estimated trajectories 
u = zeros(length(B(1,:)),T);    % control inputs
y = zeros(length(C(:,1)),T);    % measurements 
yhat = zeros(length(C(:,1)),T); % estimated output
e = zeros(length(A(:,1)),T);    % observer error
t = zeros(1,T);                 % time vector

Vf = zeros(1,T);                % terminal cost sequence
l = zeros(1,T);                 % stage cost sequence

x(:,1) = x0';

% Define MPC Control Problem

% MPC cost function
%          N-1
% V(u_N) = Sum 1/2[ x(k)'Qx(k) + u(k)'Ru(k) ] + x(N)'Sx(N) 
%          k = 0

% tuning weights
Q = 10*eye(size(A));            % state cost
R = 0.1*eye(length(B(1,:)));    % input cost

% terminal cost = unconstrained optimal cost (Lec 5 pg 6)
[S,~,~] = dare(A,B,Q,R);        % terminal cost % OLD: S = 10*eye(size(A));

% prediction horizon
N = 18; 

Qbar = kron(Q,eye(N));
Rbar = kron(R,eye(N));
Sbar = S;

LTI.A = A;
LTI.B = B;
LTI.C = C;

dim.N = N;
dim.nx = size(A,1);
dim.nu = size(B,2);
dim.ny = size(C,1);

[P,Z,W] = predmodgen(LTI,dim);
              
H = (Z'*Qbar*Z + Rbar + 2*W'*Sbar*W);
d = (x0'*P'*Qbar*Z + 2*x0'*(A^N)'*Sbar*W)';
 
%%

u_limit = 0.1;

Q_kf = 0.01*eye(8);
R_kf = 10*eye(16);

[~,Obs_eigvals,Obs_gain] = dare(A',C',R_kf,Q_kf);
Obs_gain = Obs_gain';

test = A-Obs_gain*C;

% measurement log
% 1
% 2
% 3
% 4
% 5
% 6
% 7
% 8
% 9
% 10
% 11 - OK
% 12 - OK
% 13 - OK
% 14 - OK
% 15
% 16 - 

for k = 1:1:T
    t(k) = (k-1)*h;
    
    % determine reference states based on reference input r
    x_ref = B_ref*r(:,k);
    x0_est = xhat(:,k) - x_ref;
    d = (x0_est'*P'*Qbar*Z + 2*x0_est'*(A^N)'*Sbar*W)';
    
    % compute control action
    cvx_begin quiet
        variable u_N(4*N)
        minimize ( (1/2)*quad_form(u_N,H) + d'*u_N )
        u_N >= -u_limit*ones(4*N,1);
        u_N <=  u_limit*ones(4*N,1);
    cvx_end
    
    u(:,k) = u_N(1:4); % MPC control action
    
    % apply control action on real system
    x(:,k+1) = A*x(:,k) + B*u(:,k); % + B_ref*r(:,k);
    y(:,k) = C*x(:,k);
    
    yplot(:,k) = Cplot*x(:,k);
    
    % observer 
    yhat(:,k) = C*xhat(:,k);
    xhat(:,k+1) = A*xhat(:,k) + B*u(:,k) + Obs_gain*(y(:,k)-yhat(:,k));
    
    e(:,k) = x(:,k) - xhat(:,k);
    
    % stability analysis
    Q = 10*eye(16);
    R = 0.1*eye(4);
    
    [X,eigvals,K] = dare(A,B,Q,R);
    Vf(k) = 0.5*x(:,k)'*X*x(:,k);
    l(k) = 0.5*x(:,k)'*Q*x(:,k);
end

% states_trajectory: Nx16 matrix of output trajectory of 16 states
states_trajectory = yplot';

figure(91);
clf;
% plot(t,xhat(:,1:end-1));
hold on;
% plot(t,x(:,1:end-1));
plot(t,e(:,1:end));
legend('xh1','xh2','xh3','xh4','xh5','xh6','xh7','xh8','xh9','xh10','xh11','xh12','xh13','xh14','xh15','xh16'); %,'x1','x2','x3','x4','x5','x6','x7','x8','x9','x10','x11','x12','x13','x14','x15','x16');
grid on;

%% PLOT RESULTS

% show 3D simulation
X = states_trajectory(:,[3 9 13 11 5 15 1 7]);
visualize_quadrotor_trajectory(states_trajectory(:,[3 9 13 11 5 15 1 7]),0.1);

%% Basic Plots
% plot 2D results fo state trajectories
plot_2D_plots(t, states_trajectory);

% plot the inputs
plot_inputs(t,u,u_limit);
