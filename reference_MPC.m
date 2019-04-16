%% QUADROTOR BALANCING PENDULUM MODEL PREDICTIVE CONTROL SIMULATION
%
% MATLAB simulation of the paper A Flying Inverted Pendulum by Markus Hehn 
% and Raffaello D'Andrea using a Model Predictive Controller

%% INIT
clc
clear
addpath('functions/');

disp('------------------------------------------------------------------');
disp('                        REGULATION MPC  ');
disp('');
disp('------------------------------------------------------------------');


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
simTime = 6;
h = 0.1;
fprintf('Sampling time: %f \n',h');

sysd = c2d(sysc,h);
T = simTime/h;

A = sysd.A;
B = sysd.B;
C = sysd.C;

%% MODEL PREDICTIVE CONTROL

% initial state
%     r1   r2 x1 x2 b1 b2    s1 s2 y1 y2 g1 g2   z1 z2   yaw1 yaw2
x0 = [0.01 0 0.2 0 0.1 0       0.05 0 0.12 0 0 0    0.1 0   0.12 0]';
% x0 = [0.01 0 0.01 0 0 0  0 0 0.01 0 0 0  0 0  0 0]';

% desired reference (x,y,z,yaw)
r = [zeros(1,T);     % x reference
     zeros(1,T);     % y reference
     zeros(1,T);     % z reference
     zeros(1,T)];    % yaw reference

% B_ref relates reference to states x_ref = B_ref*r
B_ref = zeros(16,4);
B_ref(3,1) = 1;
B_ref(9,2) = -1;
B_ref(13,3) = 1;
B_ref(15,4) = 1;

x = zeros(length(A(:,1)),T);    % state trajectory
u = zeros(length(B(1,:)),T);    % control inputs
y = zeros(length(C(:,1)),T);    % measurements 
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

% S = Q;
% S = 1000*eye(16);

% prediction horizon
N = 20; 
fprintf('Prediction horizon %d \n',N);

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
disp('------------------------------------------------------------------');
disp('                     Simulating MPC System');
disp('------------------------------------------------------------------');
disp('');
fprintf('Simulation time: %d seconds \n',simTime);
disp('');

u_limit = 0.1;
%            r   r    x x b b  s   s    y y g g  z z  y y
x_lim_vec = [0.2 0.5 5 5 5 5  0.2 0.5 5 5 5 5  5 5  5 5]';
x_lim_vec_full = repmat(x_lim_vec,[N 1]);

for k = 1:1:T
    t(k) = (k-1)*h;
    if ( mod(t(k),1) == 0 )
        fprintf('t = %d sec \n', t(k));
    end
    
    % determine reference states based on reference input r
    x_ref = B_ref*r(:,k);
    x0 = x(:,k) - x_ref;
    d = (x0'*P'*Qbar*Z + 2*x0'*(A^N)'*Sbar*W)';
    
    % compute control action
    cvx_begin quiet
        variable u_N(4*N)
        minimize ( (1/2)*quad_form(u_N,H) + d'*u_N )
        % input constraints
        u_N <=  u_limit*ones(4*N,1);
        u_N >= -u_limit*ones(4*N,1);
        % state constraints
        Z*u_N <= -P*x0 + x_lim_vec_full;
        Z*u_N >= -P*x0 - x_lim_vec_full; 
    cvx_end
    
    u(:,k) = u_N(1:4); % MPC control action
    
    % apply control action
    x(:,k+1) = A*x(:,k) + B*u(:,k); % + B_ref*r(:,k);
    y(:,k) = C*x(:,k);
    
    % stability analysis
%     Q = 10*eye(16);
%     R = 0.1*eye(4);
    
    [X,eigvals,K] = dare(A,B,Q,R);
    Vf(k) = 0.5*x(:,k)'*X*x(:,k);
    l(k) = 0.5*x(:,k)'*Q*x(:,k) + 0.1*u(:,k)'*R*u(:,k);
end

% states_trajectory: Nx16 matrix of trajectory of 16 states
states_trajectory = y';

%% PLOT RESULTS

% Stability plots
Vf_diff = Vf(2:end)-Vf(1:end-1);

Vfkp1 = Vf(2:end);
Vfk = Vf(1:end-1);

lQ = l(2:end);
% lQ = l(1:end-1);

kt = t(1:end-1);

% figure(123);
% clf;
% hold on;
% % stairs(kt,Vfkp1);
% % stairs(kt,Vfk);
% % stairs(kt,lQ);
% stairs(kt,Vf_diff);
% stairs(kt,-0.3*lQ);
% % stairs(kt,Vfk-lQ);
% % stairs(kt,Vfkp1-(Vfk-lQ));
% grid on;

% legend('$V_f(f(x,u)) - V_f(x)$','$-l(x,u)$','interpreter','latex');

% legend('Vf','Vf(k+1)-Vf(k)','stage l(k)','Vf - l');

% show 3D simulation
X = states_trajectory(:,[3 9 13 11 5 15 1 7]);
visualize_quadrotor_trajectory(states_trajectory(:,[3 9 13 11 5 15 1 7]));

saved_data.t = t;
saved_data.x = states_trajectory;
saved_data.u = u;

%% Basic Plots
% plot 2D results fo state trajectories
plot_2D_plots(t, states_trajectory);
% 
% % plot the inputs
plot_inputs(t,u,u_limit);

%% Comparison Plots

plot_comparison_S_different(); % 543
plot_comparison_R_different(); % 544
plot_comparison_Q_different(); % 546

plot_comparison_R_inputs();    % 589

%%
plot_comparison_MPC_LQR();     % 567
% 
% plot_comparison_horizon();      % 789
