%% QUADROTOR BALANCING PENDULUM MODEL PREDICTIVE CONTROL SIMULATION
%
% MATLAB simulation of the paper A Flying Inverted Pendulum by Markus Hehn 
% and Raffaello D'Andrea using a Model Predictive Controller
%
% Reference tracking with disturbance rejection using Optimal Target
% Selection (OTS) as in Lecture 6 page 17/21

%% INITIALIZATION
clc
clear
addpath('functions/');

disp('------------------------------------------------------------------');
disp('          OUTPUT MPC WITH DISTURBANCE REJECTION ');
disp('');
disp('------------------------------------------------------------------');


%% DEFINE CONSTANTS
g = 9.81;       % m/s^2
m = 0.5;        % kg
L = 0.565;      % meters (Length of pendulum to center of mass)
l = 0.17;       % meters (Quadrotor center to rotor center)
I_yy = 3.2e-3;  % kg m^2 (Quadrotor inertia around y-axis)
I_xx = I_yy;    % kg m^2 (Quadrotor inertia around x-axis)
I_zz = 5.5e-3;  % kg m^2 (Quadrotor inertia around z-axis)

%% DEFINE STATE SPACE SYSTEM
sysc = init_system_dynamics(g,m,L,l,I_xx,I_yy,I_zz);
check_controllability(sysc);

%% DISCRETIZE SYSTEM
simTime = 15;   % simulation time in seconds
h = 0.10;       % sampling time in seconds

sysd = c2d(sysc,h);
T = simTime/h;

A = sysd.A;
B = sysd.B;
Cplot = sysd.C;                         % use full states for plotting

C = zeros(8,16);                        % Measured outputs
% C(1,1) = 1; C(2,7) = 1;                 % Pendulum position (r=1 s=7)
% C(3,3) = 1; C(4,9) = 1; C(5,13) = 1;    % Quad position (x=3 y=9 z=13)
% C(6,5) = 1; C(7,11) = 1; C(8,15) = 1;   % Quad rotation angle rate (beta=5 gamma=11 yaw=15)

C(1,1) = 1; C(4,7) = 1;                 % Pendulum position
C(2,3) = 1; C(5,9) = 1; C(7,13) = 1;    % Quad position
C(3,6) = 1; C(6,12) = 1; C(8,15) = 1;   % Quad rotation angle rates

Oc = ctrb(A',C');
Obs_rank = rank(Oc);
disp('Rank of observability matrix');
disp(Obs_rank);

%% MODEL PREDICTIVE CONTROL

% actual initial state
%       r   r x    x beta   s   s y    y gamma z    zd  yaw yawd
x0 = [0.005 0 0.01 0 0 0  0.005 0 0.04 0  0 0  0.02 0  0.03 0]';

% observer initial state
xhat0 = zeros(1,16);
 
% desired output reference 
% y_ref_OTS = [zeros(1,T);                    % r 
%              zeros(1,2*T/5) -1*ones(1,3*T/5);   % x
%              zeros(1,T);                    % beta
%              zeros(1,T);                    % s
%              zeros(1,2*T/5)  1*ones(1,3*T/5);   % y
%              zeros(1,T);                    % gamma
%              zeros(1,T);                    % z
%              zeros(1,2*T/5)  1*ones(1,3*T/5)];  % yaw 
%          
y_ref_OTS = [zeros(1,T);                    % r 
             zeros(1,T/3) -1*ones(1,2*T/3);   % x
             zeros(1,T);                    % beta
             zeros(1,T);                    % s
             zeros(1,T/3)  1*ones(1,2*T/3);   % y
             zeros(1,T);                    % gamma
             zeros(1,T);                    % z
             zeros(1,T/3)  1*ones(1,2*T/3)];  % yaw 
          
n_d = 2;

% Optimal Target Selection reference states and inputs
x_ref_OTS = zeros(16,T);
u_ref_OTS = zeros(4,T);
dhat = zeros(n_d,T);
 
n_d = 1;
% disturbance input (x)
d_dist = [zeros(n_d,T/15) 0.01*ones(n_d,T/15) zeros(n_d,8*T/15) 0.01*ones(n_d,5*T/15);
          zeros(n_d,T/15) 0.01*ones(n_d,T/15) zeros(n_d,5*T/15) 0.01*ones(n_d,8*T/15)];  
n_d = 2;
% d_dist = [zeros(n_d,T/5) 0.01*ones(n_d,T/5) zeros(n_d,3*T/5)];

x = zeros(length(A(:,1)),T);        % state trajectory
yplot = zeros(length(A(:,1)),T);    % output to plot
xhat = zeros(length(A(:,1)),T);     % estimated trajectories 
xhaug = zeros(length(A(:,1))+n_d,T);  % augmented states (16 + 2) x + d

u = zeros(length(B(1,:)),T);        % control inputs
y = zeros(length(C(:,1)),T);        % measurements 
yhat = zeros(length(C(:,1)),T);     % estimated output

e = zeros(length(A(:,1)),T);        % observer error
t = zeros(1,T);                     % time vector

Vf = zeros(1,T);                    % terminal cost sequence
l = zeros(1,T);                     % stage cost sequence

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
N = 10; 

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
 
%% Define augmented observer

Bd = [0 0 1 0 0 0  0 0 0 0 0 0  0 0  0 0;
      0 0 0 0 0 0  0 0 1 0 0 0  0 0  0 0]';
Cd = [0   0   0    0   0   1    0    0  ;
      0   0   1    0   0   0    0    0  ]';
  
n_d = 2;

Aaug = [A           Bd;
        zeros(n_d,16) eye(n_d)];

Baug = [B; zeros(n_d,4)];

Caug = [C Cd];

test_Obs = [eye(16)-A -Bd; C Cd];
disp('Rank of Augmented System n+nd = 18 for full rank');
disp(rank(test_Obs))

Q_kf = 1*eye(8);
R_kf = 1*eye(16+n_d);

[~,Obs_eigvals,Obs_gain] = dare(Aaug',Caug',R_kf,Q_kf);
Obs_gain = Obs_gain';
disp('Eigenvalues of Kalman Filter Observer:');
disp(abs(Obs_eigvals))

%% Simulate system

u_limit = 0.1;

disp('------------------------------------------------------------------');
disp('                Simulating Output MPC System');
disp('------------------------------------------------------------------');
disp('');
fprintf('Simulation time: %d seconds\n',simTime);
disp('');

for k = 1:1:T
    t(k) = (k-1)*h;
    if ( mod(t(k),1) == 0 )
        fprintf('t = %d sec \n', t(k));
    end
    
    % Optimal Target Selector
    Q_OTS = eye(16);
    R_OTS = eye(4);
    J_OTS = blkdiag(Q_OTS,R_OTS);
    
    A_OTS = [eye(16)-A     B ;
             C         zeros(8,4)];
    b_OTS = [Bd*dhat(:,k);
             y_ref_OTS(:,k) - Cd*dhat(:,k)];
    
%     cvx_begin quiet
%         variable xr_ur(20)
%         minimize ( quad_form(xr_ur,J_OTS) )
%         A_OTS*xr_ur == b_OTS;
%     cvx_end
    
    opts = optimoptions('quadprog','Display','off');
    [xr_ur,~,exitflag] = quadprog(J_OTS,zeros(20,1),[],[],A_OTS,b_OTS,[],[],[],opts);
    
    x_ref_OTS(:,k) = xr_ur(1:16);
    u_ref_OTS(:,k) = xr_ur(17:20);
    
    % determine reference states based on reference input r
    % x_ref = B_ref*r(:,k);
    x0_est = xhaug(1:16,k) - x_ref_OTS(:,k);
    d = (x0_est'*P'*Qbar*Z + 2*x0_est'*(A^N)'*Sbar*W)';
    
    % compute control action
    cvx_begin quiet
        variable u_N(4*N)
        minimize ( (1/2)*quad_form(u_N,H) + d'*u_N )
        u_N >= -u_limit*ones(4*N,1);
        u_N <=  u_limit*ones(4*N,1);
    cvx_end
    
    u(:,k) = u_N(1:4); % MPC control action
    
    B_dist = Bd;
    
    % apply control action on real system
    x(:,k+1) = A*x(:,k) + B*u(:,k) + B_dist*d_dist(:,k); % + B_ref*r(:,k);
    y(:,k) = C*x(:,k) + Cd*d_dist(:,k);
    
    Cdplot = [0 0 1 0 0 0  0 0 0 0 0 0  0 0  0 0]';
    yplot(:,k) = Cplot*x(:,k) + Cdplot*d_dist(1,k);
    
    % observer 
    % yhat(:,k) = C*xhat(:,k);
    % xhat(:,k+1) = A*xhat(:,k) + B*u(:,k) + Obs_gain*(y(:,k)-yhat(:,k));
    
    % augmented observer
    yhat(:,k) = Caug*xhaug(:,k);
    xhaug(:,k+1) = Aaug*xhaug(:,k) + Baug*u(:,k) + Obs_gain*( y(:,k) - yhat(:,k) );
    
    dhat(:,k+1) = xhaug(17:18,k+1); 
    % e(:,k) = x(:,k) - xhat(:,k);
    
    % stability analysis
    Q = 10*eye(16);
    R = 0.1*eye(4);
    
    [X,eigvals,K] = dare(A,B,Q,R);
    Vf(k) = 0.5*x(:,k)'*X*x(:,k);
    l(k) = 0.5*x(:,k)'*Q*x(:,k);
end

% states_trajectory: Nx16 matrix of output trajectory of 16 states
states_trajectory = yplot';

saved_data.t = t;
saved_data.x = yplot;
saved_data.u = u;

%%
figure(91);
clf;
% plot(t,xhat(:,1:end-1));
hold on;
plot(t,xhaug(:,1:end-1));
% plot(t,e(:,1:end));
legend('xh1','xh2','xh3','xh4','xh5','xh6','xh7','xh8','xh9','xh10','xh11','xh12','xh13','xh14','xh15','xh16','d1','d2'); %,'x1','x2','x3','x4','x5','x6','x7','x8','x9','x10','x11','x12','x13','x14','x15','x16');
grid on;

%%

figure(93);
clf;
% plot(t,xhat(:,1:end-1));
hold on;
plot(t,x_ref_OTS(:,1:end));
ylim([-1.5 2.5]);
% plot(t,e(:,1:end));
legend('r','rd','x','xd','beta','betad','s','sd','y','yd','gamma','gammad','z','zd','yaw','yawd'); %,'x1','x2','x3','x4','x5','x6','x7','x8','x9','x10','x11','x12','x13','x14','x15','x16');
grid on;

%%

figure(92);
clf;
subplot 211
hold on;
plot(t,25*d_dist(1,:),'r-');
plot(t,25*xhaug(17,2:end),'b-');
grid on;
ylim([-0.02 0.28]);
legend('$d_1$','$d_1$ estimated','interpreter','latex','Location','southeast');
subplot 212
hold on;
plot(t,25*d_dist(2,:),'r-');
plot(t,25*xhaug(18,2:end),'b-');
grid on;
ylim([-0.02 0.31]);
legend('$d_2$','$d_2$ estimated','interpreter','latex','Location','southeast');
% title('Estimated disturbance');

%% PLOT RESULTS

% show 3D simulation
X = states_trajectory(:,[3 9 13 11 5 15 1 7]);
visualize_quadrotor_trajectory(X);

%% Basic Plots
% plot 2D results fo state trajectories

plot_2D_plots_offset(t, states_trajectory, d_dist, x_ref_OTS(3,:)); %547

% plot the inputs
plot_inputs(t,u,u_limit);
