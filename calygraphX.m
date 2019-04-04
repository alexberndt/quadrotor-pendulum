%% Determine Calygraph X Area for Linear State-Space System
%
%
%
%

%% INITIALIZATION
clc
clear
addpath('functions/');

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

%% DISCRETIZE SYSTEM
simTime = 5;
h = 0.1;

sysd = c2d(sysc,h);

Ad = sysd.A;
Bd = sysd.B;
Cd = zeros(8,16);                           % Measured outputs
Cd(1,1) = 1; Cd(2,3) = 1;                   % Pendulum position
Cd(3,5) = 1; Cd(4,7) = 1; Cd(5,9) = 1;      % Quad position
Cd(6,11) = 1; Cd(7,13) = 1; Cd(8,15) = 1;   % Quad rotation angle

%% CHECK CONTROLLABILITY AND OBSERVABILITY

K = ctrb(Ad,Bd);
rank(K);
W = ctrb(Ad',Cd');
rank(W);

%% VARIABLE DEFINITIONS

T = simTime/h;

t = zeros(1,T);         % time sequence
x = zeros(16,T);        % states
y = zeros(8,T);         % output 
u = zeros(4,T);         % inputs

%% DEFINE PREDICTION MATRICES

LTI.A = Ad;
LTI.B = Bd;
LTI.C = eye(16);

dim.N = 10;     % prediction horizon
dim.nx = 16;    
dim.ny = 16;
dim.nu = 4;

[S,Z,W] = predmodgen(LTI,dim);

Q = 10*eye(size(Ad));            % state cost
R = 0.1*eye(length(Bd(1,:)));    % input cost
Qbar = kron(Q,eye(dim.N));
Rbar = kron(R,eye(dim.N));

% terminal cost = unconstrained optimal cost (Lec 5 pg 6)
[P,~,~] = dare(Ad,Bd,Q,R);        % terminal cost 

x0 = zeros(1,16)';

H = (Z'*Qbar*Z + Rbar + 2*W'*P*W);
d = (x0'*S'*Qbar*Z + 2*x0'*(Ad^dim.N)'*P*W)';

%% DETERMINE OUTPUT ADMISSIBLE SET

[X,L,K_LQR] = dare(Ad,Bd,Q,R);

A_K = Ad-Bd*K_LQR;
eigvals_A_K = eig(A_K);
% figure(1);
% plot(real(eigvals_A_K),imag(eigvals_A_K),'*m');
% % xlim([-1 1]);
% % ylim([-1 1]);
% grid on;

% [output] = calcOutputAdmissibleSet(Ad,Bd,eye(16),K_LQR);

%% CALCULATE OUTPUT ADMISSIBLE SET

dim.nx = 16;
dim.nu = 4;

u_limit = 1; 
F = kron(eye(4),kron(ones(dim.nx,1),[1; -1]));
f = u_limit*ones(4*2*dim.nx,1);














































%% DETERMINE CALYGRAPH X

% [output, info] = calygraphX(LTI,dim)
 
%% MODEL PREDICTIVE CONTROL SIMULATION

% for k = 1:T
%     t(k) = (k-1)*h;
%     
%     
%     
%     
% end


%%




















































%.