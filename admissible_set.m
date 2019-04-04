%%
clear all
close all
clc
addpath('functions/');

%% DEFINE CONSTANTS
g = 9.81;       % m/s^2
m = 0.5;        % kg
L = 0.565;      % meters (Length of pendulum to center of mass)
l = 0.17;       % meters (Quadrotor center to rotor center)
I_yy = 3.2e-3;  % kg m^2 (Quadrotor inertia around y-axis)
I_xx = I_yy;    % kg m^2 (Quadrotor inertia around x-axis)
I_zz = 5.5e-3;  % kg m^2 (Quadrotor inertia around z-axis)

%% INITIALIZATION
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

%% MPC COST FUNCTION VALUES

Q = 10*eye(size(Ad));            % state cost
R = 0.1*eye(length(Bd(1,:)));    % input cost

[X,L,K_LQR] = dare(Ad,Bd,Q,R);

A_K = Ad-Bd*K_LQR;
eigvals_A_K = eig(A_K);

%% THEIR CODE

% Definition of the LTI system
LTI.A = [.9 1; 0 .09];  % closed loop LQR system
LTI.C = [1 1];          % input u = Kx

% Definition of system dimension
dim.nx = 2;         % state dimension
dim.ny = 1;         % output dimension

% Definition of the maximal output admissible set
% Output constraints
F=[1; -1];
f=ones(2,1);

s=size(F,1);

%%
% Algorithm implementation
exit_flag=0;
k=0;

while exit_flag==0 
    % Set the constraints for the optimization problem
     A=[];
     for t=1:k+1
         for j=1:s
            A_aux(j,:)=F(j)*LTI.C*LTI.A^(t-1);
         end
         A=[A; A_aux];
         clear A_aux
     end
     b=repmat(f,k+1,1);
     
    for i=1:s
         % Set the optimization problem: objective function and constraints
         h=F(i)*LTI.C*LTI.A^(k+1);
         
         % Solve the optimization problem with CVX
         cvx_precision best
         cvx_begin quiet
            variable x_opt(dim.nx)
            maximize(h*x_opt-f(i))
            subject to
            A*x_opt-b<=0;
        cvx_end
        
        % Save optimal value
        opt_val(i)=cvx_optval;
    end
    
    % Evaluating optimality condition
    if (sum(opt_val<=0-eps)==s && strcmp(cvx_status,'Unbounded')==0)
        exit_flag=1;
        k_star=k
        H=A
        h=b
    else
        clear opt_val A_aux A b h
        k=k+1;
    end
end






















