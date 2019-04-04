%%
clear
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
LTI.A = A_K;  % closed loop LQR system
LTI.K = K_LQR;          % input u = Kx

% Definition of system dimension
% dim.nx = 2;         % state dimension
% dim.ny = 1;         % output dimension
dim.nx = 16;
dim.nu = 4;

% Definition of the maximal output admissible set
% Output constraints
% F=[1; -1];
% f=ones(2,1);

u_limit = 0.1;

F = kron(eye(dim.nu),[1; -1]);
f = u_limit*ones(2*dim.nu,1);

s = size(F,1);

%% TEST PEDRO CODE

% [H_pedro,h_pedro] = max_out_set_pedro(A_K,K_LQR,0.1);

%%
% Algorithm implementation
exit_flag=0;
k=0;

while exit_flag==0 
    % Set the constraints for the optimization problem
     A=[];
     for t=1:k+1
         for j=1:s
            A_aux(j,:)=F(j,:)*K_LQR*A_K^(t-1);
         end
         A=[A; A_aux];
         clear A_aux
     end
     b = repmat(f,k+1,1);
     
    for i=1:s
         % Set the optimization problem: objective function and constraints
         h = F(i,:)*K_LQR*A_K^(k+1);
         
         % Solve the optimization problem with CVX
         cvx_precision best
         cvx_begin quiet
            variable x_opt(dim.nx)
            maximize(h*x_opt-f(i))
            subject to
            A*x_opt-b<=0;
         cvx_end
         disp('hello!');
        
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


%% CHECK CONSTRAINT

% initial condition
%     r    r   x   x   b b     s    s y y g g  z z  yaw y
x0 = [0.02 0   0.5 0   0 0     0.01 0 0 0 0 0  0 0  0 0]';

% check if the x-location is within X_f
inSet = all(H*x0 <= h);

Xf_set_H = H;
Xf_set_h = h;

%% TEST IF MPC LAW GUIDES TOWARDS X_f IN N-STEPS



N = 18;

dim.nx = 16;
dim.nu = 4;

% initialize MPC matrices
x = zeros(dim.nx,N+1);      % state trajectory

% initial state
x0 = [0.05 0 0.1 0 0 0  0.05 0 0.4 0 0 0  0.2 0  0.3 0]';
x(:,1) = x0;

% tuning weights
Q = 10*eye(dim.nx);         % state cost
R = 0.1*eye(dim.nu);        % input cost

% terminal cost = unconstrained optimal cost (Lec 5 pg 6)
[S,~,~] = dare(Ad,Bd,Q,R);  % terminal cost % OLD: S = 10*eye(size(A));

% prediction horizon
N = 18; 

Qbar = kron(Q,eye(N));
Rbar = kron(R,eye(N));
Sbar = S;

LTI.A = Ad;
LTI.B = Bd;
LTI.C = Cd;

dim.N = N;
dim.nx = size(Ad,1);
dim.nu = size(Bd,2);
dim.ny = size(Cd,1);

[P,Z,W] = predmodgen(LTI,dim);
              
H = (Z'*Qbar*Z + Rbar + 2*W'*Sbar*W);
d = (x0'*P'*Qbar*Z + 2*x0'*(Ad^N)'*Sbar*W)';

% Aeq = [eye(72);-eye(72)];
% beq = ones(72*2,1);
lb = -u_limit*ones(72,1);
ub = u_limit*ones(72,1);

% Define IC_test_vals as the set of initial conditions to consider

r = 1;
s = 1;

res = 31;
bnd_x = 3;
bnd_y = 2;

mat = zeros(res,res);

for initial_r = linspace(-bnd_x,bnd_x,res)
    s = 1;
    
    for initial_s = linspace(-bnd_y,bnd_y,res)

        x(:,1) = [initial_r 0 0 0 0 0  initial_s 0 0 0 0 0  0 0  0 0]';

        for k = 1:N

            x_current = x(:,k);
            d = (x_current'*P'*Qbar*Z + 2*x_current'*(Ad^N)'*Sbar*W)';

            opts = optimoptions('quadprog','Display','off');
            uopt = quadprog(H,d,[],[],[],[],lb,ub,[],opts);

            u(:,k) = uopt(1:4);

            % apply control action
            x(:,k+1) = Ad*x(:,k) + Bd*u(:,k);
        %     y(:,k) = Cd*x(:,k);

        end

        disp('Does x(N) belong to X_f?:');
        disp(initial_r);
        disp(initial_s);
        inSet = all(Xf_set_H*x(:,N+1) <= Xf_set_h);
        disp(inSet);
        
        mat(r,s) = inSet;
        
        s = s+1;
    end
    r = r+1;
end

figure(12);
imagesc(mat);

%% PLOT THE CONSTRAINT SET X_f

H_xyz = Xf_set_H(:,[3,9,13]);


figure(3);
plot3(H_xyz(:,1),H_xyz(:,2),H_xyz(:,3),'*m');
grid();


[K,V] = convhull(H_xyz(:,1),H_xyz(:,2))
hold on
plot(H_xyz(K,1),H_xyz(K,2),'r-');


% xx = -1:.05:1;
% yy = abs(sqrt(xx));
% [x,y] = pol2cart(xx,yy);
% k = convhull(x,y,x+1);
% figure(9);
% plot(x(k),y(k),'r-',x,y,'b*')













