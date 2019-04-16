%% QUADROTOR BALANCING PENDULUM LQR SIMULATION - ROTATIONAL EQUILIBRIUM
%
% MATLAB simulation of the paper A Flying Inverted Pendulum by Markus Hehn 
% and Raffaello D'Andrea 
%
% Equilibrium point rotating around a constant radius at constant rot. vel

%% INIT
clc
clear
addpath('functions/');

disp('------------------------------------------------------------------');
disp('                        ADAPTIVE MPC  ');
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

simTime = 6;   % 2 second simulation
h = 0.1;       % sampling time
fprintf('Sampling time: %f \n',h');
N = simTime/h;        

% define desired set point sequences in terms of Radius and Omega
R_radius_sequence = 2*ones(1,N+1);                  % meters (radius of turn)
Omega_sequence = 1*ones(1,N+1) + 0.007*(1:N+1);    % rad/s (rotational velocity)
OmegaAngle_prev = 0;

%% INITIALIZE VARIABLE SEQUENCES

states_cart = zeros(3,N+1);
states_pendulum = zeros(2,N+1);
states_pendulum_actual = zeros(2,N+1);
euler_angles = zeros(2,N+1);
yaw_angle = zeros(1,N+1);

beta_angle = zeros(1,N+1);
gamma_angle = zeros(1,N+1);
beta_dot_angle = zeros(1,N+1);
gamma_dot_angle = zeros(1,N+1);

reference_trajectory = zeros(2,N+1);

u_tilde = zeros(1,N+1);
v_tilde = zeros(1,N+1);
w_tilde = zeros(1,N+1);

u_actual = zeros(1,N+1);

p_tilde = zeros(1,N+1);
q_tilde = zeros(1,N+1);

p_actual = zeros(1,N+1);
q_actual = zeros(1,N+1);

mu_tilde = zeros(1,N+1);
nu_tilde = zeros(1,N+1);

mu_actual = zeros(1,N+1);
nu_actual = zeros(1,N+1);

x = zeros(12,N+1); 
u = zeros(3, N+1);
y = zeros(12, N);
t = zeros(1,N+1);

Vf = zeros(1,N+1);
P = eye(12);

% define initial conditions
%         p1  p2  q1 q2  u1 u2  v1 v2  w1 w2   mu   nu 
x(:,1) = [0.1  0  0  0   0  0   0  0   0  0   0.10  0.1 ]';

%% ITERATE THROUGH EACH TIME STEP

disp('------------------------------------------------------------------');
disp('                     Simulating MPC System');
disp('------------------------------------------------------------------');
disp('');
fprintf('Simulation time: %d seconds \n',simTime);
disp('');

u_limit = 1;

for k = 1:N
    t(:,k+1) = k*h;
    if ( mod(t(k),1) == 0 )
        fprintf('t = %d sec \n', t(k));
    end
    
    % Obtain new linearized system
    Omega = Omega_sequence(k);
    R_radius = R_radius_sequence(k);
    
    % get new linearized system matrices and Equilibrium nominal points
    [Ad,Bd,Cd,EP] = linearized_ss(L,g,Omega,R_radius,h);
    
    % Define Q and R and determine optimal LQR gain based on new system matrices
    R = 0.1*eye(3);                 
    Q = 1*eye(12);
    Q(5,5) = 10; Q(7,7) = 10; Q(9,9) = 10;
    
    % [X,eigvals,K] = dare(Ad,Bd,Q,R);
    
    % Determine control action from LQR   
    % u(:,k) = -K*x(:,k);
    
    % Define Prediction Matrices
    
    % tuning weights
    Q = 1*eye(12);
    Q(5,5) = 10; Q(7,7) = 10; Q(9,9) = 10;  % state cost
    R = 0.1*eye(length(Bd(1,:)));           % input cost
    S = 10*eye(size(Ad));                   % terminal cost
    
    [X,eigvals,K] = dare(Ad,Bd,Q,R);
    
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
    
    x0 = x(:,k);
    d = (x0'*P'*Qbar*Z + 2*x0'*(Ad^N)'*Sbar*W)';
    
    limits = [5*ones(N,1);
              1.5*ones(N,1);
              0.5*ones(N,1)];  
    
    % MPC control action
    cvx_begin quiet
        variable u_N(3*N)
        minimize ( (1/2)*quad_form(u_N,H) + d'*u_N )
        u_N >= -limits;
        u_N <=  limits;
    cvx_end
    
    u(:,k) = u_N(1:3); % MPC control action
    
    % Save variables
    u_tilde(k) = x(5,k);
    v_tilde(k) = x(7,k);
    w_tilde(k) = x(9,k);

    u_actual(k) = u_tilde(k) + EP.u_0;

    p_tilde(k) = x(1,k);
    q_tilde(k) = x(3,k);

    p_actual(k) = p_tilde(k) + EP.p_0;
    q_actual(k) = q_tilde(k) + EP.q_0;

    mu_tilde(k) = x(11,k);
    nu_tilde(k) = x(12,k);

    mu_actual(k) = mu_tilde(k) + EP.mu_0;
    nu_actual(k) = nu_tilde(k) + EP.nu_0;

    OmegaAngleT = OmegaAngle_prev + Omega*h;
    OmegaAngle_prev = OmegaAngleT;
    
    R_uvw_to_xyz = [cos(OmegaAngleT) -sin(OmegaAngleT) 0;
                    sin(OmegaAngleT)  cos(OmegaAngleT) 0;
                     0             0         1];
                 
    % determine x,y,z states of quadrotor
    states_cart(:,k) = R_uvw_to_xyz*[u_actual(k); v_tilde(k); w_tilde(k)];
    
    R_pq_to_rs = [cos(OmegaAngleT)  -sin(OmegaAngleT);
                  sin(OmegaAngleT)  cos(OmegaAngleT)];  
    R_euler = R_pq_to_rs;          
              
    % determine 
    euler_angles(:,k) = R_euler*[nu_actual(k); mu_actual(k)];
    
    yaw_angle(k) = 0;
              
    % determine r,s states of pendulum (relative to quadrotor center)
    states_pendulum(:,k) = R_pq_to_rs*[p_actual(k); q_actual(k)]; 
    
    % determine absolute pendulum position 
    states_pendulum_actual(:,k) = states_pendulum(:,k) + states_cart(1:2,k);
    
    % calculate control input in terms of omega_x,y,z
    mu = mu_actual(k);
    nu = nu_actual(k);
    
    % determine gamma, beta angles of quadrotor
    gamma = -asin(sin(OmegaAngleT)*sin(mu)*cos(nu) - cos(OmegaAngleT)*sin(nu));
    beta = asin((cos(OmegaAngleT)*sin(mu)*cos(nu) + sin(OmegaAngleT)*sin(nu))/(cos(gamma)));
    
    beta_angle(k) = beta;
    gamma_angle(k) = gamma;
    
    beta_dot_angle(k) = (R_radius*Omega^3*acos(gamma)*(tan(beta)*tan(gamma)*cos(OmegaAngleT) + acos(beta)*sin(OmegaAngleT)))/sqrt(g^2+(R_radius*Omega^2)^2);
    gamma_dot_angle(k) = (R_radius*Omega^3*acos(gamma)*cos(OmegaAngleT))/sqrt(g^2+(R_radius*Omega^2)^2);
    
    % determine reference input based on R_radius and Omega sequences
    reference_trajectory(:,k) = [cos(OmegaAngleT); sin(OmegaAngleT)]*R_radius;
    
    % Apply control and update state equations
    x(:,k+1) = Ad*x(:,k) + Bd*u(:,k);
    y(:,k) = Cd*x(:,k);
    
    % Determine Terminal Constraint
    Vf(k) = 0.5*x(:,k)'*X*x(:,k);
end

%% PLOT 3D TRAJECTORY

states_trajectory = [states_cart(1:3,:);
                     euler_angles(1,:);
                     euler_angles(2,:);   
                     yaw_angle;
                     states_pendulum(1:2,:)]';

visualize_quadrotor_trajectory_rotating(states_trajectory, reference_trajectory);
     
%% PLOT RESULTS

other.mu_actual = mu_actual;
other.beta_dot_angle = beta_dot_angle;
other.gamma_dot_angle = gamma_dot_angle;
other.beta_angle = beta_angle;
other.gamma_angle = gamma_angle;

plot_rotational(t,x,u,other);

% stability plots
% figure(123);
% Vf_dif = Vf(2:end) - Vf(1:end-1);
% stairs(Vf_dif);
% hold on
% stairs(Vf);
% grid();

%% FUNCTIONS
function Rt = R_x(y)
   Rt = [1   0      0    ;
         0 cos(y) -sin(y);
         0 sin(y)  cos(y)];
end

function Rt = R_y(y)
    Rt = [cos(y) 0  sin(y)    ;
            0    1    0;
         -sin(y) 0  cos(y)];
end
     
function Rt = R_z(y)
    Rt = [cos(y) -sin(y) 0;
          sin(y)  cos(y) 0;
            0       0    1];
end
     
function [Ad,Bd,Cd,EP] = linearized_ss(L,g,Omega,R_radius,h)

    % equilibrium constants
    q_0 = 0;
    zeta_0 = 0.5;   % initial zeta_0 estimate 
    % zeta_0 = sqrt(L^2-q_0^2-p_0^2); % zeta = relative position of pendulum along z-axis (Eqn 8)
    p_0 = -(Omega^2*R_radius)/(Omega^2+g/zeta_0);
    zeta_0 = sqrt(L^2-q_0^2-p_0^2); % correct

    z_0 = Inf; % NOT USED - constant reference altitude

    u_0 = R_radius;
    v_0 = 0;
    w_0 = 0;

    mu_0 = atan(-Omega^2*R_radius/g);        % nominal euler angle
    nu_0 = 0;                           % nominal euler angle
    a_0 = sqrt(g^2 + (R_radius*Omega^2)^2);    % nominal thrust

    C1 = zeta_0^2/L^2;
    C2 = Omega^2 + g*L^2/(zeta_0^3);
    C3 = Omega;
    C4 = -(p_0/zeta_0)*a_0*sin(mu_0) - a_0*cos(mu_0);
    C5 = (p_0/zeta_0)*cos(mu_0) - sin(mu_0);

    C6 = Omega^2 + g/zeta_0;

    %      p1  p2  q1  q2    u1  u2  v1   v2  w1  w2  mu            nu  
    Ac =[  0   1   0   0     0   0   0    0   0   0   0             0   ;  % p1
         C1*C2 0   0 2*C1*C3 0   0   0    0   0   0 C1*C4           0   ;  % p2
           0   0   0   1     0   0   0    0   0   0   0             0   ;  % q1
           0 -2*C3 C6  0     0   0   0    0   0   0   0            a_0  ;  % q2
           0   0   0   0     0   1   0    0   0   0   0             0   ;  % u1  
           0   0   0   0   C3^2  0   0   2*C3 0   0 a_0*cos(mu_0)   0   ;  % u2
           0   0   0   0     0   0   0    1   0   0   0             0   ;  % v1
           0   0   0   0     0 -2*C3 C3^2 0   0   0   0           -a_0  ;  % v2
           0   0   0   0     0   0   0    0   0   1   0             0   ;  % w1
           0   0   0   0     0   0   0    0   0   0 -a_0*sin(mu_0)  0   ;  % w2
           0   0   0   0     0   0   0    0   0   0   0             0   ;  % mu
           0   0   0   0     0   0   0    0   0   0   0             0   ]; % nu     

    %      mu_dot nu_dot a
    Bc =[  0      0      0;     % p1
           0      0    C1*C5;   % p2
           0      0      0;     % q1
           0      0      0;     % q2
           0      0      0;     % u1
           0      0  sin(mu_0); % u2
           0      0      0;     % v1
           0      0      0;     % v2
           0      0      0;     % w1
           0      0  cos(mu_0); % w2
           1      0      0;     % mu
           0      1      0 ];   % nu     
    
    % full state feedback
    Cc = eye(12);
    
    % descritize system
    sysc = ss(Ac,Bc,Cc,[]);
    
    sysd = c2d(sysc,h);

    Ad = sysd.A;
    Bd = sysd.B;
    Cd = sysd.C;
    
    % assign equilibrium parameters to struct
    EP.a_0 = a_0;
    EP.q_0 = q_0;
    EP.p_0 = p_0;
    EP.zeta_0 = zeta_0;
    
    EP.mu_0 = mu_0;
    EP.nu_0 = nu_0;
    
    EP.u_0 = u_0;
    EP.v_0 = v_0;
    EP.w_0 = w_0;
    
    EP.z_0 = z_0;
end
     
     
     
     
     
     
     
%.