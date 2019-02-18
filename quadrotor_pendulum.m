%% QUADROTOR BALANCING PENDULUM MODEL PREDICTIVE CONTROL SIMULATION
%
% MATLAB simulation of the paper A Flying Inverted Pendulum by Markus Hehn 
% and Raffaello D'Andrea 

%% INITIALIZATION 

% Constants
g = 9.81;   % m/s^2
L = 0.6;    % meters

A = [0   1 0 0 0;
     g/L 0 0 0 -g;
     0   0 0 1 0;
     0   0 0 0 g;
     0   0 0 0 0];
 
B = [0;0;0;0;1];
 
%% CHECK CONTROLLABILITY OF CONTINUOUS TIME SYSTEM
 
ctrb(A,B)

% not controllable
% use Hautus test 

eigvals = eig(A);

for idx = 1:numel(eigvals)
 lambda = eigvals(idx);
 disp('Eigenvalue: ');
 disp(lambda)
 rk = rank([(eye(size(A))*lambda-A) B]);
 disp('rank: ');
 disp(rk);
 disp('----------------');
end

%% DISCRETIZE SYSTEM



%% LINEAR QUADRATIC REGULATOR

% we have three uncontrollable modes lambda = 0,0,0

Q = eye(size(A));
R = 0.1;

[K,S,e] = lqr(A,B,Q,R,[]);

% SIMUALTION

C = eye(size(A));
% outputs we are interested are
% - r1:     the displacement of 
% - r2:     velocity of pendulum relative to quadrotor
% - x1:     x-direction displacement of quadrotor relative to inertial
%           coordinate frame O
% - x2:     x-direction velocity of quadrotor relative to inertial
%           coordinate frame O
% - beta:   pitch angle of quad (rotation angle around y-axis)

sys = ss(A-B*K,B,C,[]);

figure(1);
clf;
step(sys);
grid();
