%% QUADROTOR BALANCING PENDULUM MODEL PREDICTIVE CONTROL SIMULATION
%
% MATLAB simulation of the paper A Flying Inverted Pendulum by Markus Hehn 
% and Raffaello D'Andrea 

%% INITIALIZATION 

% Constants
g = 9.81;   % m/s^2
L = 0.6;    % meters
%123
A = [0 1 0 0 0;
     g/L 0 0 0 -g;
     0   0 0 1 0;
     0   0 0 0 g;
     0   0 0 0 0];
 
B = [0;0;0;0;1];
 
%% CHECK CONTROLLABILITY
 
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

%% LINEAR QUADRATIC REGULATOR

% we have three uncontrollable modes lambda = 0,0,0

Q = eye(size(A));
R = 0.1;
Ts = 0.1
C = eye(size(A));
SYSC = ss(A,B,C,0)
SYSD = c2d(SYSC,Ts,'zoh')
Ad = SYSD.A 
Bd = SYSD.B 
[K,S,e] = dlqr(Ad,Bd,Q,R,[]);
iterations = 1000;
x = zeros(5,length(iterations))
x0 = [1 1 0 0 1]';
x(:,1) = x0
for k = 1:iterations
    x(:,k+1) = Ad*x(:,k)-Bd*K*x(:,k);
end

stairs(linspace(1,iterations+1,iterations+1)*Ts,x(1,:))
