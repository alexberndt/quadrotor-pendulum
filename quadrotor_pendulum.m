%% Basic script for quadrotor with pendulum

% Constants
g = 9.81;   % m/s^2
L = 0.6;    % meters

% y rotation

A = [0   1 0 0 0;
     g/L 0 0 0 -g;
     0   0 0 1 0;
     0   0 0 0 g;
     0   0 0 0 0];
 
 B = [0;0;0;0;1];
 
%  A = [0 1 0; 0 0 1; 0 0 0];
%  B = [0;0;1];
 
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
 
 % we have three uncontrollable modes lambda = 0,0,0
 
 Q = eye(size(A));
 R = 1;
 
 [K,S,e] = lqr(A,B,Q,R,[]);