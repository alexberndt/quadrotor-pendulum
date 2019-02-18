%% Basic script for quadrotor with pendulum

% Constants
g = 9.81;   % m/s^2
L = 0.6;    % meters


% y rotation

A = [0   1 0 0 0;
     g/L 0 0 0 g;
     0   0 0 0 1;
     0   0 0 0 g;
     0   0 0 0 0];
 
 B = [0;0;0;0;1];
 
 ctrb(A,B)
 
 % not controllable
 % use Hautus test 
 
 eigvals = eig(A);
 
 for lambda = eigvals
     disp('Eigenvalue: ');
     disp(lambda)
    rank(eye(5)*lambda 
 end