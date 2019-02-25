function sysc = init_system_dynamics(g,m,L,l,I_xx,I_yy,I_zz)
    
    % SUBSYSTEM 1 - pitch angle dynamics (around y-axis)
    Ac1 = [  0    1 0 0  0 0 ;      % - r1:     the displacement of
             g/L  0 0 0 -g 0 ;      % - r2:     velocity of pendulum relative to quadrotor
             0    0 0 1  0 0 ;      % - x1:     x-direction displacement of quadrotor relative to inertial coordinate frame O
             0    0 0 0  g 0 ;      % - x2:     x-direction velocity of quadrotor relative to inertial coordinate frame O
             0    0 0 0  0 1 ;      % - beta1:  pitch angle of quad (rotation angle around y-axis)
             0    0 0 0  0 0 ];     % - beta2:  pitch angle rate of quad (rotation angle around y-axis)

    Bc1 = [0         0;
           0         0;
           0         0;  
           0         0;  
           0         0;
          -l/I_yy  l/I_yy];

    Cc1 = eye(size(Ac1));

    % SUBSYSTEM 2 - roll angle dynamics (around x-axis)
    Ac2 = [   0   1 0 0  0  0 ;     % - s1:     the displacement of
             g/L  0 0 0 -g  0 ;     % - s2:     velocity of pendulum relative to quadrotor
             0    0 0 1  0  0 ;     % - y1:     x-direction displacement of quadrotor relative to inertial coordinate frame O
             0    0 0 0 -g  0 ;     % - y2:     x-direction velocity of quadrotor relative to inertial coordinate frame O
             0    0 0 0  0  1 ;     % - gamma1: pitch angle of quad (rotation angle around x-axis)
             0    0 0 0  0  0 ];    % - gamma2: pitch angle rate of quad (rotation angle around x-axis)

    Bc2 = [0         0;
           0         0;
           0         0;  
           0         0;  
           0         0;
          -l/I_xx  l/I_xx];

    Cc2 = -eye(size(Ac2)); % negative to have positive y-direction

    % SUBSYSTEM 3 - verticle translational dynamics
    Ac3 = [0 1;
           0 0];

    Bc3 = [0   0   0   0   ;
           1/m 1/m 1/m 1/m];

    Cc3 = eye(size(Ac3));

    % SUBSYSTEM 4 - yaw angle of quadrotor
    Ac4 = [0 1;
           0 0];
    Bc4 = [0        0         0       0   ;
           -1/I_zz  -1/I_zz  1/I_zz 1/I_zz];
    Cc4 = eye(size(Ac4));

    % FULL SYSTEM CONCATENATION
    Ac = blkdiag(Ac1,Ac2,Ac3,Ac4);

    Bc = [Bc1        zeros(6,2);
          zeros(6,2) Bc2;
          Bc3 ;
          Bc4 ];

    Cc = blkdiag(Cc1,Cc2,Cc3,Cc4);

    % continuous-time state-space model
    sysc = ss(Ac,Bc,Cc,[]);

end