function visualize_quadrotor_trajectory(states_trajectory)
    %% VISUALIZE QUADROTOR TRAJECTORY
    %
    % plots the dynamics of the quadrotor with inverted pendulum given the 
    % provided states_trajectory which must be as specified below
    % INPUTS
    % - time:               Nx1 vector of time indices
    %       
    % - states_trajectory:  Nx8 vector of quadrotor and pendulum states
    %       1) x        displacement of COM in x-direction
    %       2) y        displacement of COM in y-direction
    %       3) z        displacement of COM in z-direction
    %       4) roll     rotation of quadrotor around x-axis
    %       5) pitch    rotation of quadrotor around y-axis
    %       6) yaw      rotation of quadrotor around z-axis
    %       7) r        displacement of pendulum COM in x-direction
    %       8) s        displacement of pendulum COM in y-direction
    % OUTPUS
    % - none 
    
    %% INIT
    
    x_r=0.3;
    y_r=0.3;
    z_r=2;
    
    % quadrotor frame and circle drawing
    l=0.3;      
    rc=0.1;
    rx=rc*cos(linspace(0,2*pi,1e1));
    rx=[rx rx(1)];
    ry=rc*sin(linspace(0,2*pi,1e1));
    ry=[ry ry(1)];
    
    X = states_trajectory(:,1:6);
    
    figure()
    clf;
    grid();

    w = 1.5;
    Ax = [-w+x_r w+x_r -w+y_r w+y_r -w+z_r w+z_r];

    [N,~] = size(X);
    
    for j = 1:N
        
        % X is the 6-states of the quadrotor
        
        % obtain relative positions
        Rt = R(X(j,4:6)); 
        R1 = Rt*([ l+rc+rx ; ry  ; zeros(size(rx)) ]) + X(j,1:3)';
        R2 = Rt*([ rx ; l+rc+ry  ; zeros(size(rx)) ]) + X(j,1:3)';
        R3 = Rt*([ -l-rc+rx ; ry ; zeros(size(rx)) ]) + X(j,1:3)';
        R4 = Rt*([ rx ; -l-rc+ry ; zeros(size(rx)) ]) + X(j,1:3)';
        
        % create plotting frames
        A1 = Rt*([-l l;0 0;0 0]) + X(j,1:3)';
        A2 = Rt*([0 0; -l l;0 0]) + X(j,1:3)';
        
        % plot to 3D plane
        plot3(x_r,y_r,z_r,'r*',A1(1,:),A1(2,:),A1(3,:),'k',A2(1,:),A2(2,:),A2(3,:),'k',R1(1,:),R1(2,:),R1(3,:),'r',R2(1,:),R2(2,:),R2(3,:),'b',R3(1,:),R3(2,:),R3(3,:),'b',R4(1,:),R4(2,:),R4(3,:),'b');
    
        
        axis(Ax);
        set(gca,'box','on')
        drawnow
    end
    
    function y = R(Xrot)
        phi = Xrot(1);
        theta = Xrot(2);
        psi = Xrot(3);

        Rpsi = [cos(psi)  -sin(psi)   0;
                sin(psi)  cos(psi)    0;
                0         0           1];

        % rotation around y with theta
        Rtheta = [cos(theta)    0       sin(theta);
                  0             1       0
                 -sin(theta)    0       cos(theta)];

        % rotation around x with phi 
        Rphi = [1       0           0;
                0       cos(phi)    -sin(phi);
                0       sin(phi)    cos(phi)];

        y=Rpsi*Rtheta*Rphi;
    end
    
end