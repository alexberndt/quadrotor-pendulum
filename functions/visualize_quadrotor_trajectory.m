function visualize_quadrotor_trajectory(states_trajectory,pause_duration)
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
    if (nargin == 1)
        pause_duration = 0;
    end
    
    % X is the 6-states of the quadrotor, 2-states of pendulum
    X = states_trajectory(:,1:8);
    [N,~] = size(X);

    % reference set point (center of plot)
    x_r = 0;
    y_r = 0;
    z_r = 0;

    % quadrotor frame and circle drawings
    l =  0.34;   
    pl = 0.756; % pendulum height
    rc = 0.1;
    rx = rc*cos(linspace(0,2*pi,20));
    rx = [rx rx(1)];
    ry = rc*sin(linspace(0,2*pi,20));
    ry = [ry ry(1)];
    
    % init figure
    figure(42);
    clf;
    
    % define plot axes limits
    w = 1.75;
    wz = 1.75;
    
    Ax = [-w+x_r w+x_r -w+y_r w+y_r -wz+z_r wz+z_r];

    % loop through trajectory inputs
    for j = 1:N

        % obtain rotational matrix for current RPY angles
        Rt = R(X(j,4:6)); 
        
        % obtain rotational matrix for pendulum relative to quad
        r = X(j,7);
        s = X(j,8);
        roll_rel = tan(s/pl);       % roll in y-direction (around x-axis)
        pitch_rel = tan(r/pl);      % pitch in x-direction (around y-axis)
        % yaw_rel = X(j,6);
        Rp = R([roll_rel, pitch_rel, 0]);
        
        % define each quadrotor circle
        R1 = Rt*([ l+rc+rx ; ry  ; zeros(size(rx)) ]) + X(j,1:3)';
        R2 = Rt*([ rx ; l+rc+ry  ; zeros(size(rx)) ]) + X(j,1:3)';
        R3 = Rt*([ -l-rc+rx ; ry ; zeros(size(rx)) ]) + X(j,1:3)';
        R4 = Rt*([ rx ; -l-rc+ry ; zeros(size(rx)) ]) + X(j,1:3)';
        
        % define black arms
        A1 = Rt*([-l l;0 0;0 0]) + X(j,1:3)';
        A2 = Rt*([0 0; -l l;0 0]) + X(j,1:3)';
        
        % define magenta inverted pendulum
        P1 = Rp*([0 0;0 0; 0 pl]) + X(j,1:3)';
        
        % verticle reference position of pendulum
        Pv = ([0 0;0 0; 0 pl]) + X(j,1:3)';
        
        % plot coordinate reference
        plot3( x_r,y_r,z_r,'r.');
        hold on
        
        % plot quadrotor frame cross and circles
        plot3( A1(1,:),A1(2,:),A1(3,:),'k',A2(1,:),A2(2,:),A2(3,:),'k');
        plot3( R1(1,:),R1(2,:),R1(3,:),'r',R2(1,:),R2(2,:),R2(3,:),'b',R3(1,:),R3(2,:),R3(3,:),'b',R4(1,:),R4(2,:),R4(3,:),'b');
        
        % plot equilibrium equilibrium (verticle upwards)
        plot3( Pv(1,:),Pv(2,:),Pv(3,:),'k--');
        % plot3( Pv(1,2),Pv(2,2),Pv(3,2),'k');
        
        % plot current pendulum position
        plot3( P1(1,:),P1(2,:),P1(3,:),'m');
        plot3( P1(1,2),P1(2,2),P1(3,2),'m.');
        
        hold off
        grid();
        
        % set axes
        axis(Ax);
        view(3);
        % view([15 25]);
        % view([0 90]);
        
        set(gca,'box','on')
        drawnow   
        if (pause_duration > 0) 
            pause(pause_duration);
        end
    end
    
    % Function to determine the rotation matrix for plotting
    function y = R(Xrot)
        
        phi = Xrot(1);
        theta = Xrot(2);
        psi = Xrot(3);

        Rpsi = [cos(psi)  -sin(psi)   0;
                sin(psi)  cos(psi)    0;
                0         0           1];

        % rotation around y with theta
        Rtheta = [cos(theta)    0       sin(theta);
                  0             1       0;
                 -sin(theta)    0       cos(theta)];

        % rotation around x with phi 
        Rphi = [1       0           0;
                0       cos(phi)    -sin(phi);
                0       sin(phi)    cos(phi)];

        y=Rpsi*Rtheta*Rphi;
    end
    
end