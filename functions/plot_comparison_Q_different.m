function plot_comparison_Q_different()

    % INPUTS = time, states_trajectory
    
    load('variables/QS_same/QS_0001_R1.mat');
    
    time = saved_traj_Q10_R01_Sdare.t;
    states_trajectory = saved_traj_Q10_R01_Sdare.x;
    
    fignum = 546;

    % Show 6 States of x-direction control
    figure(fignum);
    clf;
    % clf;
    %sgtitle('x-direction horizontal motion and pitch angles');
    sgtitle('');
    subplot 311;
    stairs(time, states_trajectory(:,1));  grid on; hold on;
    ylabel('$r$ [m]','interpreter','latex');

%     subplot 512;
%     stairs(time, states_trajectory(:,2));  grid on; hold on;
%     ylabel('$\dot{r}$ [m/s]','interpreter','latex');

    subplot 312;
    stairs(time, states_trajectory(:,3));  grid on; hold on;
    ylabel('$x$ [m]','interpreter','latex');

%     subplot 514;
%     stairs(time, states_trajectory(:,4));  grid on; hold on;
%     ylabel('$\dot{x}$ [m/s]','interpreter','latex');

    subplot 313;
    stairs(time, states_trajectory(:,5));  grid on; hold on;
    ylabel('$\beta$ [rad]','interpreter','latex');

    xlabel('Time [s]');

    hold on;
    
    load('variables/QS_same/QS_0005_R1.mat');
    
    time = saved_traj_Q10_R01_Sdare.t;
    states_trajectory = saved_traj_Q10_R01_Sdare.x;
    
    figure(fignum);
    hold on;
    subplot 311;
    stairs(time, states_trajectory(:,1));
    ylim([-0.1 0.1]);
%     subplot 512;
%     stairs(time, states_trajectory(:,2));
    subplot 312;
    stairs(time, states_trajectory(:,3));
    ylim([-0.05 0.8]);
%     subplot 514;
%     stairs(time, states_trajectory(:,4));
%     ylim([-0.7 1.0]);
    subplot 313;
    stairs(time, states_trajectory(:,5));
    ylim([-0.35 0.35]);
    
    load('variables/QS_same/QS_001_R1.mat');
    
    time = saved_traj_Q10_R01_Sdare.t;
    states_trajectory = saved_traj_Q10_R01_Sdare.x;
    
    figure(fignum);
    hold on;
    subplot 311;
    stairs(time, states_trajectory(:,1));
%     subplot 512;
%     stairs(time, states_trajectory(:,2));
    subplot 312;
    stairs(time, states_trajectory(:,3));
%     subplot 514;
%     stairs(time, states_trajectory(:,4));
    subplot 313;
    stairs(time, states_trajectory(:,5));
    
    load('variables/QS_same/QS_01_R1.mat');
    
    time = saved_traj_Q10_R01_Sdare.t;
    states_trajectory = saved_traj_Q10_R01_Sdare.x;
    
    figure(fignum);
    hold on;
    subplot 311;
    stairs(time, states_trajectory(:,1));
%     subplot 512;
%     stairs(time, states_trajectory(:,2));
    subplot 312;
    stairs(time, states_trajectory(:,3));
%     subplot 514;
%     stairs(time, states_trajectory(:,4));
    subplot 313;
    stairs(time, states_trajectory(:,5));
    
    
    
    load('variables/QS_same/QS_10_R1.mat');
    
    time = saved_traj_Q10_R01_Sdare.t;
    states_trajectory = saved_traj_Q10_R01_Sdare.x;
    
    figure(fignum);
    hold on;
    subplot 311;
    stairs(time, states_trajectory(:,1));
    ylim([-0.12 0.12]);
%     subplot 512;
%     stairs(time, states_trajectory(:,2));
%     ylim([-0.35 0.2]);
    subplot 312;
    stairs(time, states_trajectory(:,3));
    ylim([-0.05 1.1]);
%     subplot 514;
%     stairs(time, states_trajectory(:,4));
%     ylim([-0.9 1.2]);
    subplot 313;
    stairs(time, states_trajectory(:,5));
    ylim([-0.3 0.3]);
    
    legend('0.001','0.005','0.01','0.1','10');
end