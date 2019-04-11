function plot_comparison_R_different()

    % INPUTS = time, states_trajectory
    
    load('variables/Q10_R10_Sdare.mat');
    
    time = saved_traj_Q10_R01_Sdare.t;
    states_trajectory = saved_traj_Q10_R01_Sdare.x;

    % Show 6 States of x-direction control
    figure(544);
    clf;
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
    
    load('variables/Q10_R1_Sdare.mat');
    
    time = saved_traj_Q10_R01_Sdare.t;
    states_trajectory = saved_traj_Q10_R01_Sdare.x;
    
    figure(544);
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
    
    load('variables/Q10_R01_Sdare.mat');
    
    time = saved_traj_Q10_R01_Sdare.t;
    states_trajectory = saved_traj_Q10_R01_Sdare.x;
    
    figure(544);
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
    
    load('variables/Q10_R001_Sdare.mat');
    
    time = saved_traj_Q10_R01_Sdare.t;
    states_trajectory = saved_traj_Q10_R01_Sdare.x;
    
    figure(544);
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
    
    load('variables/Q10_R0001_Sdare.mat');
    
    time = saved_traj_Q10_R01_Sdare.t;
    states_trajectory = saved_traj_Q10_R01_Sdare.x;
    
    figure(544);
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
    
    legend('10','1','0.1','0.01','0.001');

%     % Show 6 States of y-direction control
%     show_y_horizontal_performance_plots = true;
%     if show_y_horizontal_performance_plots
%         figure(2);
%         clf;
%         sgtitle('y-direction horizontal motion and roll angles');
%         sgtitle('');
%         subplot 511;
%         stairs(time, states_trajectory(:,7), 'm-');  grid();
%         ylabel('$s$ [m]','interpreter','latex');
% 
%         subplot 512;
%         stairs(time, states_trajectory(:,8), 'm-');  grid();
%         ylabel('$\dot{s}$ [m/s]','interpreter','latex');
% 
%         subplot 513;
%         stairs(time, states_trajectory(:,9), 'b-');  grid();
%         ylabel('$y$ [m]','interpreter','latex');
% 
%         subplot 514;
%         stairs(time, states_trajectory(:,10), 'b-');  grid();
%         ylabel('$\dot{y}$ [m/s]','interpreter','latex');
% 
%         subplot 515;
%         stairs(time, states_trajectory(:,11), 'k-');  grid();
%         ylabel('$\gamma$ [rad]','interpreter','latex');
% 
%         xlabel('Time [s]');
%     end 
end