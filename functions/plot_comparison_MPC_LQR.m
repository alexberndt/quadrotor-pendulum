function plot_comparison_MPC_LQR()

    % INPUTS = time, states_trajectory
    
    load('variables/MPC_vs_LQR/LQR_Q_01_R_700.mat');
    fignum = 567;
    time = saved_data.t;
    states_trajectory = saved_data.x;

    % Show 6 States of x-direction control
    figure(fignum);
    clf;
    %sgtitle('x-direction horizontal motion and pitch angles');
    sgtitle('');
    subplot 311;
    stairs(time, states_trajectory(:,1));  grid on; hold on;
    ylabel('$r$ [m]','interpreter','latex');
    subplot 312;
    stairs(time, states_trajectory(:,3));  grid on; hold on;
    ylabel('$x$ [m]','interpreter','latex');
    subplot 313;
    stairs(time, states_trajectory(:,5));  grid on; hold on;
    ylabel('$\beta$ [rad]','interpreter','latex');
    xlabel('Time [s]');

    hold on;
    
    load('variables/MPC_vs_LQR/MPC_Q_1_R_1.mat');
    time = saved_data.t;
    states_trajectory = saved_data.x;
    figure(fignum);
    hold on;
    subplot 311;
    stairs(time, states_trajectory(:,1));
    subplot 312;
    stairs(time, states_trajectory(:,3));
    subplot 313;
    stairs(time, states_trajectory(:,5));
    
    legend('LQR','MPC');

end