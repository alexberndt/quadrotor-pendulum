function plot_2D_plots(time, states_trajectory)

    % Show 6 States of x-direction control
    show_x_horizontal_performance_plots = true;
    if show_x_horizontal_performance_plots
        figure(1);
        % clf;
        %sgtitle('x-direction horizontal motion and pitch angles');
        sgtitle('');
        subplot 511;
        stairs(time, states_trajectory(:,1), 'm-');  grid();
        ylabel('$r$ [m]','interpreter','latex');

        subplot 512;
        stairs(time, states_trajectory(:,2), 'm-');  grid();
        ylabel('$\dot{r}$ [m/s]','interpreter','latex');

        subplot 513;
        stairs(time, states_trajectory(:,3), 'b-');  grid();
        ylabel('$x$ [m]','interpreter','latex');

        subplot 514;
        stairs(time, states_trajectory(:,4), 'b-');  grid();
        ylabel('$\dot{x}$ [m/s]','interpreter','latex');

        subplot 515;
        stairs(time, states_trajectory(:,5), 'k-');  grid();
        ylabel('$\beta$ [rad]','interpreter','latex');

        xlabel('Time [s]');
    end

    % Show 6 States of y-direction control
    show_y_horizontal_performance_plots = true;
    if show_y_horizontal_performance_plots
        figure(2);
        clf;
        sgtitle('y-direction horizontal motion and roll angles');
        sgtitle('');
        subplot 511;
        stairs(time, states_trajectory(:,7), 'm-');  grid();
        ylabel('$s$ [m]','interpreter','latex');

        subplot 512;
        stairs(time, states_trajectory(:,8), 'm-');  grid();
        ylabel('$\dot{s}$ [m/s]','interpreter','latex');

        subplot 513;
        stairs(time, states_trajectory(:,9), 'b-');  grid();
        ylabel('$y$ [m]','interpreter','latex');

        subplot 514;
        stairs(time, states_trajectory(:,10), 'b-');  grid();
        ylabel('$\dot{y}$ [m/s]','interpreter','latex');

        subplot 515;
        stairs(time, states_trajectory(:,11), 'k-');  grid();
        ylabel('$\gamma$ [rad]','interpreter','latex');

        xlabel('Time [s]');
    end 
end