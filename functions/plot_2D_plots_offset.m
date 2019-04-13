function plot_2D_plots_offset(time, states_trajectory, offset_sequence, reference_sequence)


    
    % Show 6 States of x-direction control
    show_x_horizontal_performance_plots = true;
    if show_x_horizontal_performance_plots
        figure(547);
        clf;
        %sgtitle('x-direction horizontal motion and pitch angles');
        sgtitle('');
        subplot 311;
        stairs(time, states_trajectory(:,1), 'm-');  grid(); 
        ylabel('$r$ [m]','interpreter','latex');
        ylim([-0.055 0.030]);

%         subplot 512;
%         stairs(time, states_trajectory(:,2), 'm-');  grid();
%         ylabel('$\dot{r}$ [m/s]','interpreter','latex');

        subplot 312;
        stairs(time, states_trajectory(:,3), 'b-');  grid(); hold on;
        stairs(time, 25*offset_sequence(1,:), 'k-');
        stairs(time, reference_sequence, 'r-');
        ylabel('$x$ [m]','interpreter','latex');
        ylim([-1.05 0.34]);
        legend('x','disturbance','reference');
   

%         subplot 514;
%         stairs(time, states_trajectory(:,4), 'b-');  grid();
%         ylabel('$\dot{x}$ [m/s]','interpreter','latex');

        subplot 313;
        stairs(time, states_trajectory(:,5));  grid();
        ylabel('$\beta$ [rad]','interpreter','latex');
        ylim([-0.14 0.07]);

        xlabel('Time [s]');
    end

    % Show 6 States of y-direction control
    show_y_horizontal_performance_plots = true;
    if show_y_horizontal_performance_plots
        figure(2);
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