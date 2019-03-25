function plot_rotational(t,x,u,other)

    mu_actual = other.mu_actual;
    beta_dot_angle = other.beta_dot_angle;
    gamma_dot_angle = other.gamma_dot_angle;
    beta_angle = other.beta_angle;
    gamma_angle = other.gamma_angle;

    shouldplot = true;
    if shouldplot
        % STATES 1-5
        figure(1);
        clf;
        sgtitle('Pendulum States');
        subplot 411;
        stairs(t, x(1,:), 'm-');  grid();
        ylabel('$p$ [m]','interpreter','latex');

        subplot 412;
        stairs(t, x(2,:), 'm-');  grid();
        ylabel('$\dot{p}$ [m/s]','interpreter','latex');

        subplot 413;
        stairs(t, x(3,:), 'b-');  grid();
        ylabel('$q$ [m]','interpreter','latex');

        subplot 414;
        stairs(t, x(4,:), 'b-');  grid();
        ylabel('$\dot{q}$ [m/s]','interpreter','latex');
        xlabel('Time [s]');

        % STATES 6-10
        figure(2);
        clf;
        sgtitle('Quadrotor States');
        subplot 611;
        stairs(t, x(5,:), 'm-');  grid();
        ylabel('$u$ [m]','interpreter','latex');

        subplot 612;
        stairs(t, x(6,:), 'm-');  grid();
        ylabel('$\dot{u}$ [m/s]','interpreter','latex');

        subplot 613;
        stairs(t, x(7,:), 'b-');  grid();
        ylabel('$v$ [m]','interpreter','latex');

        subplot 614;
        stairs(t, x(8,:), 'b-');  grid();
        ylabel('$\dot{v}$ [m/s]','interpreter','latex');

        subplot 615;
        stairs(t, x(9,:), 'k-');  grid();
        ylabel('$w$ [m]','interpreter','latex');

        subplot 616;
        stairs(t, x(10,:), 'k-');  grid();
        ylabel('$\dot{w}$ [m/s]','interpreter','latex');
        xlabel('Time [s]');

        % EULER ANGLES    
        figure(3);
        clf;
        sgtitle('Euler Angles');
        subplot 211;
        stairs(t, mu_actual, 'm-');  grid();
        ylabel('$\mu$','interpreter','latex');

        subplot 212;
        stairs(t, x(12,:), 'm-');  grid();
        ylabel('$\nu$','interpreter','latex');
        xlabel('Time [s]');

        % INPUTS
        figure(4);
        clf;
        sgtitle('Inputs');
        subplot 311;
        stairs(t, u(1,:), 'm-');  grid();
        ylabel('$\dot{\mu}$','interpreter','latex');

        subplot 312;
        stairs(t, u(2,:), 'm-');  grid();
        ylabel('$\dot{\nu}$','interpreter','latex');

        subplot 313;
        stairs(t, u(3,:), 'm-');  grid();
        ylabel('$a$','interpreter','latex');
        xlabel('Time [s]');

        % INPUTS FOR ROTOR THRUST
        figure(5);
        clf;
        stairs(t,beta_dot_angle);
        hold on
        stairs(t,gamma_dot_angle);
        legend('beta dot','gamma dot');
        title('Derivative of Control Inputs Seen By Quadrotor Props');
        grid();

        figure(6);
        clf;
        stairs(t,beta_angle);
        hold on
        stairs(t,gamma_angle);
        legend('beta','gamma');
        title('Control Inputs Seen By Quadrotor Props');
        grid();
    end
end