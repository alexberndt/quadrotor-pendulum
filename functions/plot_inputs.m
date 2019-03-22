function plot_inputs(t,u,u_limit)
    
    % plot inputs over time
    figure(3);
    clf;
    stairs(t,u(1,:));
    hold on
    stairs(t,u(2,:));
    stairs(t,u(3,:));
    stairs(t,u(4,:));
    
    plot([t(1),t(end)],[u_limit u_limit],'k--');
    plot([t(1),t(end)],[-u_limit -u_limit],'k--');
    
    ylim([-1.1*u_limit 1.1*u_limit]);
    
    grid();
    % 
end