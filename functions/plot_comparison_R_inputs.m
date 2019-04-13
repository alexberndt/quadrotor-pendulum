function plot_comparison_R_inputs()
    % t,u,u_limit
    fignum = 589;
    u_limit = 0.1;
    
    % plot inputs over time
    figure(fignum);
    clf;
    hold on;
    
%     load('variables/R_inputs/QS1_R0001.mat');
%     u = saved_data.u;
%     t = saved_data.t;
%     u_norm = vecnorm(u,2);
%     u_norm = u(1,:);
%     stairs(t,u_norm);
    
    load('variables/R_inputs/QS1_R001.mat');
    u = saved_data.u;
    t = saved_data.t;
    u_norm = vecnorm(u,2);
    u_norm = u(1,:);
    stairs(t,u_norm);
    
%     load('variables/R_inputs/QS1_R01.mat');
%     u = saved_data.u;
%     t = saved_data.t;
%     u_norm = vecnorm(u,2);
%     u_norm = u(1,:);
%     stairs(t,u_norm);
    
    load('variables/R_inputs/QS1_R1.mat');
    u = saved_data.u;
    t = saved_data.t;
    u_norm = vecnorm(u,2);
    u_norm = u(1,:);
    stairs(t,u_norm);
    
%     load('variables/R_inputs/QS1_R10.mat');
%     u = saved_data.u;
%     t = saved_data.t;
%     u_norm = vecnorm(u,2);
%     u_norm = u(1,:);
%     stairs(t,u_norm);
    
    load('variables/R_inputs/QS1_R100.mat');
    u = saved_data.u;
    t = saved_data.t;
    u_norm = vecnorm(u,2);
    u_norm = u(1,:);
    stairs(t,u_norm,'m-');
    
    ylim([-1.1*u_limit 1.1*u_limit]);
%     ylim([-0.05 0.22]);
    xlim([0 3.9]);
    
    plot([t(1),t(end)],[u_limit u_limit],'k--');
    plot([t(1),t(end)],[-u_limit -u_limit],'k--');
    
    legend('0.01','1','100');
    
    
    grid();
    % 
end