% main running file for MPC project

clc
clear

%% Simulate System

answer = questdlg('Which simulation to run?', ...
	'Simulations', ...
	'Reference MPC','Trajectory MPC','Stability Analysis','Cancel');
% Handle response
switch answer
    case 'Reference MPC'
        close all
        disp(['Running: ' answer])
        run('reference_MPC.m');
    case 'Trajectory MPC'
        close all
        disp(['Running: ' answer])
        run('adaptive_LQR.m');
    case 'Output MPC'
        close all
        disp(['Running: ' answer])
        run('output_MPC.m');
    case 'Disturbance Output MPC'
        close all
        disp(['Running: ' answer])
        run('output_MPC_w_offset.m');
    case 'Stability Analysis'
        close all
        disp(['Running: ' answer])
        run('stability_analysis.m');
    case 'Cancel'
        disp('Cancelled -> No simulation run.');
end