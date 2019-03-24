% main running file for MPC project

clc
clear

%% Simulate System

answer = questdlg('Which simulation to run?', ...
	'Simulations', ...
	'Reference MPC','Trajectory MPC','Simple Quadrotor','Cancel');
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
    case 'Simple Quadrotor'
        close all
        disp(['Running: ' answer])
        run('stability_quad.m');
    case 'Cancel'
        disp('Cancelled -> No simulation run.');
end