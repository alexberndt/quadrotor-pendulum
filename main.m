% main running file for MPC project

clc
clear

%% Simulate System

close all
run('adaptive_LQR.m');

answer = questdlg('Which simulation to run?', ...
	'Simulations', ...
	'Regulation MPC','Disturbance Output MPC','Stability Analysis','Cancel'); % 'Stability Analysis',
% Handle response
switch answer
%     case 'Reference MPC'
%         close all
%         disp(['Running: ' answer]) 'Reference MPC',
%         run('reference_MPC.m');
    case 'Regulation MPC'
        close all
        disp(['Running: ' answer])
        run('reference_MPC.m');
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