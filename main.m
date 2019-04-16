% main running file for MPC project

clc
clear

%% Simulate System

% close all
% run('adaptive_LQR.m');

x = input('Which simulation to run?\n\n 1) Regulation MPC \n 2) Regulation LQR \n 3) Disturbance Output MPC \n 4) Adaptive MPC \n 5) Adaptive LQR \n 6) Stability Analysis \n \n'); 

switch x
    case 1
        close all;
        run('reference_MPC.m');
    case 2
        close all;
        run('reference_LQR.m');
    case 3
        close all;
        run('output_MPC_w_offset.m');
    case 4
        close all;
        run('adaptive_MPC.m');
    case 5
        close all;
        run('adaptive_LQR.m');
    case 6
        close all;
        run('stability_analysis.m');
    otherwise
        disp('\n\nPlease select a valid input (just a number) \n\n Run the main.m script again');
end
        

% %%
% 
% answer = questdlg('Which simulation to run?', ...
% 	'Simulations', ...
% 	'Regulation MPC','Disturbance Output MPC','Stability Analysis','Cancel'); % 'Stability Analysis',
% % Handle response
% switch answer
% %     case 'Reference MPC'
% %         close all
% %         disp(['Running: ' answer]) 'Reference MPC',
% %         run('reference_MPC.m');
%     case 'Regulation MPC'
%         close all
%         disp(['Running: ' answer])
%         run('reference_MPC.m');
%     case 'Disturbance Output MPC'
%         close all
%         disp(['Running: ' answer])
%         run('output_MPC_w_offset.m');
%     case 'Stability Analysis'
%         close all
%         disp(['Running: ' answer])
%         run('stability_analysis.m');
%     case 'Cancel'
%         disp('Cancelled -> No simulation run.');
% end