function sysd_cl = calc_lqr_system(A,B,C,K,h)
    %% DEFINE THE CLOSED LOOP SYSTEM WITH REFERENCE TRACKING INPUT

    % reference tracking input control
    B_ref = [0 0 0 0;  % r1
             0 0 0 0;  % r2
             1 0 0 0;  % x1 
             0 0 0 0;  % x2
             0 0 0 0;  % beta1
             0 0 0 0;  % beta2
             0 0 0 0;  % s1 
             0 0 0 0;  % s2
             0 1 0 0;  % y1
             0 0 0 0;  % y2
             0 0 0 0;  % gamma1
             0 0 0 0;  % gamma2
             0 0 1 0;  % z1
             0 0 0 0;  % z2
             0 0 0 1;  % 
             0 0 0 0]; % 

    % define closed loop system with LQR control law
    sysd_cl_unnormalized = ss(A-B*K,B_ref,C,[],h);

    % normalize closed-loop reference tracking gains 
    dcgain_cl = dcgain(sysd_cl_unnormalized);
    B_ref(3,1) = 1/dcgain_cl(3,1);
    B_ref(9,2) = 1/dcgain_cl(9,2);
    B_ref(13,3) = 1/dcgain_cl(13,3);
    B_ref(15,4) = 1/dcgain_cl(15,4);

    % define closed-loop system with normalized input gains
    sysd_cl = ss(A-B*K,B_ref,C,[],h);
    % dcgain(sysd_cl);

end