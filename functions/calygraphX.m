%% CALYGRAPH X DETERMINATION
%
% function which determines calygraph X given specific state space matrices
% prediction horizon and terminal set.

function [output, info] = calygraphX(LTI,dim)
    A = LTI.A;
    B = LTI.B;
    
    N_hrzn = dim.N;

end