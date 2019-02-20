function check_controllability(sysc)

    Ac = sysc.A;
    Bc = sysc.B;

    Ctrb_rank = rank(ctrb(Ac,Bc));
    disp('Number of states');
    disp(size(Ac));
    disp('Rank of controllability matrix');
    disp(Ctrb_rank);

    % Hautus test
    % hautus_test(Ac,Bc);

end