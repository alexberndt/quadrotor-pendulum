function [H,h] = max_out_set_pedro(A,K,u_lim)

    % Definition of system dimension
    nx=size(K,2);     % state dimension
    ny=size(K,1);     % output dimension

    % Definition of the maximal output admissible set
    % Output constraints
    f=u_lim*ones(2*ny,1);

    % F=[1; -1];
    % f=ones(2,1);

    s=size(f,1);

    %%
    % Algorithm implementation
    exit_flag=0;
    k=0;

    A_in=[];
    b = [];
    while exit_flag==0 

        A_in = [A_in; [K*A^k; -K*A^k]];
        b = [b; f];

        h = [K*A^(k+1); -1*K*A^(k+1)];

        for i=1:s         
    %          cvx_precision best
    %          cvx_begin 
    %             variable x_opt(nx)
    %             maximize(h(i,:)*x_opt-f(i))
    %             subject to
    %             A_in*x_opt-b<=0;
    %         cvx_end
    %         
    %         % Save optimal value
    %         opt_val(i)=cvx_optval;

            [x,fval,exit] = linprog(-h(i,:),A_in,b);
            if exit == -3
                fval = inf;
            end
            opt_val(i)=-fval-f(i);
        end
        opt_val
        % Evaluating optimality condition
    %     if (sum(opt_val<=0-eps)==s && strcmp(cvx_status,'Unbounded')==0)
    %         exit_flag=1;
    %         H=A_in;
    %         h=b;
    %     else
    %         clear opt_val A_aux A_in b h
    %         k=k+1;
    %     end

        if (sum(opt_val<=0-eps)==s && exit~=-3)
            exit_flag=1;
            k_star=k;
            H=A_in;
            h=b;
        else
            k=k+1;
        end
    end

end