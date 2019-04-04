function [Xf_set_H, Xf_set_h, k_star] = calcInvariantXf(A_K,C_aug,F,f,s,dim)
    % Algorithm implementation
    exit_flag = 0;
    k = 0;
    A = [];

    while (exit_flag == 0)  
        fprintf('%i',k);
        % Set the constraints for the optimization problem
        % A=[];
        % for t=1:k+1
        %     for j=1:s
        %         A_aux(j,:) = F(j,:)*C_aug*A_K^(k);
        %     end
        %     A=[A; A_aux];
        %     clear A_aux
        % end
        A_aux = F(:,:)*C_aug*A_K^(k);
        A = [A; A_aux];
        b = repmat(f,k+1,1);

        for i = 1:s
             % Set the optimization problem: objective function and constraints
             h = F(i,:)*C_aug*A_K^(k+1);

             % Solve the optimization problem with CVX
             cvx_precision best
             cvx_begin quiet
                variable x_opt(dim.nx)
                maximize(h*x_opt-f(i))
                subject to
                 A*x_opt - b <= 0;
             cvx_end
             % disp(k);

            % Save optimal value
            opt_val(i)=cvx_optval;
        end

        % Evaluating optimality condition
        if (sum(opt_val <= 0 - eps) == s && strcmp(cvx_status,'Unbounded') == 0)
            exit_flag = 1;
            k_star = k;
            Xf_set_H = A;
            Xf_set_h = b;
        else
            clear opt_val b h
            k=k+1;
        end
    end

end