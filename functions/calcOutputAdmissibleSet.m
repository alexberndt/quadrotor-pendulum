function [output] = calcOutputAdmissibleSet()
    % Function which determines the invariant set X_f for a given A,B,C
    % state space system with an LQR gain of K_LQR
    LTI.A = A;
    LTI.B = B;
    LTI.C = C;
    LTI.K = K_LQR;
    
    dim.nx = 16;
    dim.nu = 4;
    
    % Definition of the maximal output admissible set
    % Output constraints
    
%     F = [1; -1];
%     f = ones(2,1);
    u_limit = 1;
    
    F = kron(ones(dim.nx,1),[1; -1]);
    f = u_limit*ones(2*dim.nx,1);

    s = size(F,1);

    %%
    % Algorithm implementation
    exit_flag=0;
    k=0;

    while (exit_flag == 0) 
        
         % Set the constraints for the optimization problem
         A=[];
         
         for t = 1:k+1
             for j = 1:s
                A_aux(j,:) = F(j,:)*K_LQR*LTI.A^(t-1);
             end
             
             A=[A; A_aux];
             
             clear A_aux
         end
         b = repmat(f,k+1,1);

        for i=1:s
             % Set the optimization problem: objective function and constraints
             h=F(i)*LTI.C*LTI.A^(k+1);

             % Solve the optimization problem with CVX
             cvx_precision best
             cvx_begin 
                variable x_opt(dim.nx)
                maximize(h*x_opt-f(i))
                subject to
                A*x_opt-b<=0;
            cvx_end

            % Save optimal value
            opt_val(i)=cvx_optval;
        end

        % Evaluating optimality condition
        if (sum(opt_val<=0-eps)==s && strcmp(cvx_status,'Unbounded')==0)
            exit_flag=1;
            k_star=k
            H=A
            h=b
        else
            clear opt_val A_aux A b h
            k=k+1;
        end
    end


end