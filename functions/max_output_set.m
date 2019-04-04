function [H,h] = max_output_set(A,K,u_lim,x_lim)

% Pedro's newest, freshest code on the block

% Input constraints
f1 = u_lim*ones(2*size(K,1),1);

% State constraints
f2 = [x_lim; x_lim];

f = [f1;f2];
s = size(f,1);

%% Algorithm implementation
exit_flag=0;
k=0;

A_in=[];
b = [];

opts = optimoptions('linprog','Display','off');

while exit_flag==0 
    fprintf('\tk = %i \n',k);
    A_in = [A_in; K*A^k; -K*A^k; A^k; -A^k];
    b = [b; f];
     
    h = [K*A^(k+1); -1*K*A^(k+1); A^(k+1); -1*A^(k+1)];
     
    opt_val = zeros(1,s);
    for i=1:s                
        [~,fval,exit] = linprog(-h(i,:),A_in,b,[],[],[],[],opts);
        if exit == -3
            fval = inf;
        end
        opt_val(i)=-fval-f(i);
    end
    
    if (all(opt_val<=0-eps) && exit~=-3)
        exit_flag=1;
        H=A_in;
        h=b;
        fprintf('\tdone!\n');
    else
        k = k+1;
    end
end

end