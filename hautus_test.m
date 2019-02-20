function hautus_test(A,B)
    eigvals = eig(A);
    for idx = 1:numel(eigvals)
     lambda = eigvals(idx);
     fprintf('Eigenvalue: %f \t -> ', lambda);
     rk = rank([(eye(size(A))*lambda-A) B]);
     fprintf('rank: %i \n',rk);
    end
end