function [out] = ITAM(S, target_rank, iter)
% This function returns the covariance matrix benchmarking
%purposes

T = S;

for k = 1:iter
    [U, Sigma] = eig(T);
    [d,ind] = sort(diag(Sigma), 'descend');
    Sigma = diag(d);
    U = U(:,ind);
    
    U1 = U(:, 1:target_rank);
    Sigma1 = Sigma(1:target_rank, 1:target_rank);
    lambda_average = mean(diag(Sigma(target_rank + 1:end, target_rank + 1:end)));
    
    R_s = U1*(Sigma1 - lambda_average*eye(target_rank))*U1';
    T = toeplitzification(R_s);
end

out = T;

end


