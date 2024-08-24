function [cov_est] = LRSE_LB_UB(Y, lb, ub, rank)
% This function returns the covariance matrix and matrix M for benchmarking
% purposes using Shikhaliev et al:s method, with LB and UB

[U, D] = eig(Y * Y');
[d,ind] = sort(diag(D), 'descend');
U = U(:,ind);

num_samples = size(Y, 2);
n = length(d) - 1;
Lambda = zeros(n+1);
dBar = mean(d(rank+1:end));
for k = 1:rank
    Lambda(k, k) = max(d(k)/num_samples, lb);
end

for k = (rank + 1):(n + 1)
    Lambda(k, k) = min(max(dBar/num_samples, lb), ub);
end

cov_est = U*Lambda*U';
end


