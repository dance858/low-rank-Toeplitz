function [out] = LRSE(S, sig2, rank)

% This function returns the covariance matrix and matrix M for benchmarking
% purposes using Shikhaliev et al:s method


[U, D] = eig(S);
d = flip(diag(D));
n = length(d) - 1;
Lambda = zeros(n + 1);

for k = 1:rank
    Lambda(k, k) = max(d(k), sig2);
end

for k = (rank + 1):(n + 1)
    Lambda(k, k) = sig2;
end

cov_est = U*Lambda*U';

out = cov_est;

end


