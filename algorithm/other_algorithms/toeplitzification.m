function [T] = toeplitzification(A)
dim = size(A, 1);
n = dim - 1;
t = zeros(n+1, 1);
for k = 0:n
       t(k+1) = mean(diag(A, k)); 
end
T = toeplitz(t);
end