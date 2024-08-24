function [x, y] = adjoint(Z)
subdiag_sums = [];

n = size(Z, 1);
for k = 0:n-1
    subdiag = diag(Z, k);  
    subdiag_sum = sum(subdiag);  
    subdiag_sums = [subdiag_sums, subdiag_sum]; 
end

res_x = real(subdiag_sums)';
res_y = imag(subdiag_sums(2:end))';
x = 2 * res_x;
y = 2 * res_y;
end