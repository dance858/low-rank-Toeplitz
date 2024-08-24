function [grad, hess] = diff4(x)
% This function returns the gradient and Hessian of 
% f(x, y) = 2*(n+1)*x0

n = length(x) - 1;

grad = zeros(n + 1, 1);
grad(1) = 2*(n + 1);

hess = zeros(n + 1, n + 1);

end


