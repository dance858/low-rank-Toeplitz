function [grad_sig, hess_sig_sig, hess_x_sig, hess_y_sig] = diff5(x, y, sig2, S)
% This function returns the gradient and Hessian WITH RESPECT TO sig2 (!!!!) of 
% f(x, y) = logdet (T(x, y) + sig2*I) + Tr((T(x, y) + sig2*I)^{-1} S),
% where S is the sample covariance matrix
% Parameters:
%      x - size (n + 1) x 1
%      y - size n x 1
%      L - Cholesky factor of S,  S = L * L', lower triangular
%      n - dimension

T = toeplitz([2*x(1); x(2:end) + 1i*y]);
[grad_sig1, hess_sig_sig1, hess_x_sig1, hess_y_sig1] = trace_diff_sig2(T, sig2, S);
[grad_sig2, hess_sig_sig2, hess_x_sig2, hess_y_sig2] = logdet_diff_sig2(T, sig2);

grad_sig = grad_sig1 + grad_sig2;
hess_sig_sig = hess_sig_sig1 + hess_sig_sig2;
hess_x_sig = hess_x_sig1 + hess_x_sig2;
hess_y_sig = hess_y_sig1 + hess_y_sig2;
end