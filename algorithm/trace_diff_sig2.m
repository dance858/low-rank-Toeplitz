function [grad_sig, hess_sig_sig, hess_x_sig, hess_y_sig] = trace_diff_sig2(T, sig2, S)
dim = size(T, 1);
T_sig_2 = T + sig2 * eye(dim);
T_sig2_inv = inv(T_sig_2);

grad_sig = - real(trace(T_sig2_inv^2*S));
hess_sig_sig  = 2*real(trace(T_sig2_inv^3*S));


temp = T_sig2_inv^2 * S * T_sig2_inv;
temp = temp + temp';
[hess_x_sig, hess_y_sig] = adjoint(temp); 
end
