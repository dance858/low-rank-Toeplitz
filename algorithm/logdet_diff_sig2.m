function [grad, hess_sig_sig, hess_x_sig, hess_y_sig] = logdet_diff_sig2(T, sig2)
dim = size(T, 1);
n = dim - 1;
T_sig2_inv = inv(T + sig2 * eye(dim));
grad = trace(T_sig2_inv);
hess_sig_sig = - trace(T_sig2_inv*T_sig2_inv);

N = 2*(n+1);
% Compute DFT of columns of A
A = T_sig2_inv;
A_DFT = fft(A, N);
temp = ifft(sum(A_DFT.*conj(A_DFT), 2));
hess_x_sig = -2*real(temp(1:n+1));
hess_y_sig = 2*imag(temp(2:n+1));

end