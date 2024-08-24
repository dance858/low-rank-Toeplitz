function [grad_x, grad_y, hess_xx, hess_yy, hess_yx] = diff2(x, y, L, n, sig2)
% This function returns the gradient and Hessian of 
% f(x, y) = logdet (T(x, y) + sig2*I) + Tr((T(x, y) + sig2*I)^{-1} S),
% where S is the sample covariance matrix
% Parameters:
%      x - size (n + 1) x 1
%      y - size n x 1
%      L - Cholesky factor of S,  S = L * L', lower triangular
%      n - dimension


N = 2*(n+1);
R = chol(inv(toeplitz([2*x(1) + sig2; x(2:end) + 1i*y])), 'lower'); 

% Compute DFT of columns of R
R_DFT = fft(R, N);

% Compute DFT of columns of A
A = R*(R'*L); 
A_DFT = fft(A, N);

% Compute gradients.
temp = ifft(sum(R_DFT.*conj(R_DFT), 2));
grad1_x = 2*real(temp(1:n+1));
grad1_y = -2*imag(temp(2:n+1));
temp = ifft(sum(A_DFT.*conj(A_DFT), 2));
grad2_x = -2*real(temp(1:n+1));
grad2_y = 2*imag(temp(2:n+1));

% Compute Hessian of f1.
F = R_DFT*R_DFT';
F_F_trans = F.*transpose(F); 

temp1 = ifft(F_F_trans')';                       % (F.*transpose(F))*W/N
temp2 = transpose(ifft(transpose(F_F_trans)));   % (F.*transpose(F))*conj(W)/N;
temp3 = ifft(temp1 + temp2);

hess1_xx = -2*real(temp3(1:n+1, 1:n+1));
hess1_yx = 2*imag(temp3(2:n+1, 1:n+1));
hess1_yy = 2*real(ifft(temp2 - temp1));
hess1_yy = hess1_yy(2:n+1, 2:n+1);

% Compute Hessian of f2.
G = A_DFT*A_DFT';
temp1 = F.*transpose(G) + transpose(F).*G;
temp2 = ifft(temp1')';                         % (F.*transpose(G) + transpose(F).*G)*W/N
temp3 = transpose(ifft(temp1));                % (F.*transpose(G) + transpose(F).*G)*conj(W)/N
% IT HOLDS THAT temp3 = conj(temp2)!!! OMG.

temp4 = ifft(temp2 + temp3);
hess2_xx = 2*real(temp4(1:n+1, 1:n+1));
hess2_yx = -2*imag(temp4(2:n+1, 1:n+1));
hess2_yy = 2*real(ifft(temp2 - temp3));
hess2_yy = hess2_yy(2:n+1, 2:n+1);

grad_x = grad1_x + grad2_x;
grad_y = grad1_y + grad2_y;
hess_xx = hess1_xx + hess2_xx;
hess_yy = hess1_yy + hess2_yy;
hess_yx = hess1_yx + hess2_yx;
end