function [grad_x, grad_y, hess_xx, hess_yy, hess_yx] = diff3(x, y)
% This function returns the gradient and Hessian of 
% f(x, y) = -logdet T(x, y) 
% Parameters:
%      x - size (n + 1) x 1
%      y - size n x 1

n = length(y);
N = 2*(n+1);
R = chol(inv(toeplitz([2*x(1); x(2:end) + 1i*y])), 'lower'); 

% Compute DFT of columns of R
R_DFT = fft(R, N);

% Compute gradients.
temp = ifft(sum(R_DFT.*conj(R_DFT), 2));
grad_x = 2*real(temp(1:n+1));
grad_y = -2*imag(temp(2:n+1));

% Compute Hessian of f1.
F = R_DFT*R_DFT';
F_F_trans = F.*transpose(F); 

temp1 = ifft(F_F_trans')';                       % (F.*transpose(F))*W/N
temp2 = transpose(ifft(transpose(F_F_trans)));   % (F.*transpose(F))*conj(W)/N;
temp3 = ifft(temp1 + temp2);

hess_xx = -2*real(temp3(1:n+1, 1:n+1));
hess_yx = 2*imag(temp3(2:n+1, 1:n+1));
hess_yy = 2*real(ifft(temp2 - temp1));
hess_yy = hess_yy(2:n+1, 2:n+1);

grad_x = -grad_x;
grad_y = -grad_y;

hess_xx = -hess_xx;
hess_yx = -hess_yx;
hess_yy = -hess_yy;


end