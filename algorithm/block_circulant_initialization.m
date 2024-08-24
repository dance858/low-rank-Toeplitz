function [x0, y0] = block_circulant_initialization(X)

[n, N] = size(X);
n = n - 1;
%% Block circulant ML estimate
F = zeros(n+1, n+1);
omega = exp(2*pi*1j/(n+1));
row = omega.^(0:n); 
for k = 0:n
   F(k+1, :) = row.^(-k); 
end
U = 1/sqrt(n+1)*F;
y_DFT = U'*X;

s = 1/N*sum(y_DFT.*conj(y_DFT), 2);
R = U*diag(s)*U';
x0 = real(R(1, :));
x0(1) = x0(1)/2;
y0 = imag(R(1, 2:end));

%% same method more efficient computation.
%Z = X;
%Z_DFT = 1/sqrt(n+1)*(n+1)*ifft(Z);     % equal to y_DFT
%q = 1/N*sum(Z_DFT.*conj(Z_DFT), 2);    % equal to s
%Q = diag(q);
%R_circ = 1/(n+1)*fft(fft(Q')');

%first_col_R_circ = 1/(n+1)*fft(q);



end