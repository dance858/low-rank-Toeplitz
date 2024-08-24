function [x0, y0] = initialization(S, X, init_strategy)
% Averages the diagonals of the sample covariance S. Note that x0 
% is scaled by a factor 1/2. 

if init_strategy == 1
    [x0, y0] = block_circulant_initialization(X);
    x0 = transpose(x0);
    y0 = transpose(y0);
elseif init_strategy == 2
    n = size(S, 1) - 1;
    z0 = zeros(n+1, 1);
    for k = 0:n
       z0(k+1) = mean(diag(S, k)); 
    end
    x0 = real(z0); x0(1) = x0(1)/2;
    y0 = imag(z0(2:n+1));

    try
      R = chol(toeplitz([2*x0(1); x0(2:end) + 1i*y0]), 'lower');  
    catch
      [x0, y0] = block_circulant_initialization(X);
      x0 = transpose(x0);
      y0 = transpose(y0);
    end
end
end