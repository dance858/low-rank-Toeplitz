function [obj_value2] = obj2(x, y)
%This function calculates the value of the barrier function phi(x, y)

T = toeplitz([2*x(1); x(2:end) + 1i*y]);

%old code
%obj_value2 = -log(det(T));

L = chol(T, 'lower');
diagonalElements = diag(L);
logDiagonalElements = log(diagonalElements);
logdetT = 2*sum(logDiagonalElements);

obj_value2 = -logdetT;

end