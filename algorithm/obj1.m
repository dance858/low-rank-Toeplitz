function [obj_value1] = obj1(x, y, sig2, S)
%This function calculates the value of the objective function f(x, y)

T = toeplitz([2*x(1); x(2:end) + 1i*y]);
n = length(x) - 1;

%Previous code which do not ensure numeric stability
%obj_value1 = log(det(T + sig2*eye(n + 1))) + trace((T + sig2*eye(n + 1))\S);

%New code ensuring numeric stability for log(det(A))
A = T + sig2*eye(n + 1);

L = chol(A, 'lower');
diagonalElements = diag(L);
logDiagonalElements = log(diagonalElements);
logdetA = 2*sum(logDiagonalElements);

obj_value1 = logdetA + real(trace((T + sig2*eye(n + 1))\S));


end