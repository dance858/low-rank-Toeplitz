function [value] = barrier_function_ub_lb(x, y, sig2, S, t, reg, ub, lb)
%This function calculates the value of the entire barrier function

n = length(y);

value = t*(obj1(x, y, sig2, S) + 2*reg*(n + 1)*x(1)) + obj2(x, y) - log(sig2 - lb) - log(ub - sig2);
end