function [out] = Newton_ub_lb(x0, y0, S, sig2_0, reg, lb, ub, parameters)
%  This function implements a barrier method for solving
%  min. f(x, y, sig2) subject to T(x, y) >= 0,  ub >= sig2 >= lb 
%  where
%  f(x, y, sig2) = log det(T(x,y) + sig2 * I) + Tr((T(x,y) + sig2*I)^{-1} S) +
%             2*reg * (n+1)x0.
%  The variables are x = [x0 x1 ... xn] (size n + 1) and y = [y1 ... yn]
%  (size n).

% --------------------------------------------------
%       Parse algorithmic parameters 
% --------------------------------------------------
alpha = parameters.alpha; beta = parameters.beta; max_iter = parameters.max_iter;
verbose = parameters.verbose; n = length(x0) - 1; tol = parameters.tol;
t = parameters.t; mu = parameters.mu;

% ------------------------------------------------------------------------
%   Compute a Cholesky factor of S, taking numerics into account 
%   (we truncate negative eigenvalues of the sample covariance)
% ------------------------------------------------------------------------
x = x0;
y = y0;
sig2 = sig2_0;
try
  L = chol(S, 'lower');    % S = L*L'  
catch
  [evecs, evals] = eig(S);
  evals_pos = max(diag(evals), 1e-18); 
  L = evecs(:, (evals_pos >= 1e-11))*diag(sqrt(evals_pos(evals_pos >= 1e-11)));
end

for iter_outer = 1:max_iter
        for iter_inner = 1:max_iter
            assert(sig2>0);
            assert(norm(x) > 1e-6);
            % Calculate gradient and Hessian of f(x, y, sig2)
            [grad_f_x, grad_f_y, hess_f_xx, hess_f_yy, hess_f_yx] = diff2(x, y, L, n, sig2);
            [grad_f_sig, hess_f_sig_sig, hess_f_x_sig, hess_f_y_sig] = diff5(x, y, sig2, S);
            grad_f = [grad_f_x; grad_f_y; grad_f_sig];
            hess_f_xy = hess_f_yx';
            hess_f_sig_x = hess_f_x_sig';
            hess_f_sig_y = hess_f_y_sig';
            hess_f = [hess_f_xx, hess_f_xy, hess_f_x_sig; hess_f_yx, hess_f_yy, hess_f_y_sig; hess_f_sig_x, hess_f_sig_y, hess_f_sig_sig];

             % Calculate gradient and Hessian of phi(x, y)
             [grad_phi_x, grad_phi_y, hess_phi_xx, hess_phi_yy, hess_phi_yx] = diff3(x, y);
             grad_phi = [grad_phi_x; grad_phi_y];
             hess_phi_xy = hess_phi_yx';
             hess_phi = [hess_phi_xx, hess_phi_xy; hess_phi_yx, hess_phi_yy];

             % Calculate gradient and Hessian of regularizer
             [grad_reg_x, ~] = diff4(x);
            
             % Calculate gradient and Hessian of the entire barrier function
             grad = t*grad_f;
             hess = t*hess_f;

             grad(1:2*n + 1) = grad(1:2*n + 1) + grad_phi;
             grad(1:n + 1) = grad(1:n + 1) + t*reg*grad_reg_x;
             hess(1:2*n + 1, 1:2*n + 1) = hess(1:2*n + 1, 1:2*n + 1) + hess_phi; 
             grad(end) = grad(end) - 1/(sig2 - lb) + 1 / (ub - sig2);
             hess(end, end) = hess(end, end) + 1/((sig2 - lb)^2) + 1/((ub - sig2)^2);

             try
                 L1 = chol(hess, 'lower');
             catch
                 hess = hess - 1.1*min(min(eig(hess)), -1e-5)*eye(2*(n+1));
                 L1 = chol(hess, 'lower');
             end
           
             newton_step_u = - L1' \ (L1 \ grad);
             newton_dec = -grad'*newton_step_u;
             stop_criterion = (newton_dec) / 2;
            
            if verbose == 1
                fprintf('Inner decrement: %.10f\n', stop_criterion);
            end

            % Exit Newton's method if stop criterion is met
            if stop_criterion <= tol
                break;
            end

            %Calculate the newton step for x and y, separately
            newton_step_x = newton_step_u(1:n+1);
            newton_step_y = newton_step_u(n+2:end-1);
            newton_step_sig2 = newton_step_u(end);
            
            % Perform backtracking line search
            step_length = 1.0;
            while true
                x_new = x + step_length * newton_step_x;
                y_new = y + step_length * newton_step_y;
                sig2_new = sig2 + step_length * newton_step_sig2;

                %Making sure the toeplitz matrix is positive definite
                toep = toeplitz([2*x_new(1); x_new(2:end) + 1i*y_new]);

                try
                    not_used = chol(toep, 'lower');
                    success = true;
                catch
                    success = false;
                end

                if (sig2_new <= lb || sig2_new >= ub)
                    success = false;
                end

                if success
                    old_obj = barrier_function_ub_lb(x, y, sig2, S, t, reg, ub, lb);
                    new_obj = barrier_function_ub_lb(x_new, y_new, sig2_new, S, t, reg, ub, lb);
                    while new_obj - 1e-5 * abs(new_obj) > old_obj + alpha * step_length * grad' * newton_step_u
                        step_length = step_length * beta;
                        x_new = x + step_length * newton_step_x;
                        y_new = y + step_length * newton_step_y;
                        sig2_new = sig2 + step_length * newton_step_sig2;
                        new_obj = barrier_function_ub_lb(x_new, y_new, sig2_new, S, t, reg, ub, lb);
                    end
                    x = x_new;
                    y = y_new;
                    sig2 = sig2_new;
                    break;
                else
                    step_length = step_length / 2;
                end
            end
        end
        
        if verbose == 1
            fprintf('Outer decrement: %f\n', (n + 1) / t);
        end

        % Exit outer loop if outer decrement is less than tolerance
        if (n + 1) / t < tol %&& res1 < tol && res2 < tol
            break;
        else
            t = t * mu; % Increase t and restart Newton's method
        end
end

    
    % Return the final solution
    out.x = x;
    out.y = y;
    out.sig2 = sig2;
    out.T = toeplitz([2*x(1); x(2:end) + 1i*y]);
end
