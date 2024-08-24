function [out] = lowrankToep(x0, y0, S, sig2_0, reg0, parameters, target_rank, lb_sig2, ub_sig2, adaptive_reg)
% This function runs newton's method for different values of reg, 
% tuning it so the returned Toeplitz matrix is of the desired rank
reg = reg0;

eval_factor = 1e-2;
Newton_out = Newton_ub_lb(x0, y0, S, sig2_0, reg, lb_sig2, ub_sig2, parameters);
rank = count_rank(flip(eig(Newton_out.T)), eval_factor);

if ~adaptive_reg
    out.x = Newton_out.x; out.y = Newton_out.y; out.sig2 = Newton_out.sig2;
    out.T = Newton_out.T; out.reg = reg;
    if parameters.verbose == 2
        fprintf("chosen reg:  %f\n", reg)
    end
    return
end



interval = 0;
if rank < target_rank - interval
    lb = reg/2;
    while true
        Newton_out = Newton_ub_lb(x0, y0, S, sig2_0, lb, lb_sig2, ub_sig2, parameters);
        rank = count_rank(flip(eig(Newton_out.T)), eval_factor);
        if rank >= target_rank
            break
        else
            lb = lb/2;
        end
    end
elseif rank > target_rank + interval
    ub = reg*2;
    while true
        Newton_out = Newton_ub_lb(x0, y0, S, sig2_0, ub, lb_sig2, ub_sig2, parameters);
        rank = count_rank(flip(eig(Newton_out.T)), eval_factor);
        if rank <= target_rank
            break
        else
            ub = ub*2;
        end
    end
 else
    out.x = Newton_out.x; out.y = Newton_out.y; out.sig2 = Newton_out.sig2;
    out.T = Newton_out.T; out.reg = reg;
    if parameters.verbose == 2
        fprintf("chosen reg:  %f\n", reg)
    end
    return;
 end

while true 
    if parameters.verbose == 1
        fprintf('testing reg %d\n', reg)
    end

    Newton_out = Newton_ub_lb(x0, y0, S, sig2_0, reg, lb_sig2, ub_sig2, parameters);
    rank = count_rank(flip(eig(Newton_out.T)), eval_factor);
    
    if parameters.verbose == 1
        fprintf('rank %d\n', rank)
    end

    if abs(rank-target_rank) <= interval
        break
    elseif rank < target_rank
        ub = reg;
        reg = (ub + lb)/2;
    elseif rank > target_rank
        lb = reg;
        reg = (ub + lb)/2;
    end
    if abs(ub-lb) < 1e-6
        warning('Rank does not match target rank')
        break
    end
  
end

if parameters.verbose == 2
    fprintf("chosen reg:  %f\n", reg)
end

out.x = Newton_out.x;
out.y = Newton_out.y;
out.sig2 = Newton_out.sig2;
out.T = Newton_out.T;
out.reg = reg;
     
end



