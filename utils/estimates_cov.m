function [sample_cov, DA_cov, AML_cov, LRSE_out, ITAM_cov, lowrankToep_out, ...
         lowrankToep_zero_reg_out] =  estimates_cov(Y, sig2_0, target_rank, reg0)

[n, K] = size(Y);
n = n - 1;

lb_sig2 = 0.5*sig2_0;
ub_sig2 = 2.0*sig2_0;

% sample covariance
sample_cov = 1/K*(Y*Y');

% apply diagonal averaging
[DA_cov, ~, ~] = AAD(sample_cov);

% apply AML
if K >= n+1
    [Psi_AML] = create_Psi_AML(size(sample_cov, 1));
    [AML_cov] = AML(sample_cov, K, Psi_AML);
    AML_cov = AML_cov.estimate;
else
    AML_cov = -1;
end

% apply LRSE.
LRSE_out = LRSE_LB_UB(Y, lb_sig2, ub_sig2, target_rank);

%apply ITAM
ITAM_cov = ITAM(sample_cov, target_rank, 50);

%apply lowrankToep
init_strategy = 2; [x0, y0] = initialization(sample_cov, Y, init_strategy);
parameters.verbose = 0; parameters.max_iter = 200;
parameters.tol = 1e-7; parameters.beta = 0.8; parameters.alpha = 0.05;
parameters.t = 10; parameters.mu = 15;

Newton_out = lowrankToep(x0, y0, sample_cov, sig2_0, reg0, parameters, ...
                         target_rank, lb_sig2, ub_sig2, true);
lowrankToep_out.cov = Newton_out.T + Newton_out.sig2 * eye(n + 1);
lowrankToep_out.reg = Newton_out.reg;
lowrankToep_out.sig2 = Newton_out.sig2;

parameters.verbose=0;
Newton_out = lowrankToep(x0, y0, sample_cov, sig2_0, 0.0, parameters, ...
                         target_rank, lb_sig2, ub_sig2, false);
lowrankToep_zero_reg_out.cov = Newton_out.T + Newton_out.sig2 * eye(n + 1);
lowrankToep_zero_reg_out.sig2 = Newton_out.sig2;

end