addpath('../algorithm/', '../algorithm/other_algorithms/', '../utils/');
clear; clc; rng(1);
theta_rad = [-10, -5, 0, 5, 10]*pi/180; M = length(theta_rad);
power_source = 1; P = power_source*eye(M); wavelength = 1; d = wavelength/2;
m = 15; reg0 = 0.5;

snr = 5;       
sig2_0 = power_source*10^(-snr/10);            
samples = (15:5:50); MC_runs = 300;


% Containers for evaluating the performance
MSE_SC = zeros(1, length(samples)); MSE_low_rank_cov = zeros(1, length(samples)); 
MSE_DA = zeros(1, length(samples)); MSE_LRSE = zeros(1, length(samples));
MSE_AML = zeros(1, length(samples)); MSE_ITAM = zeros(1, length(samples)); 
MSE_low_rank_zero_reg_cov = zeros(1, length(samples));
crb_sto = zeros(1, length(samples));
crb_sto_uc = zeros(1, length(samples));

for ii = 1:length(samples)
    K = samples(ii);
    fprintf("Simulating K = %i \n", K)
    
    for run = 1:MC_runs 
        % Generate data
        [Y, true_cov] = generate_ula_data(power_source, sig2_0, d, m, M, K, ...
              wavelength, theta_rad);

        % Estimate covariance matrix.
        [sample_cov, DA_cov, AML_cov, LRSE_cov, ITAM_cov, lowrankToep_out, ...
            lowrankToep_zero_reg_out] =  estimates_cov(Y, sig2_0, M, reg0);
        low_rank_cov = lowrankToep_out.cov;
        reg0 = lowrankToep_out.reg;
        low_rank_cov_zero_reg = lowrankToep_zero_reg_out.cov;

        % Run root-music.
        doa_sample_cov = rmusic_1d(sample_cov, M, 2*pi*d/wavelength);
        doa_DA_cov = rmusic_1d(DA_cov, M, 2*pi*d/wavelength);
        doa_AML_cov = rmusic_1d(AML_cov, M, 2*pi*d/wavelength);
        doa_LRSE_cov = rmusic_1d(LRSE_cov, M, 2*pi*d/wavelength);
        doa_ITAM_cov = rmusic_1d(ITAM_cov, M, 2*pi*d/wavelength);
        doa_low_rank_cov = rmusic_1d(low_rank_cov, M, 2*pi*d/wavelength);
        doa_low_rank_zero_reg_cov = rmusic_1d(low_rank_cov_zero_reg, M, 2*pi*d/wavelength);
        
        % Compensate for different convention on angles.
        doa_sample_cov = sort(-doa_sample_cov.x_est);
        doa_DA_cov = sort(-doa_DA_cov.x_est);
        doa_AML_cov = sort(-doa_AML_cov.x_est);
        doa_LRSE_cov = sort(-doa_LRSE_cov.x_est);
        doa_ITAM_cov = sort(-doa_ITAM_cov.x_est);
        doa_low_rank_cov = sort(-doa_low_rank_cov.x_est);
        doa_low_rank_cov_zero_reg = sort(-doa_low_rank_zero_reg_cov.x_est);

        % Compute angle estimation error
        MSE_SC(ii) = MSE_SC(ii) + norm(theta_rad - doa_sample_cov)^2;
        MSE_DA(ii) = MSE_DA(ii) + norm(theta_rad - doa_DA_cov)^2;
        MSE_AML(ii) = MSE_AML(ii) + norm(theta_rad - doa_AML_cov)^2;
        MSE_LRSE(ii) = MSE_LRSE(ii) + norm(theta_rad - doa_LRSE_cov)^2;
        MSE_ITAM(ii) = MSE_ITAM(ii) + norm(theta_rad - doa_ITAM_cov)^2;
        MSE_low_rank_cov(ii) = MSE_low_rank_cov(ii) + norm(theta_rad - doa_low_rank_cov)^2;
        MSE_low_rank_zero_reg_cov(ii) = MSE_low_rank_zero_reg_cov(ii) + norm(theta_rad - doa_low_rank_cov_zero_reg)^2;

    end  
    MSE_SC(ii) = MSE_SC(ii)/(M*MC_runs);
    MSE_DA(ii) = MSE_DA(ii)/(M*MC_runs);
    MSE_AML(ii) = MSE_AML(ii)/(M*MC_runs);
    MSE_LRSE(ii) = MSE_LRSE(ii)/(M*MC_runs);
    MSE_ITAM(ii) = MSE_ITAM(ii)/(M*MC_runs);
    MSE_low_rank_cov(ii) = MSE_low_rank_cov(ii)/(M*MC_runs); 
    MSE_low_rank_zero_reg_cov(ii) = MSE_low_rank_zero_reg_cov(ii)/(M*MC_runs); 
    crb_sto(ii) = mean(diag(U_CRB(P, theta_rad, sig2_0, m, d, K)));
    crb_sto_uc(ii) = mean(diag(S_CRB(P, theta_rad, sig2_0, m, d, K)));
end

%%
figure()
% Plot the data with specified markers and line styles
h1 = semilogy(samples, MSE_SC, '-o'); hold on;
h2 = semilogy(samples, MSE_DA, '-x');
h3 = semilogy(samples, MSE_AML, '*-');
h4 = semilogy(samples, MSE_LRSE, '-v');
h5 = semilogy(samples, MSE_ITAM, '-+');
h6 = semilogy(samples, MSE_low_rank_cov, '-d');
h7 = semilogy(samples, MSE_low_rank_zero_reg_cov, '-s');
h8 = semilogy(samples, crb_sto, '--');
h9 = semilogy(samples, crb_sto_uc, '--');

% Set the colors for each plot
% Set the colors for each plot according to your specification
set(h1, 'Color', [0 0.4470 0.7410]); % Blue
set(h2, 'Color', [0.4940 0.1840 0.5560]); % Purple
set(h3, 'Color', [0.4660 0.6740 0.1880]); % Green
set(h4, 'Color', [0.3010 0.7450 0.9330]); % Light Blue (Cyan)
set(h5, 'Color', [0.6350 0.0780 0.1840]); % Red
set(h6, 'Color', [0.8500 0.3250 0.0980]); % Orange
set(h7, 'Color', [0.9290 0.6940 0.1250]); % Yellow
set(h8, 'Color', [0 0.4470 0.7410]); % Blue
set(h9, 'Color', [0.8500 0.3250 0.0980]); % Orange

hold off;
ylim([1e-4, 0.3]);
grid on;
ylabel('MSE', 'Interpreter', 'Latex', 'fontsize', 12);
xlabel('$K$', 'Interpreter', 'Latex', 'fontsize', 12)
legend('SC', 'DA', 'AML', 'LRSE', 'ITAM', 'NLRT', 'NT', 'U-CRB', 'S-CRB','Location', 'northoutside','NumColumns', 5, 'Orientation', 'horizontal');

