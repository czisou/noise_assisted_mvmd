%% Main Script for Noise-Assisted Multivariate Variatonal Mode Decomposition
%**************************************************************************
% Written by:
% Charilaos A. Zisou, charilaz@auth.gr
% Georgios K. Apostolidis, gkaposto@auth.gr
% Tested with MATLAB R2018a
%**************************************************************************
close all; clear; clc;

% Initialize constants
n_repeats = 1000;
n_atoms = 3;
snr_vals = [20 15 10];
lengths = [500; 125; 100];
centers = [250; 125; 375];
params = zeros(n_atoms, 4);
L = 500;
thresh = 0.5;
%% Simulation loop
for i = 1:length(snr_vals)
    snr = snr_vals(i);
    for r = 1:n_repeats
        % Generate random frequencies and amplitudes
        freqs = sort(rand(n_atoms, 1));
        amps  = 0.5 + rand(n_atoms, 1);
        
        % Create parameter matrix
        params = [centers, freqs*pi, lengths, amps];
        
        % Generate simulation signal
        [x, x_c, ~] = generate_signal(L, params, snr);
        
        % Run both algorithms for different K (number of output modes)
        for k = n_atoms:10
            % Original VMD
            [u_vmd, ~, omega] = VMD(x, 1000, 0, k, 0, 1, 1e-8);
            
            % Sort the output modes based on their central frequencies
            [~, sortIndex] = sort(omega(end,:));
            omega = omega(:,sortIndex);
            u_vmd = u_vmd(sortIndex, :);
            
            % Noise-Assisted MVMD
            [u_namvmd, ~, omega] = na_mvmd(x, 1000, 0, k, 1, 1, 1e-8, 1);
            u_namvmd = squeeze(u_namvmd(:, 1, :))';
            
            % Sort the output modes based on their central frequencies
            [~, sortIndex] = sort(omega(end,:));
            omega = omega(:,sortIndex);
            u_namvmd = u_namvmd(sortIndex, :);
            
            % Determine if the decomposition was successful or not
            success = is_successful(u_vmd, x_c, x, thresh);
            vmd_success(r, k, i) = success;
            
            success = is_successful(u_namvmd, x_c, x, thresh);
            namvmd_success(r, k, i) = success;
            
            % Calculate the orthogonality index
            vmd_orth(r, k, i) = orthogonality_index(u_vmd, x);             
            namvmd_orth(r, k, i) = orthogonality_index(u_namvmd, x);          
        end
    end
end
%% Create plots
% Success rate plots
vmd_success_rate = squeeze(mean(vmd_success, 1));
namvmd_success_rate = squeeze(mean(namvmd_success, 1));

figure;
plot(vmd_success_rate)
title('Success rate with respect to K for VMD')
ylabel('Success rate')
xlabel('Number of modes K')
legend('SNR=20dB', 'SNR=15dB', 'SNR=10dB', 'Location', 'northwest')

figure;
plot(namvmd_success_rate)
title('Success rate with respect to K for NA-MVMD')
ylabel('Success rate')
xlabel('Number of modes K')
legend('SNR=20dB', 'SNR=15dB', 'SNR=10dB', 'Location', 'northwest')

% Orthogonality index (OI)

% Compute OI in case of successful decomposition, for VMD
vmd_orth = vmd_success .* vmd_orth;
vmd_orth = squeeze(vmd_orth(:, 3, :));
vmd_orth(vmd_orth==0) = NaN;
vmd_std = std(vmd_orth,'omitnan');

% Remove outliers
vmd_orth(vmd_orth>=3*vmd_std) = NaN;
vmd_mean = nanmean(vmd_orth);
vmd_std = std(vmd_orth,'omitnan');

% Compute OI in case of successful decomposition, for NA-MVMD
namvmd_orth = namvmd_success .* namvmd_orth;
namvmd_orth = squeeze(namvmd_orth(:, 8, :));
namvmd_orth(namvmd_orth==0) = NaN;
namvmd_std = std(namvmd_orth,'omitnan');

% Remove outliers
namvmd_orth(namvmd_orth>=3*namvmd_std) = NaN;
namvmd_mean = nanmean(namvmd_orth);
namvmd_std = std(namvmd_orth,'omitnan');

% Create stat matrix for display purposes
stat_mat = ["NA-MVMD" 100*reshape([namvmd_mean; namvmd_std], [], 1)'; ...
            "VMD    " 100*reshape([vmd_mean; vmd_std], [], 1)'].';
        
% Print matrix
disp('        |       Orthogonality Index (OI) x 10^2      |')
formatSpec = '%s | %1.2f +- %1.2f | %1.2f +- %1.2f | %1.2f +- %1.2f | \n';
disp('SNR     |    20 dB     |    15 dB     |    10 dB     |')
fprintf(formatSpec, stat_mat)