%-------------------------------------------------%
% Generate Gaussian with Lorentzian Tail Noise
%-------------------------------------------------%
% The idea here is to generate noise that is an interpolation between
% Gaussian and Lorentzian noise
% I ensure the combined Gaussian-Lorentzian noise PSD has power of
% signal_power_total, however I do not renormalise the power after applying
% the cut-off. 


function [gauss_lor_noise, t_noise, S, S_gauss, S_lor, omega_vec, white_noise, Y, Y_filt] = ...
    generate_gauss_lor_tail_noise(t_init, t_sim, delta_t, t_noise_factor,...
    sigma, tau, signal_power_lor, signal_power_gauss, noise_corr_freq_gauss, ...
    signal_power_total, cut_off_freq)

    % CHECK INPUTS
    t_noise_final = t_init + (t_sim - t_init)*t_noise_factor;
    T_max = t_noise_final - t_init;
    t_noise = t_init:delta_t:t_noise_final;
    t_noise = t_noise(1:end-1);
    
    N = T_max/delta_t;
    omega_sampling = 2*pi./delta_t;

    % Get white noise time series
    mu_white = 0; sigma_white = 1; 
    white_noise = normrnd(mu_white, sigma_white, size(t_noise));
    
    % Discrete Fourier Transform 
    omega_vec = linspace(-omega_sampling/2, omega_sampling/2, N+1);
    omega_vec = omega_vec(1:end-1);
    Y = fftshift(fft(white_noise));
    
    % Get Power Spectrum - package this in a function
    fluctuation_amplitude_gauss = sqrt(signal_power_gauss);
    S_gauss = get_gaussian_power_spectrum(omega_vec, fluctuation_amplitude_gauss,...
    noise_corr_freq_gauss);
    S_lor = get_lorentzian_PSD(omega_vec, sigma, tau, signal_power_lor);
    S = max([S_gauss(:), S_lor(:)], [],2).';

    signal_power_curr = 1./(2*pi).*sum(S).*(omega_vec(2)-omega_vec(1));
    S = S.*(signal_power_total./signal_power_curr);
    
    S_gauss = S_gauss.*(signal_power_total./signal_power_curr);
    S_lor = S_lor.*(signal_power_total./signal_power_curr);

    % Cut-off frequencies
    cut_off_ind = omega_vec >= abs(cut_off_freq) | ...
        omega_vec <= -abs(cut_off_freq);  
    S(cut_off_ind) = 0; 

    Y_filt = sqrt(S/delta_t).*Y;  

    gauss_lor_noise = fftshift(Y_filt);
    gauss_lor_noise = ifft(gauss_lor_noise);
    
    %gauss_lor_noise = 2.*sqrt(signal_power).*gauss_lor_noise;
    t_vec = t_init:delta_t:t_sim;
    gauss_lor_noise = gauss_lor_noise(1:length(t_vec));

    % Noramlise to get the desired total signal power
   % signal_power_curr = 1./(2*pi).*sum(S).*(omega_vec(2)-omega_vec(1));
   % gauss_lor_noise = gauss_lor_noise.*sqrt(signal_power_total./signal_power_curr);

end

