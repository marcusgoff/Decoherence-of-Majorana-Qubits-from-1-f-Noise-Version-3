%--------------------------------%
% Generate Filtered Noise        %
%--------------------------------%
% This Function generates noise by filtering white noise in the frequency
% domain with the input function S. 
%
% For some simulation up to some finite time T_sim, we generate a noise
% time series for T_noise >> T_sim. 
%
% Outputs from output number 2 are more for debugging. 
%
% INPUTS:
%       S_w - (angular) freqency domain power spectral density

% NOT COMPLETED 


function [gaussian_noise, t_noise, white_noise, gaussian_spectrum, omega_vec, Y, Y_filt] = ...
    generate_gaussian_noise(t_init, t_sim, delta_t, t_noise_factor,...
    S_w)

    % CHECK INPUTS
    t_noise_final = t_init + (t_sim - t_init)*t_noise_factor;
    T_max = t_noise_final - t_init;
    t_noise = t_init:delta_t:t_noise_final;
    t_noise = t_noise(1:end-1);
    
    N = T_max/delta_t;
    omega_sampling = 2*pi./delta_t;

    % Get white noise time series
    mu = 0; sigma = 1; 
    white_noise = normrnd(mu, sigma, size(t_noise));
    
    % Discrete Fourier Transform 
    omega_vec = linspace(-omega_sampling/2, omega_sampling/2, N+1);
    omega_vec = omega_vec(1:end-1);
    Y = fftshift(fft(white_noise));
    
    % Get Gaussian Power Spectrum
    gaussian_spectrum = get_gaussian_power_spectrum(omega_vec, fluctuation_amplitude, characteristic_noise_corr_freq); 
    Y_filt = sqrt(gaussian_spectrum/delta_t).*Y;
    
    
    gaussian_noise = ifft(fftshift(Y_filt));

end






