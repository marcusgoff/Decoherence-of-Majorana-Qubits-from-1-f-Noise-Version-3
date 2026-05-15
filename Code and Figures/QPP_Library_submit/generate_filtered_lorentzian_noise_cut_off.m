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
% NOT INCLUDING MAGNITUDE OF NOISE
%
% INPUTS:
%      cut_off_freq - in angular frequency


function [lorentzian_noise, t_noise, white_noise, S_w, omega_vec, Y, Y_filt] = ...
    generate_filtered_lorentzian_noise_cut_off(t_init, t_sim, delta_t, t_noise_factor,...
    sigma, tau, signal_power, cut_off_freq)

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
    
    % Get Lorentzian Power Spectrum
    S_w = get_lorentzian_PSD(omega_vec, sigma, tau, signal_power);
    %S_w = (1/pi).*(sigma*tau./(sigma+tau).^2).*...
    %    (sigma^-1 + tau^-1)./(omega_vec.^2 + (sigma^-1+tau^-1).^2);
    %S_w = 2*pi*S_w;
    %S_w = 4.*signal_power.*S_w;  

    Y_filt = sqrt(S_w/delta_t).*Y;  
    
    % Cut-off frequencies
    cut_off_ind = omega_vec >= abs(cut_off_freq) | ...
        omega_vec <= -abs(cut_off_freq);  
    Y_filt(cut_off_ind) = 0; 

    lorentzian_noise = ifft(fftshift(Y_filt));
    %lorentzian_noise = 2.*sqrt(signal_power).*lorentzian_noise;

end






