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



function [lorentzian_noise, t_noise, white_noise, S_w, omega_vec, Y, Y_filt] = ...
    generate_filtered_lorentzian_noise(t_init, t_sim, delta_t, t_noise_factor,...
    sigma, tau, signal_power)

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
    S_w = (1/pi).*(sigma*tau./(sigma+tau).^2).*...
        (sigma^-1 + tau^-1)./(omega_vec.^2 + (sigma^-1+tau^-1).^2);
    %S_w(omega_vec == 0) = S_w(omega_vec == 0) + (sigma/(sigma+tau)).^2;
    S_w = 2*pi*S_w;

%    S_w = (1/pi).*sigma*tau/(sigma + tau).^2.*...
%        (sigma^-1 + tau^-1)./(omega_vec.^2+(sigma^-1 + tau^-1).^2);
%    S_w(omega_vec ==0) = S_w(omega_vec ==0) + (sigma/(sigma+tau)).^2;

    Y_filt = sqrt(S_w/delta_t).*Y;    
    
    lorentzian_noise = ifft(fftshift(Y_filt));
    lorentzian_noise = 2.*sqrt(signal_power).*lorentzian_noise;

end






