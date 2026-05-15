%--------------------------------%
% Generate Gaussian Noise        %
%--------------------------------%
% I don't trust this. I directly use Mishmash's equations, then I get an
% imaginary time series, because they don't allow for negative frequencies.
% I've butchered their equations to try and get something to work here, but
% I think it's best to - at least at present - just use my
% generate_gaussian_noise function. 
%
% This Function generates gaussian noise. 
% Would be worthwile to make a general noise generation function for
% different S(w) functions. 
%
% For some simulation up to some finite time T_sim, we generate a noise
% time series for T_noise >> T_sim. 
%
% Outputs from output number 2 are more for debugging. 
% UPDATES MADE:
%       01/05/24: Following suit with the definitions given in Mishmash's appendices
%                 I have removed ffshift. 
%                 Redefined angular frequency vector, omega_vec in terms of
%                 only "positive frequencies" as done in Mishmash. Note the
%                 second half of omega_vec is technically negative
%                 frequencies if you want to compare this to a continuous
%                 time FT. 
%                 Also added code to ensure that N is even. - TO DO


function [gaussian_noise, t_noise, white_noise, gaussian_spectrum, omega_vec, Y, Y_filt] = ...
    generate_gaussian_noise_mishmash_definition(t_init, t_sim, delta_t, t_noise_factor,...
    fluctuation_amplitude, characteristic_noise_corr_freq)

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
    Y = fft(white_noise);
    
    % Define the omega_k vector.
    N = T_max/delta_t;
    omega_sampling = 2*pi./delta_t; % NOT IN USE
    omega_vec_o = linspace(-omega_sampling/2, omega_sampling/2, N+1); % NOT IN USE
    omega_vec_o = omega_vec_o(1:end-1); % NOT IN USE
    
    %omega_vec = (2*pi./T_max).*(0:1:(N-1));
    omega_vec = (2*pi./T_max).*(-floor(N/2):1:(-floor(N/2))+N-1);
    
    % Get Gaussian Power Spectrum
    gaussian_spectrum = get_gaussian_power_spectrum(omega_vec, fluctuation_amplitude, characteristic_noise_corr_freq); 
    Y_filt = fftshift(sqrt(gaussian_spectrum/delta_t).*fftshift(Y));
    
    
    gaussian_noise = ifft(fftshift(Y_filt));
    
end






