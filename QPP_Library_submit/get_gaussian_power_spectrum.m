%--------------------------------%
% Get Gaussian Power Spectrum
%--------------------------------%
% Returns the power spectrum function for Gaussian noise. This corresponds
% to equation 31 in Mishmas, Oppen, Alicea, Phys Rev. B 101 075404 (2020).


function S = get_gaussian_power_spectrum(omega_vec, fluctuation_amplitude,...
    characteristic_noise_corr_freq)

    D = fluctuation_amplitude;
    kappa = characteristic_noise_corr_freq;

    S = D.^2*sqrt(4*pi/kappa.^2)*exp(-(omega_vec./kappa).^2);
    

end