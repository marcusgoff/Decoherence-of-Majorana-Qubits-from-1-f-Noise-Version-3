%-------------------------------------%
% Get Lorentzian Power Spectral Density
%-------------------------------------%
% Get PSD function for Lorentzan function. 
% Used for filtering white noise. Note: you will need to scale this by
% 1/delta_t when multiplying by the fft of a time domain signal. 
% 

function lorentzian = get_lorentzian_PSD(omega_vec, sigma, tau, signal_power)

    lorentzian = (1/pi).*(sigma*tau./(sigma+tau).^2).*...
        (sigma^-1 + tau^-1)./(omega_vec.^2 + (sigma^-1+tau^-1).^2);
    lorentzian = 2*pi*lorentzian;
    lorentzian = 4.*signal_power.*lorentzian;  

end 