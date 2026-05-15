%--------------------------------------------%
%  Generate Two-Parameter Random Noise
%--------------------------------------------%
% Following the Method of Machlup: J. Appl. Phys. 25, 341–343 (1954).
%
% INPUTS
%   sigma - 1 mean lifetime
%   tau   - 0 mean life time
%


function [state_vec, t_vec] = ...
        generate_two_parameter_noise(sigma, tau, t_init, t_final, delta_t, signal_power, remove_dc)
    %% Setup
    t_vec = (t_init:delta_t:t_final).';
    N = length(t_vec);
    
    state_vec = zeros(N,1);

    % Get the initial state
    random_num = rand(1); 
    if random_num < sigma/(sigma+tau)
        init_state = 1;
    else
        init_state = 0; 
    end
    state_vec(1) = init_state;
    
    %% Run simulation
    
    for ii = 1:(N-1)
        if state_vec(ii) == 0
            trans_prob = delta_t/tau;
        elseif state_vec(ii) == 1
            trans_prob = delta_t/sigma;
        else
            error('Invalid State');
        end
        
        random_num = rand(1);
        if trans_prob > random_num
            state_vec(ii+1) = mod(state_vec(ii) + 1,2);
        else 
            state_vec(ii+1) = state_vec(ii);
        end
       
    end

    switch remove_dc
        case 'remove_dc'
            state_vec = state_vec - 0.5;
            state_vec = 2*sqrt(signal_power).*state_vec;
        case 'do_not_remove_dc'
            state_vec = sqrt(signal_power).*sqrt((sigma + tau)/sigma); % this is wrong
        otherwise
            error('Please enter command for remove_dc in generate_two_parameter_noise');
end