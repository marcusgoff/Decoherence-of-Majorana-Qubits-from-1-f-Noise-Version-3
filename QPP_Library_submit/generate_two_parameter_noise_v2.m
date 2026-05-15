%--------------------------------------------%
%  Generate Two-Parameter Random Noise
%--------------------------------------------%
% Following the Method of Machlup: J. Appl. Phys. 25, 341–343 (1954).
%
% INPUTS
%   sigma - 1 mean lifetime
%   tau   - 0 mean life time
%


function [state_vec, t_vec ,jump_ind] = ...
        generate_two_parameter_noise_v2(sigma, tau, t_init, t_final, delta_t, signal_power, remove_dc, init_state_str)
    %% Setup
    t_vec = (t_init:delta_t:t_final).';
    N = length(t_vec);
    jump_ind = zeros(size(t_vec)); % vector with time steps at the time of each jump (these 
                  % are the times at which I'll evaluate the H evolved
                  % simulation)

    state_vec = zeros(N,1);

    % Get the initial state
    switch init_state_str
        case 'init_low'
            init_state = 0;
        case 'init_high'
            init_state = 1;
        case 'init_random'
            random_num = rand(1); 
            if random_num < sigma/(sigma+tau)
                init_state = 1;
            else
                init_state = 0; 
            end
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
            jump_ind(ii+1) = 1;
        else 
            state_vec(ii+1) = state_vec(ii);
        end
       
    end
    jump_ind = jump_ind ==1; %Make elements logical values

    switch remove_dc
        case 'remove_dc'
            state_vec = state_vec - 0.5;
            state_vec = 2*sqrt(signal_power).*state_vec;
        case 'do_not_remove_dc'
            % Not sure if this is correct in the spirit of Machlup, but
            % this is what I need:
            state_vec = 2*sqrt(signal_power).*state_vec;
        otherwise
            error('Please enter command for remove_dc in generate_two_parameter_noise');
end