%--------------------------------------%
% Charge Noise Engine - General - SHORT
%--------------------------------------%
% This is a shorter version of charge_noise_engine_inst_basis_setup. With
% some rarely used functionality removed. 
%
%
% INPUTS:
%   noise_type: A string denoting the noise type. Allowed inputs:
%               'gaussian', 'two_state_exp', or 'custom'
%   noise_parameters: An array or matrix, depending on the the entry on the
%               noise_type. For: 
%               Gaussian: [t_noise_factor, fluctuation_amplitude, characteristic_noise_corr_freq]
%                       two_state_exp: [sigma, tau, fluctuation_amplitude]
%               custom: matrix of size "num_time_steps x 2*num_trials", whose
%                       columns contain independent noise trajectories. The
%                       first num_trial columns are for the top kitaev
%                       chain and the following num_trial columns are for
%                       the bottom Kitaev chain. 
%
%   time_evolution_method: Choose method for obtaining the unitary U_t generating
%                          time evolution at each time step, from the noisy
%                          Hamiltonian, H_t, at that time. 
%               'time_ev_diagonalise' - Diagonalise H_t to construct U_t.
%               'time_ev_exp' - Use matlab function "expm" to obtain U_t
%               'time_ev_1st_order' - Use a 1st order Taylor series
%                              expansion to obtain U_t.
%               'time_ev_Magnus' - use the Magnus expansion - TO BE IMPLEMENTED.  
%
%   save_output = 'save_output' or 'do_not_save_output'
%
%   plot_results = 'plot_all_results', 'plot_X_exp_only' or 'do_not_plot_results'
%
%   noise_trajectories = 1st element of varargin. 2*num_trial x
%       num_time_steps_noise matrix, where each row is a noise realisation. 
%       num_time_steps_noise = (t_final-t_init)/delta_t. First num_trial rows
%       are for the top Kitaev chain, the next num_trial rows are for the 
%       bottom Kitaev chain. If this input is not provided, then a new noise 
%       trajectory is generated. 
%
%   init_condition - Select whether at t=0 the system is initialised 
%                       according to the clean Hamiltonian or the first time 
%                       step of the noise-trajectory. 
%                       Allowed values are:
%                       'initial_state_no_noise', 'initial_state_noise'.
%
%   basis_choice - Select whether to define the Pauli operators based of the 
%                  Hamiltonian at t = 0 or the instanteneous Hamiltonian
%                  H(t). Allowed values are:
%                  'init_hamil_basis', 'instant_hamil_basis'.
%   delta_t_vec - either a scalar or vector of delta_t values. 
%
% OUTPUTS
%   e_vals_mat - this contains the eigenvalues only for the first trajectory. 
%                I don't want to pass all trajectory's eigenvalues out. 
%                The (i, j)'th element corresponds to the j'th eigenvalue
%                at the t_i time step. 
%
% To do - change above plot parameter to distinguish, the full plots of
% individual trajectories and everythign, or just the X_exp mean results
% and fit.
%
% INIT Eigenvectors are both the initial eigenvectors and the "averaged" Hamiltonian
% % eigenvectors - ie. with no noise.
% 
    
function [X_exp, Z_exp, Y_exp, ...
         P_y1y2y3y4, overlap_with_init_state, overlap_fixed_basis_states, ...
         t_vec, comparison_determinants_mat, quasi_occ_t, ...
         overlap_inst_basis_states] = ...
    charge_noise_engine_ibs_short(mu_mean, mu_offset, w, delta, N, BC, t_init, t_final,...
     delta_t, noise_parameters, time_evolution_method, init_condition, basis_choice, init_state)  

     %% Process Input Parameters
     mu = mu_mean + mu_offset*[+1, -1];
    
     %% Get Initial Hamiltonian - if defined by clean Hamiltonian
    if strcmp(init_condition, 'initial_state_no_noise')
        H_tetron_init = get_tetron_BdG_Hamiltonian(mu, w, delta, N, BC, 'chain_1_chain_2');
        [e_vecs_init, e_vals_init, majorana_zero_modes_ref] = ... 
            diagonalise_uncoupled_tetron_via_Kitaev_chains(H_tetron_init(1:2*N,1:2*N), ...
            H_tetron_init((2*N+1):end, (2*N+1):end), 'dirac_zero_modes');
    
        [~, ~, corr_qp_init] = ...
            get_tetron_init_cov_mat_general(e_vecs_init, init_state); 

    end

    %% Run Simulations
    if isscalar(delta_t)
        t_vec = t_init:delta_t:t_final;
    else
        t_vec = t_init + [t_init, cumsum(delta_t(:)).'];
        if abs((t_vec(end) - t_final)./t_final) >10*eps
            warning('delta_t vector does not match t_init and t_final');
            t_final = t_vec(end);
        end
    end

    tStart = tic;

    %% Get Noise Trajectories for top and bottom chains
    noise_top = noise_parameters(:, 1).';
    noise_bot = noise_parameters(:, 2).';

    % Store the top/bot noise trajectories in a matrix
    mu_vec = zeros(2, length(t_vec));
    mu_vec(1,:) = mu_mean + mu_offset +  ...
        noise_top;
    
    mu_vec(2,:) = mu_mean - mu_offset +  ...
        noise_bot;

    % Get initial Hamiltonian - if defining by first step of noise
    % trajectory
    if strcmp(init_condition, 'initial_state_noise')
        H_tetron_init = get_tetron_BdG_Hamiltonian(mu_vec(:,1), w, delta, N, BC, 'chain_1_chain_2');
        [e_vecs_init, e_vals_init, majorana_zero_modes_ref] = ... 
            diagonalise_uncoupled_tetron_via_Kitaev_chains(H_tetron_init(1:2*N,1:2*N), ...
            H_tetron_init((2*N+1):end, (2*N+1):end), 'dirac_zero_modes');
      
        [~, ~, corr_qp_init] = ...
             get_tetron_init_cov_mat_general(e_vecs_init, init_state);     
    end

    % Do time evolution
    tetron_params = {w, delta, N, BC};

    % Dodgey handling of the method I used to calculate
    % overlap_inst_basis_states
    [Z_exp, X_exp, Y_exp, P_y1y2y3y4, ~, overlap_with_init_state, ...
        overlap_fixed_basis_states, ~,~,~, comparison_determinants_mat, ...
        e_vals_mat_curr, ~, quasi_occ_t, ...
        overlap_inst_basis_states] =...
        calc_noisy_tetron_ev_inst_basis_setup(...
        corr_qp_init, e_vecs_init, e_vals_init, tetron_params,...
        mu_vec, delta_t, t_init, t_final, time_evolution_method, ...
        basis_choice, majorana_zero_modes_ref, mu);


    %comparison_determinants_mat_top(trial_no, :) = comparison_determinants_mat(1, :);
    %comparison_determinants_mat_bot(trial_no, :) = comparison_determinants_mat(2, :);

    % Z_exp_mat(trial_no, :) = Z_exp;
    % X_exp_mat(trial_no, :) = X_exp;
    % Y_exp_mat(trial_no, :) = Y_exp;
    % P_y1y2y3y4_v1_mat(trial_no, :) = P_y1y2y3y4_v1;
    % P_y1y2y3y4_v2_mat(trial_no, :) = P_y1y2y3y4_v2;        
    % overlap_with_init_state_mat(trial_no, :) = overlap_with_init_state;
    % overlap_fixed_basis_states_mat(trial_no, :) = overlap_fixed_basis_states;
    % overlap_inst_basis_states_mat(trial_no, :) = overlap_inst_basis_states;        
    % 


 
    %% Save Output

    % Format a nice string which gives all the parameters and saves them in
    % this cell. 

end


