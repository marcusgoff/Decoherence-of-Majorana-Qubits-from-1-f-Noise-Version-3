%--------------------------------------%
% Charge Noise Engine - General - SHORT
%--------------------------------------%
% Adapated from charge_noise_engine_ibs_short for specifically a two level system. 
% This will either only diagonalise the upper and lower states' Hamiltonians 
% once, each. Or alternatively the user may provide the diagonalised
% Hamiltonians. Doing so will drastically decrease the run times. 
%
% Requires the initial state to be "low"
%
% Note: I deleted the option "initial_state_noise" which was left over from
% the source of this function.
%
% IMPORTANT: for this function "mu_vec" replaces noise_traj in previous
% interations. Here "mu_vec" is just the vector of mu values (not
% delta_mu values) - with no offset between the chains.  
%
% INPUTS:
%   noise_type: A string denoting the noise type. Allowed inputs:
%               'gaussian', 'two_state_exp', or 'custom'
%   noise_parameters: A single noise trajectory of size "num_time_steps x 1". 
%
%   time_evolution_method: Set to 'time_ev_diagonalise'.  
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
         comparison_determinants_mat, quasi_occ_t, ...
         overlap_inst_basis_states] = ...
    charge_noise_engine_ibs_tls(mu_low, mu_high, mu_offset, w, delta, N, BC, t_init, t_final,...
     delta_t, mu_vec_no_off, basis_choice, init_state)  

   
     %% Process Inputs
    mu_vec = mu_vec_no_off(:).' + mu_offset*[+1; -1];
    if size(mu_vec_no_off, 2) > 1
        error('For this function, mu_vec_no_off must be a single column vector');
    end

    if mu_vec_no_off(1) ~= mu_low
        error('Initial State must be mu_low');
    end
    %% Run Simulations
    % if isscalar(delta_t)
    %     t_vec = t_init:delta_t:t_final;
    % else
    %     t_vec = t_init + [t_init, cumsum(delta_t(:)).'];
    %     if abs((t_vec(end) - t_final)./t_final) >1000*eps
    %         %warning('delta_t vector does not match t_init and t_final');
    %         t_final = t_vec(end);
    %     end
    % end

    tStart = tic;

    %% Get Noise Trajectories for top and bottom chains

    mu_val_1 = mu_vec(:,1);
    
    % Get initial Hamiltonian
    H_tetron_1 = get_tetron_BdG_Hamiltonian(mu_val_1, w, delta, N, BC, 'chain_1_chain_2');
    [e_vecs_1, e_vals_1, majorana_zero_modes_ref] = ... 
        diagonalise_uncoupled_tetron_via_Kitaev_chains(H_tetron_1(1:2*N,1:2*N), ...
        H_tetron_1((2*N+1):end, (2*N+1):end), 'dirac_zero_modes'); 
    [~, ~, corr_qp_1] = ...
         get_tetron_init_cov_mat_general(e_vecs_1, init_state);     

    mu_val_2 = mu_high + mu_offset*[+1; -1];
    H_tetron_2 = get_tetron_BdG_Hamiltonian(mu_val_2, w, delta, N, BC, 'chain_1_chain_2');
    [e_vecs_2, e_vals_2] = ... 
        diagonalise_uncoupled_tetron_via_Kitaev_chains_no_optim(H_tetron_2(1:2*N,1:2*N), ...
        H_tetron_2((2*N+1):end, (2*N+1):end), 'dirac_zero_modes', majorana_zero_modes_ref); 

    % Put this in the initial basis
    H_tetron_2_init_basis = e_vecs_1'*H_tetron_2*e_vecs_1;     
    H_tetron_2_init_basis = (H_tetron_2_init_basis + H_tetron_2_init_basis')/2;
    [e_vecs_2_init_b, e_vals_2_init_b] = eig(H_tetron_2_init_basis);
    e_vals_2_init_b = diag(e_vals_2_init_b);
    

    % Dodgey handling of the method I used to calculate
    % overlap_inst_basis_states
    mu_vec_binary = ((mu_vec_no_off-mu_low)./mu_high) == 1;

    [Z_exp, X_exp, Y_exp, P_y1y2y3y4, ~, overlap_with_init_state, ...
        overlap_fixed_basis_states, ~, comparison_determinants_mat, ...
        e_vals_mat_curr, ~, quasi_occ_t, ...
        overlap_inst_basis_states] =...
        calc_noisy_tetron_ev_ibs_tls(...
        corr_qp_1, e_vecs_1, e_vals_1, e_vecs_2, e_vals_2, e_vecs_2_init_b, e_vals_2_init_b,...
        mu_vec_binary, delta_t, t_init, t_final, basis_choice, majorana_zero_modes_ref);


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


    %I've changed this so that the initial state has to be "low": uncomment
    %this if you want the option to start of with "high"
    % if mu_vec_no_off(1) == mu_low
    %     mu_init_flag = 0;
    %     mu_val_2 = mu_high + mu_offset*[+1; -1];
    % elseif mu_vec_no_off(1) == mu_high
    %     mu_init_flag = 1;
    %     mu_val_2 = mu_low + mu_offset*[+1; -1];
    % else
    %     error('Invalid mu_vec_no_off(1) value');
    % end