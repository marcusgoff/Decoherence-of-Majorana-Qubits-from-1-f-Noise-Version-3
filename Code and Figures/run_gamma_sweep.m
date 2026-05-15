%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: Run Gamma Sweep
%
% Author: Marcus C. Goffage
% Date: 12-Aug-2025
% Affiliation: University of New South Wales
%
% Paper: "Decoherence in Majorana Qubits by 1/f Noise"
% Paper Authors: A. Alase^1, M. C. Goffage^2, M. C. Cassidy^2, 
%                S. N. Coppersmith^{2*}
% Affiliations:  ^1 University of Sydney
%                ^2 University of New South Wales
%                *  Corresponding Author
%
% -------------------------------------------------------------------------
% ABOUT THIS FUNCTION
% -------------------------------------------------------------------------
% Runs a TLS (two level system) simulation for various gamma points.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [L_even_mat, L_odd_mat, X_exp_mat, Z_exp_mat, ...
          L_even_cell, L_even_final_mat, ...
          L_even_final_vec, L_odd_final_vec, gamma_vec, t_vec_noise_cell] = ...
    run_gamma_sweep(mu_low, mu_high, mu_offset, w, delta, N, BC, ...
                       num_gamma_points, gamma_min, gamma_max, ...
                       t_final, num_trials, init_state)

    gamma_vec = logspace(log10(gamma_min), log10(gamma_max), num_gamma_points);
    % --- Diagonalise Hamiltonians ---
    disp('Diagonalising Hamiltonians');
    [e_vecs_1, e_vals_1, e_vecs_2, e_vals_2, ...
     e_vecs_2_init_b, e_vals_2_init_b, ...
     majorana_zero_modes_ref] = ...
        diagonalise_TLS_Hamiltonians(mu_low, mu_high, mu_offset, w, delta, N, BC);
    disp('Diagonalising Complete');

    % --- Prepare final outputs ---
    L_even_cell = cell(length(gamma_vec), 1);
    L_even_final_mat = zeros(length(gamma_vec), num_trials);
    L_even_final_vec = zeros(length(gamma_vec), 1);
    L_odd_final_vec  = zeros(length(gamma_vec), 1);
    t_vec_noise_cell = cell(length(gamma_vec), 1);

    % --- Loop over gamma points ---
    tic
    fprintf('Running Gamma point # of %i:\n', num_gamma_points);
    for gamma_ind = 1:length(gamma_vec)

        fprintf('%i ', gamma_ind);
        if mod(gamma_ind, 10) == 0
            fprintf('\n')
        end

        gamma_tlf = gamma_vec(gamma_ind);
        delta_t_noise = (1 ./ gamma_tlf) ./ 1000;

        % --- Simulation Time ---
        t_vec_noise = (0:delta_t_noise:t_final).';
        t_vec_noise_cell{gamma_ind} = t_vec_noise;
        L_even_mat = zeros(length(t_vec_noise), num_trials);
        L_odd_mat  = zeros(length(t_vec_noise), num_trials);
        X_exp_mat  = zeros(length(t_vec_noise), num_trials);
        Z_exp_mat  = zeros(length(t_vec_noise), num_trials);

        signal_power = (mu_high - mu_low)^2 ./ 4;
        sigma = 1 ./ (2 * gamma_tlf);
        tau   = 1 ./ (2 * gamma_tlf);
        t_init = 0;

        % --- Trial loop ---
        tic
        for trial_no = 1:num_trials
            % Get Noise Trajectory
            [state_vec, t_vec_noise , jump_ind] = ...
                generate_two_parameter_noise_v2(sigma, tau, t_init, t_final, ...
                delta_t_noise, signal_power, 'do_not_remove_dc', 'init_low');

            % Extract noise trajectory only at transitions
            jump_ind_backspace = [jump_ind(2:end); false];
            jump_ind_both = jump_ind_backspace | jump_ind;
            t_vec_jump = [0; t_vec_noise(jump_ind_both); t_vec_noise(end)];
            state_vec_jump = [state_vec(1); state_vec(jump_ind_both); state_vec(end)];

            % Pass Through the Charge Noise Engine
            basis = 'instant_hamil_basis';
            noise_traj_mat = [state_vec_jump(:), state_vec_jump(:)];
            delta_t_vec = diff(t_vec_jump);

            mu_vec_no_off = state_vec_jump + mu_low;
            [X_exp_mean, Z_exp_mean, Y_exp, ...
             P_y1y2y3y4, ~, ~, ...
             ~, ~, ...
             overlap_inst_basis_states] = ...
                charge_noise_engine_ibs_tls_mat_inputs(e_vecs_1, e_vals_1, e_vecs_2, ...
                    e_vals_2, e_vecs_2_init_b, e_vals_2_init_b, ...
                    majorana_zero_modes_ref, ...
                    mu_vec_no_off, mu_offset, mu_low, mu_high, init_state, delta_t_vec, ...
                    t_init, t_final, basis);

            L_odd = 0.5 * (1 - P_y1y2y3y4);
            L_gs  = 1 - overlap_inst_basis_states;
            L_even = L_gs - L_odd;

            % --- Up-Sample ---
            L_even_ups = zeros(length(t_vec_noise), 1); 
            L_odd_ups  = zeros(length(t_vec_noise), 1); 
            X_mean_ups = zeros(length(t_vec_noise), 1); 
            Z_mean_ups = zeros(length(t_vec_noise), 1); 

            for kk = 1:length(t_vec_noise)
                [~, jump_ind_close] = min(abs(t_vec_noise(kk) - t_vec_jump));
                L_even_ups(kk) = L_even(jump_ind_close);
                L_odd_ups(kk)  = L_odd(jump_ind_close);
                X_mean_ups(kk) = X_exp_mean(jump_ind_close);
                Z_mean_ups(kk) = Z_exp_mean(jump_ind_close);
            end

            % Store trial results
            L_even_mat(:, trial_no) = L_even_ups;
            L_odd_mat(:, trial_no)  = L_odd_ups;
            X_exp_mat(:, trial_no)  = X_mean_ups;
            Z_exp_mat(:, trial_no)  = Z_mean_ups;
        end
        toc

        % Store per-gamma results
        L_even_cell{gamma_ind} = L_even_mat;
        L_even_final_mat(gamma_ind, :) = L_even_mat(end,:);
        L_even_final_vec(gamma_ind) = mean(L_even_mat(end,:));
        L_odd_final_vec(gamma_ind)  = mean(L_odd_mat(end,:));
    end
end