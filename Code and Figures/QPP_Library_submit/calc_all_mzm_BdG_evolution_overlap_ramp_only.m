%---------------------------------------------%
% Calc Single BdG Evolution Overlap Ramp Only
%---------------------------------------------%
%
% RAMP ONLY - only used for a ramped chemical potential
%
% This function evolves a single BdG eigenvector in time, as per in the
% input noise_traj, and returns the overlap of that time evolved BdG e-vec
% with each of the instantaneous eigenvectors with positie eigenvalues of
% the top Kitaev Chain in tetron.
%
%
% This function is currently hardcoded to evolve the left MZM of the top
% Kitaev chain. Further extensions would be to add an input which selects a
% BdG e-vec to time evolve. 
%
% Notes on outputs:
%           overlaps_3D_mat - includes overlaps with negative and positive
%                            eigenvalued eigenvectors
%           overlas_exc_total, overlaps_exc_each_mzm - only include the sum 
%                            over eigenvectors with positive eigenvectors. 
%
% Notes on Inputs:
%      noise_traj - this function only accepts a single trajectory, so this
%                   length(t_vec) x 2 matrix.

function [overlaps_exc_total, overlaps_exc_each_mzm, overlaps_3D_mat,...
    mzm_vecs_3D_mat, t_vec, overlaps_exc_modes_all_mzms] = ...
calc_all_mzm_BdG_evolution_overlap_ramp_only(mu_init,...
    mu_offset, w, delta, N, BC, t_init, ramp_rate, delta_t, ramp_height)

    % Get Tetron
    mu_init_full = mu_init + mu_offset.*[+1 -1];
    H_tetron_init =  get_tetron_BdG_Hamiltonian(mu_init_full, w, delta, N, BC, 'chain_1_chain_2');
    
    % Diagonalise Tetron 
    [e_vecs_init, e_vals_init, mzm_vecs_init] = ... 
        diagonalise_uncoupled_tetron_via_Kitaev_chains(H_tetron_init(1:2*N,1:2*N), ...
        H_tetron_init((2*N+1):end, (2*N+1):end), 'dirac_zero_modes');
        
    % % Take left MZM on chain 1 BdG eigenvector
    % init_bdg_state = mzm_vecs_init(:, chosen_ind); 
     
    % Run a for loop, time-evolving the initial e-vector 
    %       At each time-step, take the overlap with each eigsenstate.
    ramp_time = (ramp_height - mu_init)./ramp_rate;
    t_final = ramp_time + t_init;
    t_vec = t_init:delta_t:t_final;
    disp('Check mu_init in calc_all_mzm...');
    
    % Define mu_vec
    mu_vec = (ramp_height - mu_init).*t_vec./ramp_time;
    %mu_vec = mu_vec + mu_init + mu_offset.*[+1; -1];
    
    % psi_mat = zeros(4*N, length(t_vec));
    % psi_mat(:,1) = init_bdg_state;

    mzm_vecs_3D_mat = zeros(4*N, 4, length(t_vec)); 
    mzm_vecs_3D_mat(:,:, 1) = mzm_vecs_init;
        
    overlaps_3D_mat = zeros(4*N, 4, length(t_vec));

    % Repeat for the for other four MZMs. Output each of their outputs with
    % ALL instantaneous BdG eigenvectors. And output the sum of their
    % overlaps with each chain. And output the combined sum. 

    % Collect Overlaps for excited states on first chain
    % overlaps = zeros(length(t_vec), N);
    % 
    % for jj = 1:N
    %     overlaps(1, jj) =  abs(e_vecs_init(:, jj)'*psi_mat(:,1)).^2;
    % end
    
    for mzm_ind = 1:4
        for jj = 1:(4*N)
            overlaps_3D_mat(jj, mzm_ind, 1) = abs(e_vecs_init(:, jj)'*mzm_vecs_3D_mat(:, mzm_ind, 1)).^2;
        end
    end

    % Do Time Evolution
    for ii = 2:length(t_vec)
 
        H_curr =  get_tetron_BdG_Hamiltonian(mu_vec(:, ii), w, delta, N, BC, 'chain_1_chain_2');
       % psi_mat(:, ii) = expm(-1i*H_curr*delta_t)*psi_mat(:, ii-1); % you should also just diagonalise directly
        H_curr = (H_curr + H_curr')./2;
        [V, D] = eig(H_curr);
        D = real(D); 
        U_delta_t = V*diag(exp(-1i*diag(D)*delta_t))*V';
        % psi_mat(:, ii) = U_delta_t*psi_mat(:, ii-1);

        mzm_vecs_3D_mat(:,:, ii) = U_delta_t*mzm_vecs_3D_mat(:,:, ii-1); 

        [e_vecs_curr, e_vals_curr] = ... 
            diagonalise_uncoupled_tetron_via_Kitaev_chains_no_optim(H_curr(1:2*N,1:2*N), ...
            H_curr((2*N+1):end, (2*N+1):end), 'dirac_zero_modes', mzm_vecs_init);
        % 
        % for jj = 1:N
        %     overlaps(ii, jj) =  abs(e_vecs_curr(:, jj)'*psi_mat(:,ii)).^2;
        % end

        for mzm_ind = 1:4
            for jj = 1:(4*N)
                overlaps_3D_mat(jj, mzm_ind, ii) = abs(e_vecs_curr(:, jj)'*mzm_vecs_3D_mat(:, mzm_ind, ii)).^2;
            end
        end
    end

   exc_state_ind = [1:N, (2*N+1):3*N]; %1:4*N;
   exc_state_ind(exc_state_ind == N | exc_state_ind == (N+1) |... 
       exc_state_ind == 3*N | exc_state_ind == (3*N + 1) ) = [];

   overlaps_exc_total = sum(sum(overlaps_3D_mat(exc_state_ind,:,:),1),2);
   overlaps_exc_total = overlaps_exc_total(:);
   
   overlaps_exc_each_mzm_temp = sum(overlaps_3D_mat(exc_state_ind,:,:),1);

   overlaps_exc_each_mzm = zeros(4, length(t_vec));
   for ind = 1:4
        overlaps_exc_each_mzm(ind, :) = overlaps_exc_each_mzm_temp(:, ind, :);
   end

   % Get the BdG leakage for each excited mode, summed up across all 4
   % initial MZM eigenvectors
   overlaps_final_exc = overlaps_3D_mat(exc_state_ind,:,end);
   overlaps_exc_modes_all_mzms = sum(overlaps_final_exc, 2);

end

