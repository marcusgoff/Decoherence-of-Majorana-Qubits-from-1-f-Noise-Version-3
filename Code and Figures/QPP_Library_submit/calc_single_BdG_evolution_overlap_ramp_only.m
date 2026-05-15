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
% Notes on Inputs:
%      noise_traj - this function only accepts a single trajectory, so this
%                   length(t_vec) x 2 matrix.

function [overlaps, psi_mat, t_vec] = calc_single_BdG_evolution_overlap_ramp_only(mu_init,...
    mu_offset, w, delta, N, BC, t_init, ramp_rate, delta_t, ramp_height)

    % Get Tetron
    mu_init_full = mu_init + mu_offset.*[+1 -1];
    H_tetron_init =  get_tetron_BdG_Hamiltonian(mu_init_full, w, delta, N, BC, 'chain_1_chain_2');
    
    % Diagonalise Tetron 
    [e_vecs_init, e_vals_init, mzm_vecs_init] = ... 
        diagonalise_uncoupled_tetron_via_Kitaev_chains(H_tetron_init(1:2*N,1:2*N), ...
        H_tetron_init((2*N+1):end, (2*N+1):end), 'dirac_zero_modes');
        
    % Take left MZM on chain 1 BdG eigenvector
    init_bdg_state = mzm_vecs_init(:,1); 
        
    % Run a for loop, time-evolving the initial e-vector 
    %       At each time-step, take the overlap with each eigsenstate.
    ramp_time = (ramp_height - mu_init)./ramp_rate;
    t_final = ramp_time + t_init;
    t_vec = t_init:delta_t:t_final;
    
    % Define mu_vec
    mu_vec = (ramp_height - mu_init).*t_vec./ramp_time;
    mu_vec = mu_vec + mu_init + mu_offset.*[+1; -1];
    
    psi_mat = zeros(4*N, length(t_vec));
    psi_mat(:,1) = init_bdg_state;
    

    % Collect Overlaps for excited states on first chain
    overlaps = zeros(length(t_vec), N);
  
    for jj = 1:N
        overlaps(1, jj) =  abs(e_vecs_init(:, jj)'*psi_mat(:,1)).^2;
    end
    
    % Do Time Evolution
    for ii = 2:length(t_vec)
    
        H_curr =  get_tetron_BdG_Hamiltonian(mu_vec(:, ii), w, delta, N, BC, 'chain_1_chain_2');
       % psi_mat(:, ii) = expm(-1i*H_curr*delta_t)*psi_mat(:, ii-1); % you should also just diagonalise directly
        H_curr = (H_curr + H_curr')./2;
        [V, D] = eig(H_curr);
        D = real(D); 
        U_delta_t = V*diag(exp(-1i*diag(D)*delta_t))*V';
        psi_mat(:, ii) = U_delta_t*psi_mat(:, ii-1);
    
        [e_vecs_curr, e_vals_curr] = ... 
            diagonalise_uncoupled_tetron_via_Kitaev_chains_no_optim(H_curr(1:2*N,1:2*N), ...
            H_curr((2*N+1):end, (2*N+1):end), 'dirac_zero_modes', mzm_vecs_init);
    
        for jj = 1:N
            overlaps(ii, jj) =  abs(e_vecs_curr(:, jj)'*psi_mat(:,ii)).^2;
        end
    end
end