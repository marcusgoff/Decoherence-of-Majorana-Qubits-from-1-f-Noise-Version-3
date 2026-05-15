% -------------------------------------------- %
% L_p (L_odd) Sudden Limit Approximation
% -------------------------------------------- %
%
% April, 2025
% This is the function get_parity_quench_limit with extended functionality
% to:
% Calculate the sudden limit approximation for L_p (L_odd) in the sudden 
% limit, using the updated approximation formula which Abhijeet Alase
% wrote in the supplementary material in March, 2025. 
%
% - - - - - - - - - - - - - - - - - - 
% ORIGINAL FUNCTION DOCUMENTATION:
% - - - - - - - - - - - - - - - - - - 
% This function calculates the parity of the Kiteav-tetron immediately
% following a quench, using the following expression:
%
% <P> = <psi_init | (-i2)^2 gamma_1 gamma_2 gamma_3 gamma_4| psi_init>
% where, |psi_init> is the state immediately before the quench and
%        gamma_i is the i'th MZM immedate after the quench.
%
% According to my notes on overleaf, this is just calculated by taking the
% pfaffian of the covariance matrix for psi_init in the basis of
% Kitaev-tetron after the quench, restricted to the 4x4 matrix corresponding
% % to the operators gamma_i, for i \in [1,4].
%
% This is all done in the instantaneous basis. 

% I should extend this for mu_mean_init and mu_mean_final as 2x1 vectors.

function [P, P_approx, gap_init, gap_final] = get_parity_quench_limit_with_approx(mu_mean_init, mu_mean_final, ...
    mu_offset, w, delta, N, BC)
    %warning('get_parity_quench_limit implented only for same delta_mu on both chains');

    %------------------------------%
    % Get initial basis eigenvectors. Store in U_1.

    % Express |+> state as a covariance matrix in the basis of the initial
    % basis eigenvectors. Store in M_1.

    % Get final basis eigenvectors. store in U_2.

    % Rotate M_1 to basis of new eigenvectors using U_2. Store in M_2.

    % Restrict M_2 to gamma_1 to gamma_4 subspace in the new basis. Use
    % whatever you did in your calculations of the parity before. 

    % Use calc_pfaffian_4x4(M_2_restr) to get the parity of the initial state in
    % the final basis. Store result in P. 
    %-----------------------------%
    mu_init = mu_mean_init(:) + mu_offset*[+1, -1].';
    mu_final = mu_mean_final(:) + mu_offset*[+1, -1].';


    H_tetron_init = get_tetron_BdG_Hamiltonian(mu_init, w, delta, N, BC, 'chain_1_chain_2');
    [e_vecs_init, e_vals_init, mzm_mat_init] = ... 
        diagonalise_uncoupled_tetron_via_Kitaev_chains(H_tetron_init(1:2*N,1:2*N), ...
        H_tetron_init((2*N+1):end, (2*N+1):end), 'dirac_zero_modes');   
    gap_init = sort(abs(e_vals_init)); 
    gap_init = gap_init(5);


    [~, corr_mat_site_init_basis, corr_mat_qp_init_basis] = ...
         get_tetron_init_cov_mat_general(e_vecs_init, "+"); 
    %[cov_site, corr_site, corr_qp] =...
        %   get_tetron_init_cov_mat_general(e_vecs, state)

    H_tetron_final = get_tetron_BdG_Hamiltonian(mu_final, w, delta, N, BC, 'chain_1_chain_2');
    % In the original version of this function (without the approximation)
    % I simply used diagonalise_uncoupled_tetron_via_Kitaev_chains_no_optim
    % in the below line. However I think it's important for the
    % approximation to use the localised MZMs (if I want to use the formula
    % directly). 
   % [e_vecs_final, e_vals_final, mzm_mat_final_v2, comparison_determinants] = ... 
   % diagonalise_uncoupled_tetron_via_Kitaev_chains_no_optim(H_tetron_final(1:2*N,1:2*N), ...
   %    H_tetron_final((2*N+1):end, (2*N+1):end), 'dirac_zero_modes', ...
   %    mzm_mat_init);

    [e_vecs_final, e_vals_final, mzm_mat_final] = ...
        diagonalise_uncoupled_tetron_via_Kitaev_chains(H_tetron_final(1:2*N,1:2*N), ...
            H_tetron_final((2*N+1):end, (2*N+1):end), 'dirac_zero_modes');    

    gap_final = sort(abs(e_vals_final)); 
    gap_final = gap_final(5);

    % Now rotate corr_mat_site_init_basis to cov_mat_qp_final_basis. Same
    % state different basis.

    omega_total = get_corr_cov_conversion_matrix(N);
    corr_mat_qp_final_basis = e_vecs_final.'*corr_mat_site_init_basis*conj(e_vecs_final);
    cov_mat_qp_curr = -1i.*omega_total*(2*corr_mat_qp_final_basis - eye(4*N))*omega_total'; 


    % Restrict to the MZM 4x4 matrix and take the Pfaffian
    mzm_indices_cov_qp = [1,N+1, 1+2*N, 1+3*N];
    P = calc_pfaffian_4x4(cov_mat_qp_curr(mzm_indices_cov_qp, mzm_indices_cov_qp));
    P = remove_negligible_imag_parts(P);


    % Get Approximate Parity Leakage 
    P_approx_a = diag(mzm_mat_final'*mzm_mat_init);
    P_approx_a = prod(remove_negligible_imag_parts(P_approx_a));

   % P_approx_b = diag(mzm_mat_final_v2'*mzm_mat_init);
   % P_approx_b = prod(remove_negligible_imag_parts(P_approx_b));

    % This is the old approximation formula, and should be very close to
    % the new one. 
   % comp_det_prod = remove_negligible_imag_parts(comparison_determinants(1)*comparison_determinants(2));

   % P_approx = [P_approx_a, P_approx_b, comp_det_prod];
    % I'd rather just output P_approx_a or P_approx_b. But it's nice to
    % have all of these as a sanity check. 
    P_approx = P_approx_a;
end



%% Calculate 4x4 Pfaffian 
% INPUT: M - must be a 4x4 covariance matrix in majorana basis
% (quasiparticle or site basis is good). 

function pf = calc_pfaffian_4x4(M)
        pf = M(1,2)*M(3,4) - M(1,3)*M(2,4) + M(2,3)*M(1,4); 
end

%% Get omega_total 

function omega_total = get_corr_cov_conversion_matrix(N)
    omega_s = 1/sqrt(2)*[eye(N),eye(N); 1i*eye(N), -1i*eye(N)]; 
    lam = [fliplr(eye(N)), zeros(N); zeros(N), eye(N)]; 
    omega_total = kron(eye(2),omega_s*lam);
end

%% Remove negligible imag parts

function y = remove_negligible_imag_parts(x)

    if max(abs(imag(x))) < 1e-12 % & ...
           % norm(imag(x))/norm(real(x)) < 1e-8
        y = real(x);
    else
        y = x;
    end  
    
end