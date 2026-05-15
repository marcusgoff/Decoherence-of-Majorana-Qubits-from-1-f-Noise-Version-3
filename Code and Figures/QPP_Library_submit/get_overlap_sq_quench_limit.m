%--------------------------------%
% Get Overlap Sq Quench Limit
%--------------------------------%
% This function calculates the overlap-mod-squaed of the Kiteav-tetron immediately
% following a quench, using the following expression:
%
% overlap_sq = |<0_f|psi_init>|.^2 + |<1_f|psi_init>|.^2.
%
% To do this I use covariance marices to represent all the states and use
% the formula for calculating overlaps (in Bravyi I believe, and used
% multiple places throughout my code). 
%
% I use this to calculate L_c0 = 1 - overlap_sq. 
%
% This is all done in the instantaneous basis. 

% I should extend this for mu_mean_init and mu_mean_final as 2x1 vectors.

function [overlap_sq, gap_init, gap_final] = get_overlap_sq_quench_limit(mu_mean_init, mu_mean_final, ...
    mu_offset, w, delta, N, BC)
    %warning('get_parity_quench_limit implented only for same delta_mu on both chains');

    %------------------------------%
    % From Parity version of this code
    % Get initial basis eigenvectors. Store in U_1.

    % Express |+> state as a covariance matrix in the basis of the initial
    % basis eigenvectors. Store in M_1.

    % Get final basis eigenvectors. store in U_2.
    
    % Get covariance matrices for |0_f> and |1_f> in final basis (this is
    % just their generic form!). 

    % Rotate M_1 to basis of new eigenvectors using U_2. Store in M_2.

    % Use overlap formula to calculate |<0_f|psi_init>|.^2 and
    % |<1_f|psi_init>|.^2.

    % Return overlap_q = |<0_f|psi_init>|.^2 + |<1_f|psi_init>|.^2.
    %-----------------------------%
    overlap_sq = 0; % STUB

    mu_init = mu_mean_init(:) + mu_offset*[+1, -1].';
    mu_final = mu_mean_final(:) + mu_offset*[+1, -1].';


    H_tetron_init = get_tetron_BdG_Hamiltonian(mu_init, w, delta, N, BC, 'chain_1_chain_2');
    [e_vecs_init, e_vals_init, majorana_zero_modes_ref_init] = ... 
        diagonalise_uncoupled_tetron_via_Kitaev_chains(H_tetron_init(1:2*N,1:2*N), ...
        H_tetron_init((2*N+1):end, (2*N+1):end), 'dirac_zero_modes');   
    gap_init = sort(abs(e_vals_init)); 
    gap_init = gap_init(5);


    [~, corr_mat_site_init_basis, corr_mat_qp_init_basis] = ...
         get_tetron_init_cov_mat_general(e_vecs_init, "+"); 


    H_tetron_final = get_tetron_BdG_Hamiltonian(mu_final, w, delta, N, BC, 'chain_1_chain_2');
    % [e_vecs_final, e_vals_final, majorana_zero_modes_ref_final] = ... 
    %     diagonalise_uncoupled_tetron_via_Kitaev_chains(H_tetron_final(1:2*N,1:2*N), ...
    %     H_tetron_final((2*N+1):end, (2*N+1):end), 'dirac_zero_modes'); 
    %NOTE: you may want to replace above line with the delocalised MZM guy, but
    %using the reference majoranas. That might be important tbh. 
    [e_vecs_final, e_vals_final, majorana_zero_modes_ref_final, comparison_determinants] = ... 
    diagonalise_uncoupled_tetron_via_Kitaev_chains_no_optim(H_tetron_final(1:2*N,1:2*N), ...
        H_tetron_final((2*N+1):end, (2*N+1):end), 'dirac_zero_modes', ...
        majorana_zero_modes_ref_init);
    gap_final = sort(abs(e_vals_final)); 
    gap_final = gap_final(5);

    % Get Covariance matrices for |0_f> and |1_f> 
    omega_total = get_corr_cov_conversion_matrix(N);
    [~, ~, state_0_final_corr_qp] = ...
         get_tetron_init_cov_mat_general(e_vecs_init, "0"); 
    [~, ~, state_1_final_corr_qp] = ...
         get_tetron_init_cov_mat_general(e_vecs_init, "1"); 
    state_0_final_cov_qp = convert_corr_qp_to_cov_qp(state_0_final_corr_qp, omega_total);
    state_1_final_cov_qp = convert_corr_qp_to_cov_qp(state_1_final_corr_qp, omega_total);

    % Now rotate corr_mat_site_init_basis to cov_mat_qp_final_basis. Same
    % state different basis.
    corr_mat_qp_final_basis = e_vecs_final.'*corr_mat_site_init_basis*conj(e_vecs_final);
    cov_mat_qp_curr = -1i.*omega_total*(2*corr_mat_qp_final_basis - eye(4*N))*omega_total'; % Could use the function for this.


    % Now Calculate the Mod Overlap Squared
    overlap_sq = calc_overlap_cov_mats(state_0_final_cov_qp, cov_mat_qp_curr,2*N)  + ...
        calc_overlap_cov_mats(state_1_final_cov_qp, cov_mat_qp_curr,2*N) ;
   
     overlap_sq = remove_negligible_imag_parts(overlap_sq);

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


%% Calculate Gaussian State Overlap from Covariance Matrices
% You may want to put this inline, to speed up the code. 

function overlap = calc_overlap_cov_mats(A,B,N) 
    overlap = 2^(-N)*sqrt(det(A + B)); 
end

%% Convert corr_qp to cov_qp
function cov_mat_qp = convert_corr_qp_to_cov_qp(corr_mat_qp, omega_total)
    N = size(corr_mat_qp,1)./4;
    
    % %Conversion matrices - inefficient recreating these all the time.
    % % would be better to hold in memory.
    % omega_s = 1/sqrt(2)*[eye(N),eye(N); 1i*eye(N), -1i*eye(N)]; 
    % lam = [fliplr(eye(N)), zeros(N); zeros(N), eye(N)]; 
    % omega_total = kron(eye(2),omega_s*lam);
    
    % Get QP covariance matrix
    cov_mat_qp = -1i.*omega_total*(2*corr_mat_qp - eye(4*N))*omega_total'; 

end