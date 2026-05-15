%------------------------------------------------------------------------%
% Diagonalise TLS Hamiltonians
% This Function Diagonalises the relevant TLS Hamiltonianians and saves them
% locally
%------------------------------------------------------------------------%



 
function [e_vecs_1, e_vals_1, e_vecs_2, e_vals_2, e_vecs_2_init_b, e_vals_2_init_b, ...
    majorana_zero_modes_ref] = ...
    diagonalise_TLS_Hamiltonians(mu_low, mu_high, mu_offset, w, delta, N, BC)  

  
    mu_val_1 = mu_low +  mu_offset*[+1; -1];
    
    % Get initial Hamiltonian
    H_tetron_1 = get_tetron_BdG_Hamiltonian(mu_val_1, w, delta, N, BC, 'chain_1_chain_2');
    [e_vecs_1, e_vals_1, majorana_zero_modes_ref] = ... 
        diagonalise_uncoupled_tetron_via_Kitaev_chains(H_tetron_1(1:2*N,1:2*N), ...
        H_tetron_1((2*N+1):end, (2*N+1):end), 'dirac_zero_modes');  

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
           

    %% Save Output
    % To do
end