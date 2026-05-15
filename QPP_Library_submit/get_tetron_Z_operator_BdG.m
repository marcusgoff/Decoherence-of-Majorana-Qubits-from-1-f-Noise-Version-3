%--------------------------%
% Get Tetron X Operator BdG
%--------------------------%
% I'm susipicious that this could be off by a factor of 2. Double check.
% Same with X_BdG. 
%
% Returns the tetron Z operator as BdG matrix, defined as
% <-2i gamma_1 gamma_2 > in the
% chain_1-chain_2 site-local dirac fermion basis (see note below). Note
% gamma_1 and gamma_2 are defined by the quasiparticles, they are not in
% generally entirely isolated to the ends of the chain. 
% INPUTS:
%          e_vecs - tetron eigenvectors, ordered by chain-1 chain-2 first, 
%                   and by descending eigenvalues second. 
%                   This is achieved by using diagonalise_uncoupled_tetron_via_Kitaev_chains
%                   function. Make sure you use this ordering, or this
%                   function will return an error. 
%          e_vals - tetron eigenvalues, ordered by descending eigenvalues  
%
% OUTPUTS:
%         X_BdG - tetron X operator in dirac fermion site
%         basis with chain-1-chain-2 ordering
%         X_BdG_QP - tetron X operator in the dirac QP basis (ordered 
%                       by descending eigenvalue!
% To do: add check for ordering of e_vals and e_vecs. 


function [Z_BdG, Z_BdG_QP] = get_tetron_Z_operator_BdG(e_vecs, e_vals)

    check_inputs(e_vecs, e_vals);
    N = length(e_vals)./4;

    A = zeros(2*N); 
    A(N:N+1,N:N+1) = [-1i/2, 0; 0, 1i/2]; 
    % Define X_BdG in quasiparticle dirac basis (without -2i factor):
    % ordered  by eigenvalues
    gamma_1_gamma_2_BdG_QP = [A, zeros(2*N); zeros(2*N), zeros(2*N)];
    % Rotate to site basis and multiply by 2i:
    Z_BdG =  - 2i*e_vecs*gamma_1_gamma_2_BdG_QP*e_vecs'; 
    Z_BdG_QP = -2i*gamma_1_gamma_2_BdG_QP;
     
end
% sort_eigenvalues_and_eigenvectors

% Local Functions
function check_inputs(e_vecs, e_vals)

    if size(e_vecs,1) ~= size(e_vecs,2)
        error('Input e_vecs must be a 4N x 4N matrix');
    end
    if size(e_vals,1) ~= 1 && size(e_vals,2) ~= 1
       error('e_vals must be a 4N x 1 vector'); 
    end
    if size(e_vecs,1) ~= length(e_vals)
        error('e_vecs must be 4N x 4N and e_vals must be 4N x 1')
    end
    if mod(length(e_vals),4)
        error('e_vals must be a 4N x 1 vector');
    end     
    N = length(e_vals)./4;
    
%     % Check ordering of e-vals as per the doc-string. - only valid if the
%     two Kitaev chains are identical. 
%     if norm(e_vals(1:2*N) - e_vals(2*N+1:end)) > 1e-10
%         error('Invalid ordering of e_vals, see doc-string');
%     end
%     e_val_order_check = sort(e_vals(1:2*N), 'descend');
%     if norm(e_val_order_check - e_vals(1:2*N)) >1e-10
%         error('Invalid ordering of e_vals, see doc-string');
%     end
end
    

    
    
    
    