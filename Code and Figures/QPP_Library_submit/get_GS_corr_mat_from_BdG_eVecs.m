%----------------------------------------------%
% Get GS Correlation Matrix from BdG Evecs
%----------------------------------------------%
% This function returns the ground state correlation matrix given
% eigenvectors and eigenvalues of a BdG Hamiltonian. Be careful if there
% are degenerate ground states, in which case this code may need to be
% modified. 
% INPUTS:
%       eVecs_BdG: 2N x 2N matrix where the j'th column is the j'th
%       eigenvector of a BdG Hamiltonian.
%       eVals_BdG: length 2N vector where the j'th entry is the j'th
%       eigenvalue of the BdG Hamiltonian. 
% OUTPUTS:
%       corr_mat: 2N x 2N correlation matrix corresponding to the ground
%       state of the BdG Hamiltonian.
% Note: the eigenvalues do not need to be ordered, they only need to
% corresponding to the appropriate columns of eVecs_BdG.

function [corr_mat, corr_GS_QP_basis] = get_GS_corr_mat_from_BdG_eVecs(eVecs_BdG, eVals_BdG)
    % Check inputs:
    if (size(eVecs_BdG,1) ~= size(eVecs_BdG,2)) || mod(size(eVecs_BdG,1),2)
        error('eVecs_BdG must be a 2N x 2N matrix');
    end
    N = size(eVecs_BdG,1)./2;
    if (size(eVals_BdG,1) ~= 1 && size(eVals_BdG,2) ~= 1) || ...
            length(eVals_BdG) ~= 2*N
        error('eVals_BdG must be a length 2N vector');
    end

    [eVecs_BdG, e_vals_BdG] = sort_eigenvectors_and_eigenvalues(eVecs_BdG, eVals_BdG, 'descend');

    % Here I follow the method of Surace (2022) arXiv:2111.08343v2,
    % specically see page 15 in equations 62 and 63. I need to do an extra
    % transformation since I have ordered my "C" or "alpha" vector with
    % annihilation operators followed by creation operators (Surace does
    % the opposite way around). 
    X = [0 1; 1 0];
    S = kron(X,eye(N));
    U = eVecs_BdG;
    corr_GS_QP_basis = [zeros(N), zeros(N); zeros(N), eye(N)]; 
    corr_BdG_dirac_deprec = S*U*S*corr_GS_QP_basis*S*U'*S; % dirac fermion basis - % you can delete this
       
    corr_BdG_dirac_new_method = conj(U)*corr_GS_QP_basis*U.';
    
    %These methods should give the same result ONLY because there is a
    %degenerate ground state and the mistake implicit in the depreciated
    %method doesn't come into effect since most of the elements are zero.

    deprec_check = norm(corr_BdG_dirac_new_method - corr_BdG_dirac_deprec); % you can delete this.
    
    if deprec_check > 1e-10
        warning('SERIOUS WARNING: I have a misunderstanding!');
    end
    
    corr_mat = corr_BdG_dirac_new_method;
    %This should be equivalent: corr_BdG_dirac = conj(U)*corr_GS_QP_basis*U.'
    
    % Double check correlation matrix calculation using 
    
end





