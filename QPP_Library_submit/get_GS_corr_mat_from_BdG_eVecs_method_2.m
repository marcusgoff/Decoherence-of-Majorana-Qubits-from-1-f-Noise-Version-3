%----------------------------------------------%
% Get GS Correlation Matrix from BdG Evecs - Method 2
%----------------------------------------------%
% Using "method 2" - I'm just using this as a check.
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

function corr_mat = get_GS_corr_mat_from_BdG_eVecs_method_2(eVecs_BdG, eVals_BdG)
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
    alpha = eVecs_BdG(1:N,(N+1):end);
    beta =  eVecs_BdG((N+1):end,(N+1):end);
    
    % Double check correlation matrix calculation using 
    corr_mat = zeros(2*N, 2*N);
    
    for ii = 1:N
        for jj = 1:N
            for m = 1:N
                corr_mat(ii,jj) = corr_mat(ii,jj) + alpha(ii,m)'*alpha(jj,m);
                corr_mat(ii, jj+N) = corr_mat(ii, jj+N) + alpha(ii,m)'*beta(jj,m);
                corr_mat(ii+N, jj) = corr_mat(ii+N, jj) + beta(ii,m)'*alpha(jj,m);
                corr_mat(ii+N, jj+N) = corr_mat(ii+N, jj+N) + beta(ii,m)'*beta(jj,m);                
             end
        end
    end
end





