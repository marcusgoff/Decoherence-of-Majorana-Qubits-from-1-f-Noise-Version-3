%-----------------------------------------%
% Get Covariance Matrix from Fock Space
%-----------------------------------------%
% Date: 15.03.24
% Author: Marcus Goffage
%
% Constructs the covariance matrix from a state vector in 2^N x 2^N Fock
% space. FUTURE: add code to optionally take a density matrix rho, instead
% of state-vector psi. 
% Inputs:
% psi - 2^N length vector describing the state in Fock space
% majorana_op -    a size 2*N cell array containg 2^N x 2^N matrices
%                  describing 2*N majorana operators, which satisfy the 
%                  majorana anticommutation relations. 


function cov_mat = get_covariance_matrix_from_fock_space(psi, majorana_op)  
    % cov_mat = i Tr[rho[r,r]]
    N = length(majorana_op)/2; %careful if you ever modify to take in dirac fermions
    cov_mat = zeros(2*N, 2*N);
    for ii = 1:(2*N)
        for jj = 1:(2*N)
            %cov_mat(ii,jj) = 1i*psi'*calc_commutator(majorana_op{ii}, majorana_op{jj})*psi;
            cov_mat(ii,jj) = 1i*psi'*calc_commutator(majorana_op{ii}, majorana_op{jj})*psi;
        end
    end  
end


function result = calc_commutator(A,B)
    result = A*B - B*A;
end

