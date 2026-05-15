%---------------------------------------------------------------%
% Get Quasiparticle Operators from BdG Hamiltonian Eigenvectors
%---------------------------------------------------------------%
%
% Input: Matrix V and D obtained from [V, D] = eig(H, 'diag')
%        where H is  2N x 2N BdG
% Hamiltonian.
% Output: qp_operators - cell array containing all N quasiparticle annhilation
% operators.
%         qp_energies - vector containing the energies of all quasiparticles.
%
% The basis for d_all is denoted by all N length binary strings (in ascending
% order, eg. 000, 001, 010, 011, 100, 101, 110, 111). This formally
% represents the state:
%   (a_1^dag)^alpha_1 ... (a_N^dag)^alpha_N|vac>
% Where alpha_j is the value of the j't bit in the string. 
% Eg 011 represents alpha_1 = 0, alpha_2 = 1, alpha_3 = 1, and so
%        is the state: a_2^dag a_3^dag |vac>. 
%
%
% For debugging you may want to extract more variables out:
%[qp_op_from_BdG, qp_e_pos, qp_op_from_BdG_neg, qp_e_neg] = get_...
%

function [qp_operators, qp_energies] = get_QP_operators_from_BdG_Hamiltonian_eigenvectors(V, D)
    N = length(D)./2;
    
    % Ensure eigenvalues and eigenvectors are sorted:
    [V, D] = sort_eigenvectors_and_eigenvalues(V,D, 'descend');
    fermion_operators = get_N_body_fermionic_annihilation_operators(N);

%   % OPTION 1 - use positive energies (and implicit PH transformation)
%     qp_weights = conj(V(:,D>0)); %each column corresponds to a different 
%                                   % quasiparticle.
%     qp_energies = D(D>0);
%     
%     qp_operators = cell(N,1);
%     for k = 1:N
%         qp_op_k = zeros(2^N); 
%         for m = 1:N
%             qp_op_k = qp_op_k + qp_weights(m,k).*fermion_operators{m} +...
%                 qp_weights(m+N,k).*fermion_operators{m}';
%         end
%         qp_operators{k} = qp_op_k; 
%     end
    
    % OPTION 2 - use negative energies:
    qp_weights_neg = V(:,D<0); % each column corresponds to a different 
                                  % quasiparticle.
    qp_energies_neg = -D(D<0);
    
    qp_operators_neg = cell(N,1);
    for k = 1:N
        qp_op_k = zeros(2^N); 
        for m = 1:N
            qp_op_k = qp_op_k + qp_weights_neg(m,k).*fermion_operators{m}' +...
                qp_weights_neg(m+N,k).*fermion_operators{m};
        end
        qp_operators_neg{k} = qp_op_k; %these are ordered in opposite direction to the positive ones.
    end
    qp_operators = qp_operators_neg;
    qp_energies = qp_energies_neg; 
      
    %NOTE:
    % The two options may give different qp_operators by a factor of a
    % global phase. This is simply because when Matlab diagonalises H each
    % eigenvector can have any arbitrary global phase. This isn't
    % particularly an issue. 
end



