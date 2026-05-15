%---------------------------------------------------------------%
% Get Quasiparticle Operators from BdG Vectors
%---------------------------------------------------------------%
%
% Input: Matrix V with 2N rows and no_vecs columns. Each
% column is represents a vector in the 2N x 2N BdG space. THe vectors do
% not need to be eigenvectors, unlike
% get_QP_operators_from_BdG_Hamiltonian_eigenvectors. AT PRESENT I only
% trust inputs which are BdG vectors with negative energy (corresponding to
% annihilation operators). 
%
% Hamiltonian.
% Output: qp_operators - cell array containing all no_vecs quasiparticle
% operators. BdG vectors with negative energies map to annihilation
% operators and BdG vectors with positive energies map to creation
% operators. If no_vecs == 1, then qp_operators is the N^2 x N^2 matrix for
% the single vectors. 
%      
%
% The basis for d_all is denoted by all N length binary strings (in ascending
% order, eg. 000, 001, 010, 011, 100, 101, 110, 111). This formally
% represents the state:
%   (a_1^dag)^alpha_1 ... (a_N^dag)^alpha_N|vac>
% Where alpha_j is the value of the j't bit in the string. 
% Eg 011 represents alpha_1 = 0, alpha_2 = 1, alpha_3 = 1, and so
%        is the state: a_2^dag a_3^dag |vac>. 
%

function qp_operators = get_QP_operators_from_BdG_vectors(V)
    if mod(size(V,1),2) ~= 0 
        error('V must have 2N rows.');
    end   
    N = size(V,1)./2;
    no_vecs = size(V,2); 
    fermion_operators = get_N_body_fermionic_annihilation_operators(N);
    
    qp_weights = V; 
    qp_operators = cell(size(V,2),1);
    for k = 1:no_vecs
        qp_op_k = zeros(2^N); 
        for m = 1:N
            qp_op_k = qp_op_k + qp_weights(m,k).*fermion_operators{m}' +...
                qp_weights(m+N,k).*fermion_operators{m};
        end
        qp_operators{k} = qp_op_k; %these are ordered in opposite direction to the positive ones.
    end
    if no_vecs == 1
       qp_operators = qp_operators{1}; 
    end
    
end
