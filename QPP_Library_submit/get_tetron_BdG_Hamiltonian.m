%-------------------------------%
% Get Tetron BdG Hamiltonian
%-------------------------------%
% Author: Marcus
% Date: 28.03.24
%
% This function returns the BdG Hamiltonian for the tetron qubit, modelled
% as two uncoupled Kitaev chains.
%
% Inputs: mu, w, delta, N - KC parameters:
%        BC - Boundary conditions. Valid inputs are 
%             "OBC" (open b.c.) and "PBC" (periodic b.c.).
%        Note: N is the number of sites in one of the Kitaev chains, so in
%        total there are 2N sites and the BdG Hamiltonian is 4N x 4N. 
%        basis_str - 'chain_1_chain_2' or 'c_c_dag"
%        mu - to specify different values for top/bottom chain use a 2x1
%        vector. 
%
% Future additions:
%        Modify code to allow for mu, w, and delta to be length 2 vectors
%        where the two values give the parameters for the each of the two
%        uncoupled Kitaev chains. 
%        Let basis_str be an optional input, and by default do
%        "chain_1_chain_2"

function H = get_tetron_BdG_Hamiltonian(mu, w, delta, N, BC, basis_str)
    if length(mu) == 2
        mu_1 = mu(1);
        mu_2 = mu(2); 
    elseif length(mu) == 1;
        mu_1 = mu; mu_2 = mu;
    else
        error('mu must be 1x1 or 2x1 vector');
    end
    
%     A_1 = -mu*diag(ones(N,1)) - w*diag(ones(N-1,1),1) + ...
%         -w*diag(ones(N-1,1),-1);
%     B_1 = delta*diag(ones(N-1,1),1) + -delta*diag(ones(N-1,1),-1);
%     
%     if BC == "PBC"
%         A_1(1,end) = - w; A_1(end,1) = -w;
%         B_1(1,end) = -delta; B_1(end,1) = delta;
%     end
% 
%     A_2 = A_1; %Change when you make future alteration
%     B_2 = B_1; 

    [A_1, B_1] = get_A_B_matrices(mu_1, w, delta, N, BC);
    [A_2, B_2] = get_A_B_matrices(mu_2, w, delta, N, BC);
    
    if strcmp(basis_str,'c_c_dag')
        A = [A_1, zeros(N); zeros(N), A_2];
        B = [B_1, zeros(N); zeros(N), B_2]; 
        H = [A, -conj(B); B, -conj(A)];
    elseif strcmp(basis_str, 'chain_1_chain_2')
        H_KC_1 = [A_1, -conj(B_1); B_1, -conj(A_1)];
        H_KC_2 = [A_2, -conj(B_2); B_2, -conj(A_2)];
        H = [H_KC_1, zeros(2*N); zeros(2*N), H_KC_2];
    else 
        error('Invalid input for: basis');
    end
        
end

function [A, B] = get_A_B_matrices(mu, w, delta, N, BC)

    A = -mu*diag(ones(N,1)) - w*diag(ones(N-1,1),1) + ...
        -w*diag(ones(N-1,1),-1);
    B = delta*diag(ones(N-1,1),1) + -delta*diag(ones(N-1,1),-1);
    
    if BC == "PBC"
        A(1,end) = - w; A(end,1) = -w;
        B(1,end) = -delta; B(end,1) = delta;
    end
    
end