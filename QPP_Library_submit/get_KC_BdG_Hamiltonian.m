%-----------------------------------%
% Get Kitaev Chain Hamiltonian      %
%-----------------------------------%
% This function returns the BdG Hamiltonian for the 
% Kitaev chain.
%
% Inputs: mu, w, delta - KC parameters:
%        BC - Boundary conditions. Valid inputs are 
%             "OBC" (open b.c.) and "PBC" (periodic b.c.).
%
% Future alteration: change 

function [H, A, B] = get_KC_BdG_Hamiltonian(mu, w, delta, N, BC)
    if ~isreal(mu)
       error('Input parameter mu must be real.') 
    end
        
    A = -mu*diag(ones(N,1)) - w*diag(ones(N-1,1),1) + ...
        -w*diag(ones(N-1,1),-1);
    B = delta*diag(ones(N-1,1),1) + -delta*diag(ones(N-1,1),-1);

    if BC == "PBC"
        A(1,end) = - w; A(end,1) = -w;
        B(1,end) = -delta; B(end,1) = delta;
    end
    
    H = [A, -conj(B); B, -conj(A)];

end