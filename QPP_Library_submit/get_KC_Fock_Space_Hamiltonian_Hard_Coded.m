%---------------------------------------------------%
% Get KC Full Fock Space Hamiltonian - Hard Coded
%---------------------------------------------------%
% This function returns the BdG Hamiltonian for the 
% Kitaev chain.
%
% Inputs: mu, w, delta - KC parameters:
%        PB - Boundary conditions. Valid inputs are 
%             "OBC" (open b.c.) and "PBC" (periodic b.c.).
%

function H = get_KC_Fock_Space_Hamiltonian_Hard_Coded(mu, w, delta, N, BC)

    f = [0 1; 0 0];
    Z = [1 0; 0 -1];
    I = eye(2);
    
    if N == 2
        c1 = kron(f, I);
        c2 = kron(Z, f);
        
        H = -mu*(c1'*c1 + c2'*c2) - w*c1'*c2 - w*c2'*c1 +...
            delta*c1*c2 + conj(delta)*c2'*c1';       
    elseif N == 3
        c1 = kron(kron(f, eye(2)), eye(2));
        c2 = kron(kron(Z, f), eye(2));
        c3 = kron(kron(Z, Z), f);
        
        on_site = -mu*(c1'*c1 + c2'*c2 + c3'*c3);
        hopping_LR =  c1'*c2 + c2'*c3;
        hopping_RL =  c2'*c1 + c3'*c2;
        pairing_LR = c1*c2 + c2*c3;
        pairing_RL = c2'*c1' + c3'*c2';
        H = on_site + -w*hopping_LR + - w*hopping_RL +...
            delta*pairing_LR + conj(delta).*pairing_RL;
    else
        warning('Solution for this N value not hardcoded'); 
    end
    
    if BC == "PBC"
        warning("PBC not yet implemented");
    end

    %Check Hermitian
    if norm(H - H') > 1e-10
        warning('Error: Fock space Hamiltonian not Hermitian');
    end
end

%function

