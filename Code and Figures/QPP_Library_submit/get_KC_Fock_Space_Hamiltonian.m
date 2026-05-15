%---------------------------------------------------%
% Get KC Full Fock Space Hamiltonian - Generalised (any length Chain)
%---------------------------------------------------%
% This function returns the BdG Hamiltonian for the 
% Kitaev chain.
% Coded fo
%
% Inputs: mu, w, delta - KC parameters:
%        PB - Boundary conditions. Valid inputs are 
%             "OBC" (open b.c.) and "PBC" (periodic b.c.).
%
% Consistency checks performed: This gives identical output as
% get_KC_Fock_space_Hamiltonian_Hard_Coded for N = 2 and N = 3. No other
% cases have been checked. 
%
function H = get_KC_Fock_Space_Hamiltonian(mu, w, delta, N, BC)

    c_all = get_N_body_fermionic_annihilation_operators(N);
    
    % Get on_site term
    on_site = zeros(2^N);
    for j = 1:N
       on_site = on_site - mu*c_all{j}'*c_all{j}; 
    end
    
    % Get hopping left to right term
    hopping = zeros(2^N);
    for j = 1:(N-1)
        hopping = hopping + -w*c_all{j}'*c_all{j + 1} + ...
            -w*c_all{j+1}'*c_all{j};
    end
    
    % Get pairing term
    pairing = zeros(2^N);
    for j = 1:(N-1)
       pairing = pairing + delta.*c_all{j}*c_all{j+1} + ...
           conj(delta).*c_all{j+1}'*c_all{j}';
    end
    H = on_site + hopping + pairing;
    
    if BC == 'PBC' 
        boundary_hopping = -w*c_all{end}'*c_all{1} + ...
            -w*c_all{1}'*c_all{end}; 
        boundary_pairing = delta.*c_all{end}*c_all{1} + ...
           conj(delta).*c_all{1}'*c_all{end}';
       H = H + boundary_hopping + boundary_pairing; 
    end
    
    %Check Hermitian
    if norm(H - H') > 1e-10
        warning('Error: Fock space Hamiltonian not Hermitian');
    end
end

%function

