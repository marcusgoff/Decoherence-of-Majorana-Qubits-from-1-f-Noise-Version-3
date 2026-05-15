%-----------------------------------------%
% Get Many-Body Excitation Energies from  %
% singe particle e-vals                   %
%-----------------------------------------%
% Given an array of single-body excitation energies (E_excitation_single), this
% function returns all excitation energies in
% the many-body space.
%
% For example, I'm using feeding into this function all positive
% eigenergies of the BdG Hamiltonian and checking this against the full
% Fock space energies.

function E_many = get_Many_Body_Energies(E_excitation_single)

    %We need to compute all possible combinations of E_single
    v = E_excitation_single;
    comb_total = zeros(2.^(length(v)), length(v)); 
    row_count = 1;

    for ii = 1:length(v)
        c_temp = nchoosek(v, ii);
        %Pad with zeros
        c_temp_pad = zeros(size(c_temp,1),length(v));
        c_temp_pad(:, 1:size(c_temp,2)) = c_temp;
        row_count_next = size(c_temp,1) + row_count;
        comb_total(row_count:row_count_next-1,:) = c_temp_pad;
        row_count = row_count_next;
    end
    %This leaves last row as all zeros. This is valid and corresponds to no
    % excitations
    E_many = sort(sum(comb_total, 2));
    
end