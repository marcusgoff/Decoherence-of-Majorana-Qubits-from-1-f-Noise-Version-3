%-----------------------------------%
% Check Anticommutation Rules in BdG
%-----------------------------------%
% Inputs: 
%       vec_1 - length 2N vector in "BdG" Space 
%       vec_2 - length 2N vector in "BdG" Space
%       fermion_type - either "majorana_fermions" or "dirac_fermions"
% Outputs:
%       anticommutator - anticommutator of the operators defined by BdG 
%       vecotrs, vec_1 and vec_2.
%       I.e. anticommutator = {vec_1_hat, vec_2_hat}
%
% Notes: I'm not so proud of this at present. Should have better handling
% for cases of daggers and so forth. At present this only works for the
% eigenvectros directly returned from diagonalising a BdG Hamiltonian.

function anticommutator = calculate_anticommutation_in_BdG_space(vec_1, vec_2)
    % Check valid inputs
%     if ~(strcmp(fermion_type, 'majorana_fermions') || ...
%             strcmp(fermion_type, 'dirac_fermions'))
%         error('check_anticommutation_in_BdG_space');      
%     end
    if ~(size(vec_1,1) ~= 1 || size(vec_1,2) ~= 1) || mod(length(vec_1),2)
        error('A must be a vector of even length');
    end
    if ~(size(vec_2,1) ~= 1 || size(vec_2,2) ~= 1)|| mod(length(vec_2),2)
        error('A must be a vector of even length');
    end
    if length(vec_1) ~= length(vec_2)
        error('A and B must be equal length vectors');
    end
    
    N = length(vec_1)./2;
    anticommutator = 0;
    for m = 1:(N)
        anticommutator = anticommutator +...
                vec_1(m)*vec_2(N+m) + vec_1(N+m)*vec_2(m);
    end
    
end