%---------------------------------------------%
% Get N-body Fermionic Annhilation Operators  %
%---------------------------------------------%
% Inputs: Number of fermionic modes
% Outputs: c_all - cell array containing the N annihilation operators,
% represented using the Jordan Wigner transformation. 
%
%
% The basis used corresponds to listing out all binary strings of length N
% where a 0 or 1 in the j'th place denotes whether c_j^dag is acted or not.
% For example:
% 101 --> c_1^dag c_3^dag|vac>
% 0111 --> c_2^dag c_3^dag c_4^dag |vac>
% Due to anticommutation of fermionic operators, the order is important. 

function c_all = get_N_body_fermionic_annihilation_operators(N)

    f = [0 1; 0 0];
    Z = [1 0; 0 -1];
    I = eye(2);
    
    % Store all c_j matrices in a cell array
    c_all = cell(N,1);
    for j = 1:N % create each c_j
       if j == 1
           c_j = f;
       else 
           c_j = Z;
       end
       for ii = 2:N
           if ii < j
               on_site_mat = Z;
           elseif ii == j
               on_site_mat = f;
           elseif ii > j
               on_site_mat = I;
           end
           c_j = kron(c_j, on_site_mat);
       end     
       c_all{j,1} = c_j;
    end

end