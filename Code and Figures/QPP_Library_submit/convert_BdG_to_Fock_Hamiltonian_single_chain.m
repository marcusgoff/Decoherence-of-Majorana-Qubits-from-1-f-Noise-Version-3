%---------------------------------------%
% Convert BdG Hamiltonian to Fock Space
%---------------------------------------%
% Given a 2N x 2N BdG Hamiltonian in the site basis, this function returns
% the corresponding 2^N x 2^N Hamiltonian in Fock Space.
%
% NOTE: Currently implemented only in the chain-1 chain-2 basis. You should
% code up control parameters so this works in the chain-1 only basis. 

function H = convert_BdG_to_Fock_Hamiltonian_single_chain(H_BdG)

     %CHECK FORM OF BDG HAMILTONIAN

     % Add input checks
     N = size(H_BdG,1)./2;
     H = zeros(2^N, 2^N);
     
     c_ann = get_N_body_fermionic_annihilation_operators(N);
     
     % I'm going to sweep through every element. This is redundant since H
     % should be Hermitian. 

     for ii = 1:2*N
        for jj = 1:2*N
            if ii <=N
                op_1 = c_ann{ii}';
            else
                op_1 = c_ann{ii-N};
            end
            if jj <=N
                op_2 = c_ann{jj};
            else
                op_2 = c_ann{jj-N}';
            end
            
            H = H + 0.5*op_1*H_BdG(ii,jj)*op_2;        
        end
     end
        
     
end
