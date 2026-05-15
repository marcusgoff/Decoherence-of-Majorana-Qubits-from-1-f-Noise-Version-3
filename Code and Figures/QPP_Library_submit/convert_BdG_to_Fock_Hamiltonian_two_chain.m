%---------------------------------------%
% Convert BdG Hamiltonian to Fock Space
%---------------------------------------%
% Given a 2N x 2N BdG Hamiltonian in the site basis, this function returns
% the corresponding 2^N x 2^N Hamiltonian in Fock Space.
%
% INPUT BdG Hamiltonian must be in the chain-1 chain-2 basis. 

function H = convert_BdG_to_Fock_Hamiltonian_two_chain(H_BdG)
    %CHECK FORM OF BDG HAMILTONIAN

     % Add input checks
     N = size(H_BdG,1)./4;
     H = zeros(2^(2*N), 2^(2*N));
     
     c_ann = get_N_body_fermionic_annihilation_operators(2*N);
     c_ann_1 = c_ann(1:N); 
     c_ann_2 = c_ann(N+1:end);
     
     % I'm going to sweep through every element. This is redundant since H
     % should be Hermitian. 

     for ii = 1:4*N
        for jj = 1:4*N
            if ii <=N
                op_1 = c_ann_1{ii}';
            elseif ii >= N+1 && ii <= 2*N
                op_1 = c_ann_1{ii-N};
            elseif ii >= 2*N+1 && ii <= 3*N
                op_1 = c_ann_2{ii - 2*N}';                                
            elseif ii >= 3*N+1 && ii <= 4*N
                op_1 = c_ann_2{ii-3*N};                
            else
                error('Bug in code');
            end          
            
            if jj <=N
                op_2 = c_ann_1{jj};
            elseif jj>= N+1 && jj <= 2*N
                op_2 = c_ann_1{jj-N}';
            elseif jj>= 2*N+1 && jj <= 3*N
                op_2 = c_ann_2{jj - 2*N};                                
            elseif jj>= 3*N+1 && jj <= 4*N
                op_2 = c_ann_2{jj-3*N}';                
            else
                error('Bug in code');
            end
            
            H = H + 0.5*op_1*H_BdG(ii,jj)*op_2;        
        end
     end
        
     
end
