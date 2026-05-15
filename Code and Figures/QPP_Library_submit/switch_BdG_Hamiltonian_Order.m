%------------------------------------%
% Switch BdG Hamiltonian Order
%------------------------------------%
% This function switches between the two conventions of ordering rows and 
% columns of the BdG Hamiltonian.
%
% Note: So far I have only coded the case of:
% (c1,..cN,c1^dag, ... cN^dag) -> (c1, c1^dag, ..., cN, cN^dag)
% Next I must include a control parameter to encode the opposite direction.


function Hb = switch_BdG_Hamiltonian_Order(Ha)

    % Let's at first do a really lazy way, to make sure I'm sure of what's
    % going on. Should be able to do this in one line of code or with a
    % unitary matrix.
    
    if size(Ha,1) ~= size(Ha,2)
        error('Ha must be a square matrix'); % should also be Hermitian
    end
    N = size(Ha,1)./2;
    Hb = zeros(2*N, 2*N);
    for ib = 1:2*N
        if mod(ib,2) == 1 % odd row
            ia = (ib+1)./2;
        else % odd row
            ia = ib/2 + N;
        end
        for jb = 1:2*N
            if mod(jb,2) == 1 % odd col
                ja = (jb+1)./2;
            else % odd col
                ja = jb/2 + N;
            end             
            Hb(ib,jb) = Ha(ia,ja);
        end
    end
end
