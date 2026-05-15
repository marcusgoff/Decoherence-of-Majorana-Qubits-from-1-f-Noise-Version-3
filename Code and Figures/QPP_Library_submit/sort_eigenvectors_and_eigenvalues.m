%-----------------------------------%
% Sort Eigenvectors and Eigenvalues
%-----------------------------------%
% Small utility function to sort the eigenvectors and eigenvalues in 
% ascending order of the eigenvalues
%
% Inputs: V and D obtained from [V,D] = eig(A,'vector'), where A is some
% matrix. 
%         order - either 'ascend', 'descend', 'ascend_mag_pos_neg'. 
%                 'ascend_mag_pos_neg' can only be used for BdG vectors. 



function [V, D] = sort_eigenvectors_and_eigenvalues(V, D, order)
    if strcmp(order, 'ascend') || strcmp(order, 'descend')
        [D, ind] = sort(D, order);
        V = V(:,ind);
    elseif strcmp(order, 'ascend_mag_pos_neg') 
        N = length(D)/2;
        if mod(N,2) ~=0
            error('Input order = ascend_mag_pos_neg is only for BdG eigenvectors');
        end
        [D, ind] = sort(D, 'descend');
        V = V(:,ind);
        D(1:N) = fliplr(D(1:N));
        ind_2 = N:-1:1; 
        V(:,1:N) = V(:,ind_2);     
        % I can use this method since the BdG has PHS. I have done it this
        % weird way, to ensure that completely zero eigenvalues also get
        % split up in the expected way for this ordering. 
        %
        % Check energies come out as expected. 
        warning('Function Option Incomplete');
    else
        error('Invalid input for ''order''');
    end
end