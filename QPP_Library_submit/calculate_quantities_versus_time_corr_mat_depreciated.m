%-----------------------------------------------%
% Calculate Quantities Versus Time Correlation Matrix
%-----------------------------------------------%
% This function needs to be cleaned up
%
% Local utility function (accessed by a few scripts) to calculate the time
% evolution of specific quantities using Fock space time evolution. 
% Note that here H_evolve is a 2N x 2N BdG Hamiltonian/
%
% Recommended updates: take in H_evolve in diagonal form and its
% eigenvectors? Is that better than using expm()?
%
% Check definition of time evolution. Your swapping of the C_vec may change
% signs or something?
%
% Surely you could define time evolution in the QP basis instead. 
%
% Inputs: corr_mat_init - must be in dirac fermion basis
% Returns <d_1^dag d_1>

function [overlap, number_operator_exp, t_vec, corr_mat_final] =...
    calculate_quantities_versus_time_corr_mat(corr_mat_init,...
    eVecs_BdG_init, eVals_BdG_init, H_evolve, delta_t, t_init, t_final, qp_index)

    % Go through Surace derivation in Surace. Pretty sure you need to drop
    % the factor of 2 in eq. 110 since you define H_BdG with a factor of
    % 1/2. Go through derivation yourself. 
    % ^I cleared this up and this is correct. Also note that this calculation agreed with
    % the Fock space calculations. 
    disp('Please tidy up calculate_quantities_versus_time_corr_mat');
    % i.e. put checks on allowed inputs. Make the code look nice. 
    U_delta_t = expm(1i*H_evolve.*delta_t);
    [eVecs_BdG_init, eVals_BdG_init] = ...
        sort_eigenvectors_and_eigenvalues(eVecs_BdG_init, eVals_BdG_init, 'descend');
    N = size(corr_mat_init,1)./2;
    
    % Run time evolution
    t_vec = t_init:delta_t:t_final;
    corr_mat_curr = corr_mat_init;
    overlap = zeros(size(t_vec));
    number_operator_exp = zeros(size(t_vec));  
    
    
    %Define QP to site dirac fermions transformation
    X = [0 1; 1 0];
    S = kron(X,eye(N));
    U_QP_site = S*eVecs_BdG_init*S; % DEPRECIATED METHOD
   
    
    corr_mat_init_QP = U_QP_site'*corr_mat_init*U_QP_site; % DEPRECIATED METHOD
    
    number_operator_exp(1) = corr_mat_init_QP(qp_index,qp_index) ;
    cov_mat_init = convert_correlation_to_covariance_mat(corr_mat_init,...
         'input_basis_dirac_fermions');
    overlap(1) = calc_overlap_cov_mats(cov_mat_init, cov_mat_init,N);
        
    for ii = 2:length(t_vec)
        corr_mat_next = U_delta_t*corr_mat_curr*U_delta_t';
        corr_mat_curr = corr_mat_next; 
        
        corr_mat_curr_QP = U_QP_site'*corr_mat_curr*U_QP_site; % DEPRECIATED METHOD
              
        number_operator_exp(ii) = corr_mat_curr_QP(qp_index,qp_index);
        
        cov_mat_curr = convert_correlation_to_covariance_mat(corr_mat_curr,...
            'input_basis_dirac_fermions');
        overlap(ii) = calc_overlap_cov_mats(cov_mat_curr, cov_mat_init,N);        
    end
    corr_mat_final = corr_mat_curr;
    
    %Dodgey handling of negligble imaginary components
    if max(imag(number_operator_exp)) < 1e-14
        number_operator_exp = real(number_operator_exp);
    end
    if max(imag(overlap)) < 1e-14
        overlap = real(overlap);
    end    
    
end


function overlap = calc_overlap_cov_mats(A,B,N)
    overlap = 2^(-N)*sqrt(det(A + B)); 
end