%------------------------------------------%
% Convert Correlation to Covariance Matrix
%------------------------------------------%
% Takes the correlation matrix (either in dirac fermion or majorana fermion
% basis) and returns the covariance matrix. These matrices are of the form
% in Surace, 2022 - Fermionic Gaussian states: an introduction to numerical
% approaches. 
% INPUTS: 
%        corr_mat: correlation matrix, size 2N x 2N
%        gamma_fermion_or_majorana: valid inputs are
%        "input_basis_dirac_fermions" or "input_basis_majorana_fermions"
% OUTPUTS:
%        cov_mat: covariance matrix, size 2N x 2N as given by eq. 45 in
%        Surace. This in the majorana basis.


function cov_mat = convert_correlation_to_covariance_mat(corr_mat, input_basis)
    % Check valid inputs
    if ~(strcmp(input_basis, 'input_basis_dirac_fermions') || ...
            strcmp(input_basis, 'input_basis_majorana_fermions'))
        error('Invalid input for input_basis in get_correlation_matrix_from_fock_space');      
    end
    if size(corr_mat,1) ~= size(corr_mat,2)
        error('corr_mat must be a square matrix');
    end
    if mod(size(corr_mat,1), 2) ~= 0
        error('corr_mat must be an even dimensioned (2N x 2N) square matrix');
    end
    
    N = size(corr_mat,1)/2;
    % Note omega_s = omega kron(X, eye(N)), corresponds to omega given in
    % Surace. 
    omega_s = 1/sqrt(2)*[eye(N),eye(N); 1i*eye(N), -1i*eye(N)]; 
    %X = [0 1; 1 0];
    %U = kron(X, eye(N));
    
    if strcmp(input_basis,'input_basis_dirac_fermions')
        cov_mat = -1i*omega_s*(2*corr_mat - eye(2*N))*omega_s';
        
    elseif strcmp(input_basis,'input_basis_majorana_fermions')
        cov_mat = -1i*(2*corr_mat - eye(2*N));    
    end

end