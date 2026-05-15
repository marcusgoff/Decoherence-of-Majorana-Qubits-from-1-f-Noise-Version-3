%-----------------------------------------------%
% Get Correlation Matrix From Fock Space        %
%-----------------------------------------------%
% Constructs the correlation matrix from a state vector in 2^N x 2^N Fock
% space. FUTURE: add code to optionally take a density matrix rho, instead
% of state-vector psi. 
% Inputs:
% psi - 2^N length vector describing the state in Fock space
% annihilation_op - size N cell array containing 2^N x 2^N matrices
%                  describing N annihilation operators for all N fermionic
%                  modes satisfying CARs. These can be the fermionic site
%                  operators, quasiparticle operators, or other sets of N
%                  fermionic operators satisfying the CARs. 
%                  OR a size 2N cell array containg 2^N x 2^N matrices
%                  describing 2N majorana operators. 
% input_basis - describes basis of annihilation_op. 
%               either input_basis_dirac_fermions, or
%               input_basis_majorana_fermions.
%               FUNCTIONALITY ONLY WRITTEN FOR input_basis_dirac_fermions
%               so far. 
% output_basis - either output_basis_dirac_fermions, or
%               output_basis_majorana_fermions
%
% input_basis and output_basis should be changed to optional functions. 
% Returns:
% gamma - correlation matrix. Consistent with Surace, 2022, Fermionic
%         Gaussian state: an introduction to numerical approaches. 
%         In dirac fermion basis, this is defined as 
%         <c_vec^dag*c_vec> = <alpha_vec*alpha_vec^dag>
%         with c_vec = (c1,..,cN, c1_dag,...,cN_dag) and alpha_vec =
%         c_vec^dag.
%         In Majorana basis this is defined as: <r_vec*r_vec^dag>, 
%         where r_vec = (x0, ..., xN, p0, ... pN), 
%                 x_j = (c_j + c_j^dag)/sqrt(2)
%                 p_j = (c_j - c_j^dag)/(i*sqrt(2));        
%
% Ideas for additional inputs:
% basis - string either as "dirac_fermions" or "majorana_fermions". Maybe
% call this variable something else like dirac_or_majorana, since there are
% various dirac fermion or majorana bases you could use, which are actually
% specified by annihilation op. 
% could allow state to be a density matrix instead
%
% % qp_op - quasiparticle operators in Fock space (2^N x 2^N matrices)
%         - uh oh these are 2^N x 2^N !!!. Allow for other fermionic
%         operators to define this in other bases. 


function gamma = get_correlation_matrix_from_fock_space(psi, annihilation_op, input_basis, output_basis)
    % Implement option for density matrices as an input
    %if min(size(state)) == 1
    % psi = state % pure state as a vector
    % else 
    % rho = state; % density matrix
    % end
    %psi = state;
    
    % Check valid inputs
    if ~(strcmp(input_basis, 'input_basis_dirac_fermions') || ...
            strcmp(input_basis, 'input_basis_majorana_fermions'))
        error('Invalid input for input_basis in get_correlation_matrix_from_fock_space');      
    end
    if ~(strcmp(output_basis, 'output_basis_dirac_fermions') || ...
            strcmp(output_basis, 'output_basis_majorana_fermions'))
        error('Invalid input for output_basis in get_correlation_matrix_from_fock_space');      
    end   
    % Write code to handle allowed input_basis and output_basis strings.    
    N = log2(length(psi)); 
    omega = 1/sqrt(2)*[eye(N),eye(N); -1i*eye(N), 1i*eye(N)];
    X = [0 1; 1 0];
    
    if strcmp(input_basis,'input_basis_dirac_fermions') %change to strcmp?
        gamma_cr_ann = zeros(N);
        gamma_cr_cr = zeros(N);
        gamma_ann_ann = zeros(N);
        gamma_ann_cr = zeros(N); 

        for i = 1:N
            for j = 1:N
                gamma_cr_ann(i,j) = psi'*annihilation_op{i}'*annihilation_op{j}*psi;
                gamma_cr_cr(i,j) = psi'*annihilation_op{i}'*annihilation_op{j}'*psi;            
                gamma_ann_ann(i,j) = psi'*annihilation_op{i}*annihilation_op{j}*psi;     
                gamma_ann_cr(i,j) = psi'*annihilation_op{i}*annihilation_op{j}'*psi;            
            end
        end

        gamma = [gamma_cr_ann, gamma_cr_cr; gamma_ann_ann, gamma_ann_cr];

        if strcmp(output_basis, "output_basis_majorana_fermions") 
            %dirac gamma to majorana gamma
            U = kron(X, eye(N));
            gamma = omega*U*gamma*U*omega'; 
        end
        
    elseif strcmp(input_basis, 'input_basis_majorana_fermions')
        gamma = zeros(2*N);
        majorana_op = annihilation_op; %input given is Majoranas, renaming to avoid confusion.
        for i = 1:(2*N)
            for j = 1:(2*N)
                gamma(i,j) = psi'*majorana_op{i}*majorana_op{j}'*psi;
            end
        end
        
        if strcmp(output_basis, 'output_basis_dirac_fermions')
            %majorana gamma to dirac gamma
            U = kron(X, eye(N));
            gamma = U*omega'*gamma*omega*U; 
        end
        %Put option to convert back to dirac fermion basis
    else
        error('Invalid input for input_basis in get_correlation_matrix_from_fock_space');
    end
 end



