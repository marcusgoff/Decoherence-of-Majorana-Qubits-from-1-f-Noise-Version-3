%-------------------------------------------------%
% Tetron Time Evolution - QP Correlation - with Noise - With Various time
% evolution methods. 
%-------------------------------------------------%
% Adapted from calculate_noisy_tetron_evolution_dirac_qp_corr_mat to
% include the input time_evolution_method. 
%
%
% Notes on inputs (not comprehensive): 
% INPUTS:
%          e_vecs - tetron eigenvectors, ordered by chain-1 chain-2 first, 
%                   and by descending eigenvalues second. 
%                   This is achieved by using diagonalise_uncoupled_tetron_via_Kitaev_chains
%                   function. Make sure you use this ordering, or this
%                   function will return an error. 
%          tetron_params = {mu_offset, w, delta, N, BC}
%          mu_vec - [2x num_time_steps] value of mu (noisy) at each time step. 
%                   The first row is the for the top Kitaev chain and the
%                   second row is for the bottom Kitaev Chain.
%          length as t_init:delta_t:t_final.
%          time_evolution_method: Choose method for obtaining the unitary U_t generating
%                          time evolution at each time step, from the noisy
%                          Hamiltonian, H_t, at that time. 
%               'time_ev_diagonalise' - Diagonalise H_t to construct U_t.
%               'time_ev_exp' - Use matlab function "expm" to obtain U_t
%               'time_ev_1st_order' - Use a 1st order Taylor series
%                              expansion to obtain U_t.
%               'time_ev_Magnus' - use the Magnus expansion - TO BE IMPLEMENTED.  
%          basis_choice - Select whether to define the Pauli operators based of the 
%                  Hamiltonian at t = 0 or the instanteneous Hamiltonian
%                  H(t). Allowed values are:
%                  'init_hamil_basis', 'instant_hamil_basis'.
%
% TO DO: CHANGE e_vects_init an e_vals_init to e_vecs_av e_vals_av or
% noiseless. THey are not necessarily the initial e-vecs and e-vals, they are the that
% of the tetron without noise - or the long time averaged tetron.

function [Z_exp, X_exp, Y_exp, P_y1y2y3y4, overlap_with_init_state, ...
        overlap_fixed_basis_states, t_vec, first_site_occ, mu_vec] =...
    calc_noisy_tetron_ev_dirac_qp_corr_mat_general_method_basis(...
    corr_mat_init, e_vecs_init, e_vals_init, tetron_params, mu_vec, ...
    delta_t, t_init, t_final, time_evolution_method, basis_choice)

    if length(tetron_params) ~= 4
        error('tetron_params input must be a length 4 vector');
    end
    w = tetron_params{1};
    delta = tetron_params{2}; 
    N = tetron_params{3}; 
    BC = tetron_params{4};
    
    % You basically initialise with the first eigenstate. 
    % Get H_evolve at each eigentstate. This is just 
    % mu_curr = mu_mean_vec(ii) + mu_offset*[+1, -1];
    % H_evolve = get_tetron_BdG_Hamiltonian(mu_curr, w, delta, N, BC, 'chain_1_chain_2');
 
    %mu_curr = mu_vec(:,1); %2 x 1 array (top and bottom chain chemical potentials)
    %H_tetron = get_tetron_BdG_Hamiltonian(mu_curr, w, delta, N, BC, 'chain_1_chain_2'); 
    
    %U_delta_t = expm(1i*H_evolve.*delta_t);  

    N = length(e_vals_init)./4;
    
    % Run time evolution
    t_vec = t_init:delta_t:t_final;
    Z_exp = zeros(size(t_vec));
    X_exp = zeros(size(t_vec));
    Y_exp = zeros(size(t_vec));
    P_y1y2y3y4 = zeros(size(t_vec));
    overlap_with_init_state = zeros(size(t_vec));
    overlap_fixed_basis_states = zeros(size(t_vec));
    first_site_occ = zeros(size(t_vec));

    omega_total = get_corr_cov_conversion_matrix(N);
    [~, ~, state_0_fixed_corr_qp] = ...
         get_tetron_init_cov_mat_general(e_vecs_init, "0");
    [~, ~, state_1_fixed_corr_qp] = ...
         get_tetron_init_cov_mat_general(e_vecs_init, "1");  
    state_0_fixed_cov_qp = convert_corr_qp_to_cov_qp(state_0_fixed_corr_qp, omega_total);
    state_1_fixed_cov_qp = convert_corr_qp_to_cov_qp(state_1_fixed_corr_qp, omega_total);


    corr_mat_curr = corr_mat_init; % In QP basis - DEFINE INITIAL CONDITION HERE
    cov_mat_qp_init = -1i.*omega_total*(2*corr_mat_curr - eye(4*N))*omega_total'; 
    cov_mat_qp_curr = cov_mat_qp_init;
    %cov_mat_qp_curr = convert_corr_qp_to_cov_qp(corr_mat_curr, omega_total);   
    Z_exp(1) = cov_mat_qp_curr(1,1+N);
    X_exp(1) = cov_mat_qp_curr(1,1+2*N);
    Y_exp(1) = -cov_mat_qp_curr(1,1+3*N); 
    % Leakage:
    mzm_indices_cov_qp = [1,N+1, 1+2*N, 1+3*N];
    P_y1y2y3y4(1) = calc_pfaffian_4x4(cov_mat_qp_curr(mzm_indices_cov_qp, mzm_indices_cov_qp));
    % Overlaps: 
    overlap_with_init_state(1) = calc_overlap_cov_mats(cov_mat_qp_init,cov_mat_qp_curr, 2*N);
    overlap_fixed_basis_states(1) = calc_overlap_cov_mats(cov_mat_qp_curr, state_0_fixed_cov_qp, 2*N) ...
        + calc_overlap_cov_mats(cov_mat_qp_curr, state_1_fixed_cov_qp, 2*N);

    % DEBUG
    % abs(calc_overlap_cov_mats(cov_mat_qp_curr, state_0_fixed_cov_qp, 2*N)).^2
    % abs(calc_overlap_cov_mats(cov_mat_qp_curr, state_1_fixed_cov_qp, 2*N)).^2
    % 
    % [~, ~, state_p_fixed_corr_qp] = ...
    %     get_tetron_init_cov_mat_general(e_vecs_init, "+");  
    % state_p_fixed_cov_qp = convert_corr_qp_to_cov_qp(state_p_fixed_corr_qp, omega_total);
    % [~, ~, state_m_fixed_corr_qp] = ...
    %     get_tetron_init_cov_mat_general(e_vecs_init, "-");  
    % state_m_fixed_cov_qp = convert_corr_qp_to_cov_qp(state_m_fixed_corr_qp, omega_total);
    % 
    % abs(calc_overlap_cov_mats(state_p_fixed_cov_qp, state_0_fixed_cov_qp, 2*N)).^2 
    % abs(calc_overlap_cov_mats(state_m_fixed_cov_qp, state_0_fixed_cov_qp, 2*N)).^2        
    % abs(calc_overlap_cov_mats(state_p_fixed_cov_qp, state_m_fixed_cov_qp, 2*N)).^2   

   %---------------

    if nargout >= 5 % THIS SHOULD BE NARGOUT and >=5
        corr_mat_site_curr = conj(e_vecs_init)*corr_mat_curr*e_vecs_init.';
        first_site_occ(1) = corr_mat_site_curr(1,1);
    end
    

    for ii = 2:length(t_vec) % THIS LOOP IS WHERE YOU SHOULD FOCUS YOUR EFFORTS ON SPEED-UP
        % Get current time evolution unitary
        mu_curr = mu_vec(:,ii);
        H_tetron = get_tetron_BdG_Hamiltonian(mu_curr, w, delta, N, BC, 'chain_1_chain_2'); 
        H_evolve = e_vecs_init'*H_tetron*e_vecs_init; %in initial QP basis   
        
        switch time_evolution_method
            case 'time_ev_exp'
                U_delta_t = expm(1i*H_evolve.*delta_t); 
            case 'time_ev_diagonalise'
                H_e = (H_evolve + H_evolve')/2;
                [V_QP, D_evolve] = eig(H_e);
                D_evolve = real(D_evolve);
                U_delta_t = V_QP*diag(exp(1i*diag(D_evolve).*delta_t))*V_QP';  
            case 'time_ev_diagonalise_debug'
                [V, D_evolve] =...
                    diagonalise_uncoupled_tetron_via_Kitaev_chains(H_tetron(1:2*N,1:2*N), ...
                H_tetron((2*N+1):end, (2*N+1):end), 'dirac_zero_modes');
                D_evolve = diag(D_evolve);
                V_QP = e_vecs_init'*V;
                
                D_evolve = real(D_evolve);
                U_delta_t = V_QP*diag(exp(1i*diag(D_evolve).*delta_t))*V_QP'; 
            case 'time_ev_1st_order'
                H_e = (H_evolve + H_evolve')/2;
                U_delta_t = eye(4*N) + 1i*H_e*delta_t;
            otherwise
                error('Invalid value for time_evolution_method');
        end

        %Time evolution
        corr_mat_next = U_delta_t*corr_mat_curr*U_delta_t';
        corr_mat_curr = corr_mat_next; 

        %[e_vecs, e_vals, majorana_zero_modes] = ... 
        %diagonalise_uncoupled_tetron_via_Kitaev_chains(H_1, H_2, zero_mode_type)

        % Calculate quantities
        if strcmp(basis_choice, 'init_hamil_basis') 
            %Use basis of initial Hamiltonain
            cov_mat_qp_curr = -1i.*omega_total*(2*corr_mat_curr - eye(4*N))*omega_total'; 

        elseif strcmp(basis_choice, 'instant_hamil_basis')
            % Get covariance matrix in basis of instantaneous Hamiltonian
            % eigenvectors
            % I still do time evolution using the initial basis. 

            % Convert correlation matrix in QP basis to site basis
            corr_mat_curr_site = conj(e_vecs_init)*corr_mat_curr*e_vecs_init.'; 
            if ~strcmp(time_evolution_method, 'time_ev_diagonalise')
                [V_QP, D_evolve] = eig(H_e); %In initial Hamiltonian eigenbasis
            end
            % ORDER e_vecs_instant or V_QP properly - decreasing
            % eigenvalues 
            disp('ORDER e_vecs_instant properly');
            % Convert diagonalised H_e back to site basis
            e_vecs_instant = e_vecs_init*V_QP;  
            % This^ has same eigenvalues, D_evolve

            %Debug
            if norm(e_vecs_instant*D_evolve*e_vecs_instant' - H_tetron) >1e-10
                warning('Ive made a mistake');
            end
        
            % Put correlation matrix in instanteneous Hamiltonian basis
            % THIS IS ONLY "CORRECT" IF e_vecs_instant has been ordered
            % properly. 
            corr_mat_curr_inst_qp = e_vecs_instant.'*corr_mat_curr_site*conj(e_vecs_instant);
             % DOUBLE CHECK NEXT LINE IS CORRECT
            cov_mat_qp_curr = convert_corr_qp_to_cov_qp(corr_mat_curr_inst_qp, omega_total);
            % Convert corr_mat_curr_inst_qp to covariance matrix
            % THEN DO CHECKS - make sure a noiseless tetron gets same
            % results for initial and instanteneous bases. 
            error('instant_hamil_basis code not complete');
        else
            error('Invalid Input for basis_choice');
        end


        %cov_mat_qp_curr = convert_corr_qp_to_cov_qp(corr_mat_curr, omega_total);   
        Z_exp(ii) = cov_mat_qp_curr(1,1+N); 
        X_exp(ii) = cov_mat_qp_curr(1,1+2*N);
        Y_exp(ii) = -cov_mat_qp_curr(1,1+3*N);
        P_y1y2y3y4(ii) = calc_pfaffian_4x4(cov_mat_qp_curr(mzm_indices_cov_qp, mzm_indices_cov_qp));     
        overlap_with_init_state(ii) = calc_overlap_cov_mats(cov_mat_qp_init,cov_mat_qp_curr, 2*N);
        overlap_fixed_basis_states(ii) = calc_overlap_cov_mats(cov_mat_qp_curr, state_0_fixed_cov_qp, 2*N) ...
            + calc_overlap_cov_mats(cov_mat_qp_curr, state_1_fixed_cov_qp, 2*N);

        if nargout >= 5
            corr_mat_site_curr = conj(e_vecs_init)*corr_mat_curr*e_vecs_init.';
            first_site_occ(ii) = corr_mat_site_curr(1,1);
        end    
    end
    %corr_mat_final = corr_mat_curr;    
    % warning('Check minus sign in converting corr to cov mat - should be fine.');
    %Dodgely remove imaginary part to Z_exp if its negligble
    
%     if max(abs(imag(Z_exp))) < 1e-10 & ...
%             norm(imag(Z_exp))/norm(real(Z_exp)) < 1e-8
%         Z_exp = real(Z_exp);
%     end
    Z_exp = remove_negligible_imag_parts(Z_exp);
    X_exp = remove_negligible_imag_parts(X_exp);
    Y_exp = remove_negligible_imag_parts(Y_exp);
    P_y1y2y3y4 = remove_negligible_imag_parts(P_y1y2y3y4);
    overlap_with_init_state = remove_negligible_imag_parts(overlap_with_init_state);
    overlap_fixed_basis_states = remove_negligible_imag_parts(overlap_fixed_basis_states);
end

function omega_total = get_corr_cov_conversion_matrix(N)
    omega_s = 1/sqrt(2)*[eye(N),eye(N); 1i*eye(N), -1i*eye(N)]; 
    lam = [fliplr(eye(N)), zeros(N); zeros(N), eye(N)]; 
    omega_total = kron(eye(2),omega_s*lam);
end


function cov_mat_qp = convert_corr_qp_to_cov_qp(corr_mat_qp, omega_total)
    N = size(corr_mat_qp,1)./4;
    
    % %Conversion matrices - inefficient recreating these all the time.
    % % would be better to hold in memory.
    % omega_s = 1/sqrt(2)*[eye(N),eye(N); 1i*eye(N), -1i*eye(N)]; 
    % lam = [fliplr(eye(N)), zeros(N); zeros(N), eye(N)]; 
    % omega_total = kron(eye(2),omega_s*lam);
    
    % Get QP covariance matrix
    cov_mat_qp = -1i.*omega_total*(2*corr_mat_qp - eye(4*N))*omega_total'; 

end


function cov_mat_qp = convert_corr_site_to_cov_qp_deprecated(corr_mat_site, e_vecs)
    N = size(corr_mat_site,1)./4;
    
    %Conversion matrices - inefficient recreating these all the time.
    % would be better to hold in memory.
    omega_s = 1/sqrt(2)*[eye(N),eye(N); 1i*eye(N), -1i*eye(N)]; 
    lam = [fliplr(eye(N)), zeros(N); zeros(N), eye(N)]; 
    omega_total = kron(eye(2),omega_s*lam);
    %Note: lam converts between QPs ordered by descending eigenvalue to the
    %usual way, ie.(dN...d1 d1^dag ... dN^dag)-->(d1...dN d1^dag ... dN^dag)  
    X = [0 1; 1 0];
    S = kron(X, eye(N));
    T = kron(eye(2),S)*e_vecs*kron(eye(2),S);

    % Get QP corr matrix
    corr_mat_qp = T'*corr_mat_site*T;
    
    % Get QP covariance matrix
    cov_mat_qp = -1i.*omega_total*(2*corr_mat_qp - eye(4*N))*omega_total'; 

end


% Local Functions
function cov_mat = convert_corr_to_cov_local_func(corr_site)
   N = size(corr_site,1)./4;
   omega_s = 1/sqrt(2)*[eye(N),eye(N); 1i*eye(N), -1i*eye(N)]; 
   cov_mat = -1i*kron(eye(2),omega_s)*(2*corr_site - eye(4*N))*kron(eye(2),omega_s');
end

function y = remove_negligible_imag_parts(x)

    if max(abs(imag(x))) < 1e-12 % & ...
           % norm(imag(x))/norm(real(x)) < 1e-8
        y = real(x);
    else
        y = x;
    end  
    
end

% NOT SURE ABOUT THIS FUNCTION: THIS IS ONLY VALID FOR WHEN THE TETRONS ARE
% IDENTICAL
function check_inputs(e_vecs, e_vals)

    N = length(e_vals)./4;
    % Check ordering of e-vals as per the doc-string. 
    if norm(e_vals(1:2*N) - e_vals(2*N+1:end)) > 1e-10
        error('Invalid ordering of e_vals, see doc-string');
    end
    e_val_order_check = sort(e_vals(1:2*N), 'descend');
    if norm(e_val_order_check - e_vals(1:2*N)) >1e-10
        error('Invalid ordering of e_vals, see doc-string');
    end
end


%% Calculate 4x4 Pfaffian 
% INPUT: M - must be a 4x4 covariance matrix in majorana basis
% (quasiparticle or site basis is good). 

function pf = calc_pfaffian_4x4(M)
        pf = M(1,2)*M(3,4) - M(1,3)*M(2,4) + M(2,3)*M(1,4); 
end

%% Calculate Gaussian State Overlap from Covariance Matrices
% You may want to put this inline, to speed up the code. 

function overlap = calc_overlap_cov_mats(A,B,N) 
    overlap = 2^(-N)*sqrt(det(A + B)); 
end