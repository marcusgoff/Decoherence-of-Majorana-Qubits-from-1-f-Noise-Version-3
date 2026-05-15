%-------------------------------------------------%
% Tetron Time Evolution - QP Correlation - with Noise - With Various time
% evolution methods. 
%
% DO NOT USE THIS FOR ACTUAL SIMULATIONS. ONLY USED FOR CHECKS IN SETUP.
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
%           delta_t - either a scalar or vector of delta_t values
%
% OUTPUTS: 
%           P_y1y2y3y4_v1 - using localised majoranas
%           P_y1y2y3y4_v2 - using delocalised majoranas
%            ...
%           e_vals_mat - the (i,j)'th element is the j'th instantaneous 
%                             eigenvalue at time step t_i. This matrix is
%                             only calculated if basis_choice =
%                             'instant_hamil_basis'. Otherwise it is a
%                             matrix of NaN's. 
%           quasi_occ_final - N x 2 matrix of the occupations of the Quasiparticles
%                         at t = t_final. Column 1 is for KC 1 and Column 2
%                         is for KC 2. The ordering for the rows is the same as stored
%                         in e_vecs, i.e ordered by DECREASING eigenvalue.
%                         i.e [d_n, ..., d_1] ordering. 
%           quasi_occ_t - 2*N x length(t_vec) containing quasiparticles as
%           a function of time. quasi_occ_t(:,end) = quasi_occ_final(:).
%           Don't output this massive matrix unless you need to!
%
%           Note that the "fixed" basis states are the computational basis % 
%           states defined using the eigenvectors of H at t=0.
%           This is as opposed to the instantaneous basis states. 
%
% TO DO: CHANGE e_vects_init an e_vals_init to e_vecs_av e_vals_av or
% noiseless. THey are not necessarily the initial e-vecs and e-vals, they are the that
% of the tetron without noise - or the long time averaged tetron.
%
% BIG CHANCE TO SPEED UP: Looking over this code, I'm pretty sure in the
% case where basis_choice = 'instant_hamil_basis' and  time_evolution_method =
% 'time_ev_diagonalise', I'm uncessarily diagonalising the Hamiltonian
% twice. I should stop that from happening to be honest. I could even write
% a streamlined version of this code where 'time_ev_diagonalise' is the
% only method used, since in practice I'm not using any of the others. 
%
% Note that "overlaps" outputted are actually the abs(overlap).^2.
%
% TO FIX!!!!:
%       %overlap_inst_basis_states - you need to code up calculating the
%       instantaneous "0" and "1" states. At the moment this variable just
%       outputs the same values as overlap_fixed_basis_states.



function [Z_exp, X_exp, Y_exp, P_y1y2y3y4_v1, P_y1y2y3y4_v2, overlap_with_init_state, ...
        overlap_fixed_basis_states, t_vec, first_site_occ, mu_vec,...
        comparison_determinants_mat, ...
        e_vals_mat, quasi_occ_final, quasi_occ_t, overlap_inst_basis_states] =...
    calc_noisy_tetron_ev_inst_basis_setup(...
    corr_mat_init, e_vecs_init, e_vals_init, tetron_params, mu_vec, ...
    delta_t, t_init, t_final, time_evolution_method, basis_choice, majorana_zero_modes_ref, mu)


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
    if isscalar(delta_t)
        t_vec = t_init:delta_t:t_final;
    else
        t_vec = t_init + [t_init, cumsum(delta_t(:)).'];
        if t_vec(end) ~= t_final
            warning('delta_t vector does not match t_init and t_final');
            t_final = t_vec(end);
        end
    end
    Z_exp = zeros(size(t_vec));
    X_exp = zeros(size(t_vec));
    Y_exp = zeros(size(t_vec));
    P_y1y2y3y4_v1 = zeros(size(t_vec));
    overlap_with_init_state = zeros(size(t_vec));
    overlap_fixed_basis_states = zeros(size(t_vec));
    first_site_occ = zeros(size(t_vec));
    quasi_occ_final = zeros(N, 2);
    quasi_occ_t = zeros(2*N, length(t_vec));
    
    omega_total = get_corr_cov_conversion_matrix(N);
    [~, ~, state_0_fixed_corr_qp] = ...
         get_tetron_init_cov_mat_general(e_vecs_init, "0"); 
    [~, ~, state_1_fixed_corr_qp] = ...
         get_tetron_init_cov_mat_general(e_vecs_init, "1");  
    % Comment: computation of state_0_fixed_corr_qp does not depend on the fixed basis,
    % eigenvectors as it is using a quasiparticle basis. 
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
    P_y1y2y3y4_v1(1) = calc_pfaffian_4x4(cov_mat_qp_curr(mzm_indices_cov_qp, mzm_indices_cov_qp));
    % Overlaps: 
    overlap_with_init_state(1) = calc_overlap_cov_mats(cov_mat_qp_init, cov_mat_qp_curr, 2*N);

    overlap_fixed_basis_states(1) = calc_overlap_cov_mats(cov_mat_qp_curr, state_0_fixed_cov_qp, 2*N) ...
        + calc_overlap_cov_mats(cov_mat_qp_curr, state_1_fixed_cov_qp, 2*N);
    
    % Instantaneous and initial bases same at first time step:
    overlap_inst_basis_states = overlap_fixed_basis_states; 
    overlap_inst_basis_states_deloc_MZM_debug = overlap_inst_basis_states; % DEBUG vector

    %Delocalised MZM code
    %--------%
    mu_init = mu_vec(:,1);
    H_tetron = get_tetron_BdG_Hamiltonian(mu, w, delta, N, BC, 'chain_1_chain_2'); 

    [e_vecs_deloc_init, e_vals_deloc_init, ~, comparison_determinants_init] = ... 
    diagonalise_uncoupled_tetron_via_Kitaev_chains_no_optim(H_tetron(1:2*N,1:2*N), ...
    H_tetron((2*N+1):end, (2*N+1):end), 'dirac_zero_modes', majorana_zero_modes_ref);
    
    comparison_determinants_mat = zeros(2, length(t_vec));
    comparison_determinants_mat(:,1) = comparison_determinants_init(:);

    % Store Eigenvalues: 'init_hamil_basis', 'instant_hamil_basis'.
    switch basis_choice
        case 'instant_hamil_basis'
            e_vals_mat = zeros(length(t_vec), 4*N);
            e_vals_mat(1,:) = e_vals_deloc_init(:).';
        case 'init_hamil_basis'
            e_vals_mat = NaN;
    end

    % Get initial quasiparticle occupation. Note that here instantaneous=
    % initial basis. 
    quasi_occ_t(1:N, 1) = diag(corr_mat_curr(1:N, 1:N));
    quasi_occ_t((1:N)+N, 1) = diag(corr_mat_curr((1:N)+2*N, (1:N)+2*N));  

    % Note: If time_evolution_method = 'time_ev_diagonalise', then I do
    % compute the eigenvalues even when basis_choice =
    % 'init_hamil_basis'. However their ordering may not be realiable,
    % since I just use "eig" in that case. 
            

    % if strcmp(basis_choice, 'init_hamil_basis') % is this right?
    %     comparison_determinants_mat(1,:) = comparison_determinants_init(1);
    %     comparison_determinants_mat(2,:) = comparison_determinants_init(2);        
    % end

    % End delocalised MZM code
    %-------%

    %[e_vecs_deloc_init, e_vals_deloc_init] = ... 
    %diagonalise_uncoupled_tetron_via_Kitaev_chains(H_tetron(1:2*N,1:2*N), ...
    %H_tetron((2*N+1):end, (2*N+1):end), 'dirac_zero_modes');

    corr_mat_curr_site = conj(e_vecs_init)*corr_mat_curr*e_vecs_init.'; 
    corr_mat_curr_deloc = e_vecs_deloc_init.'*corr_mat_curr_site*conj(e_vecs_deloc_init);
    cov_mat_curr_deloc = -1i.*omega_total*(2*corr_mat_curr_deloc - eye(4*N))*omega_total'; 
    P_y1y2y3y4_v2 = zeros(size(t_vec));
    P_y1y2y3y4_v2(1) = calc_pfaffian_4x4(cov_mat_curr_deloc(mzm_indices_cov_qp, mzm_indices_cov_qp));


    % % DEBUG
    % calc_overlap_cov_mats(cov_mat_qp_curr, state_0_fixed_cov_qp, 2*N)
    % calc_overlap_cov_mats(cov_mat_qp_curr, state_1_fixed_cov_qp, 2*N)
    % 
    % [~, ~, state_p_fixed_corr_qp] = ...
    %     get_tetron_init_cov_mat_general(e_vecs_init, "+");  
    % state_p_fixed_cov_qp = convert_corr_qp_to_cov_qp(state_p_fixed_corr_qp, omega_total);
    % [~, ~, state_m_fixed_corr_qp] = ...
    %     get_tetron_init_cov_mat_general(e_vecs_init, "-");  
    % state_m_fixed_cov_qp = convert_corr_qp_to_cov_qp(state_m_fixed_corr_qp, omega_total);
    % 
    % calc_overlap_cov_mats(state_p_fixed_cov_qp, state_0_fixed_cov_qp, 2*N)
    % calc_overlap_cov_mats(state_m_fixed_cov_qp, state_0_fixed_cov_qp, 2*N)       
    % calc_overlap_cov_mats(state_p_fixed_cov_qp, state_m_fixed_cov_qp, 2*N)   
    % calc_overlap_cov_mats(state_1_fixed_cov_qp, state_0_fixed_cov_qp, 2*N)
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
        
        if isscalar(delta_t)
            delta_t_curr = delta_t;
        else
            delta_t_curr = delta_t(ii-1);
        end

        switch time_evolution_method
            case 'time_ev_exp'
                U_delta_t = expm(1i*H_evolve.*delta_t_curr); 
            case 'time_ev_diagonalise'
                H_e = (H_evolve + H_evolve')/2;
                [V_QP, D_evolve] = eig(H_e);
                D_evolve = real(D_evolve);
                U_delta_t = V_QP*diag(exp(1i*diag(D_evolve).*delta_t_curr))*V_QP';  
                %U_delta_t = V_QP*diag(exp(1i*diag(D_evolve).*delta_t_curr))*V_QP';  
            case 'time_ev_diagonalise_debug'
                error('Option not implemented correctly');
                [V, D_evolve] = ... 
                    diagonalise_uncoupled_tetron_via_Kitaev_chains_no_optim(H_tetron(1:2*N,1:2*N), ...
                    H_tetron((2*N+1):end, (2*N+1):end), 'dirac_zero_modes', majorana_zero_modes_ref);
              %  [V, D_evolve] =...
              %      diagonalise_uncoupled_tetron_via_Kitaev_chains(H_tetron(1:2*N,1:2*N), ...
              %  H_tetron((2*N+1):end, (2*N+1):end), 'dirac_zero_modes');
                D_evolve = diag(D_evolve);
                V_QP = e_vecs_init'*V;
                
                D_evolve = real(D_evolve);
                U_delta_t = V_QP*diag(exp(1i*diag(D_evolve).*delta_t_curr))*V_QP'; 
            case 'time_ev_1st_order'
                H_e = (H_evolve + H_evolve')/2;
                U_delta_t = eye(4*N) + 1i*H_e*delta_t_curr;
            otherwise
                error('Invalid value for time_evolution_method');
        end

        %Time evolution
        corr_mat_next = U_delta_t*corr_mat_curr*U_delta_t';
        corr_mat_curr = corr_mat_next; % Recall this is in the initial QP basis
       
        %[e_vecs, e_vals, majorana_zero_modes] = ... 
        %diagonalise_uncoupled_tetron_via_Kitaev_chains(H_1, H_2, zero_mode_type)

        % Calculate quantities
        if strcmp(basis_choice, 'init_hamil_basis') 
            %Use basis of initial Hamiltonain
            cov_mat_qp_curr = -1i.*omega_total*(2*corr_mat_curr - eye(4*N))*omega_total'; 

            % Here I can extract my QP initial basis populations from:
            % corr_mat_curr
            quasi_occ_t(1:N, ii) = diag(corr_mat_curr(1:N, 1:N));
            quasi_occ_t((1:N)+N,ii) = diag(corr_mat_curr((1:N)+2*N, (1:N)+2*N));            

            % Ad-hoc fix for getting instantaneous eigenvalues:
             % could copy deloc majorana zero mode code. 

        elseif strcmp(basis_choice, 'instant_hamil_basis')
            % Get covariance matrix in basis of instantaneous Hamiltonian
            % eigenvectors
            % I still do time evolution using the initial basis. 

            % Convert correlation matrix in initial QP basis to site basis
            corr_mat_curr_site = conj(e_vecs_init)*corr_mat_curr*e_vecs_init.'; 
            if ~strcmp(time_evolution_method, 'time_ev_diagonalise')
                [V_QP, D_evolve] = eig(H_e); %In initial Hamiltonian eigenbasis
            end

            % Convert from the site basis to the instantaneous QP basis
            [e_vecs_deloc_curr, e_vals_deloc_curr, ~, comparison_determinants_curr] = ... 
                diagonalise_uncoupled_tetron_via_Kitaev_chains_no_optim(H_tetron(1:2*N,1:2*N), ...
                H_tetron((2*N+1):end, (2*N+1):end), 'dirac_zero_modes', majorana_zero_modes_ref);
            comparison_determinants_mat(:, ii) = comparison_determinants_curr(:);

            corr_mat_curr_deloc = e_vecs_deloc_curr.'*corr_mat_curr_site*conj(e_vecs_deloc_curr);
            cov_mat_qp_curr = -1i.*omega_total*(2*corr_mat_curr_deloc - eye(4*N))*omega_total'; 
    
            % Here I can extract my instantaneous QP basis populations from
            % corr_mat_curr_deloc
            quasi_occ_t(1:N, ii) = diag(corr_mat_curr_deloc(1:N, 1:N));
            quasi_occ_t((1:N)+N,ii) = diag(corr_mat_curr_deloc((1:N)+2*N, (1:N)+2*N));   

            % OLD CODE - to be more efficient you should avoid
            % diagonalising H_tetron twice. I didn't end up completing
            % that approach. 
            % % ORDER e_vecs_instant or V_QP properly - decreasing
            % % eigenvalues 
            % disp('ORDER e_vecs_instant properly');
            % % Convert diagonalised H_e back to site basis
            % e_vecs_instant = e_vecs_init*V_QP;  
            % % This^ has same eigenvalues, D_evolve
            % 
            % %Debug
            % if norm(e_vecs_instant*D_evolve*e_vecs_instant' - H_tetron) > 1e-10
            %     warning('Ive made a mistake');
            % end
            % Put correlation matrix in instanteneous Hamiltonian basis
            % THIS IS ONLY "CORRECT" IF e_vecs_instant has been ordered
            % properly - THAT IS SOMETHING I HAVE NOT DONE! 
            %corr_mat_curr_inst_qp = e_vecs_instant.'*corr_mat_curr_site*conj(e_vecs_instant);
             % DOUBLE CHECK NEXT LINE IS CORRECT
            %cov_mat_qp_curr = convert_corr_qp_to_cov_qp(corr_mat_curr_inst_qp, omega_total);
            

            % GO FROM HERE:----------------
            % Convert corr_mat_curr_inst_qp to covariance matrix
            % THEN DO CHECKS - make sure a noiseless tetron gets same
            % results for initial and instanteneous bases. 

        else
            error('Invalid Input for basis_choice');
        end

        %cov_mat_qp_curr = convert_corr_qp_to_cov_qp(corr_mat_curr, omega_total);   
        Z_exp(ii) = cov_mat_qp_curr(1,1+N); 
        X_exp(ii) = cov_mat_qp_curr(1,1+2*N);
        Y_exp(ii) = -cov_mat_qp_curr(1,1+3*N);
        P_y1y2y3y4_v1(ii) = calc_pfaffian_4x4(cov_mat_qp_curr(mzm_indices_cov_qp, mzm_indices_cov_qp));  

        % Error detected: if using instantaneous basis, then
        % cov_mat_qp_curr is in the instantaneous basis but cov_mat_qp_init
        % is in the initial basis. 
        overlap_with_init_state(ii) = calc_overlap_cov_mats(cov_mat_qp_init,cov_mat_qp_curr, 2*N);

        % Warning: you originaly were not careful about the bases used here -
        % It was equivalent to doing this in the instantaneous basis with 
        % delocalised MZMs (which I shouldn't be doing)
        % I beleive I fixed this. But remain vigilant.
        cov_mat_qp_init_basis_curr = -1i.*omega_total*(2*corr_mat_curr - eye(4*N))*omega_total'; 
        overlap_fixed_basis_states(ii) = calc_overlap_cov_mats(cov_mat_qp_init_basis_curr, state_0_fixed_cov_qp, 2*N) ...
            + calc_overlap_cov_mats(cov_mat_qp_init_basis_curr, state_1_fixed_cov_qp, 2*N);
        % overlap_fixed_basis_states(ii) = calc_overlap_cov_mats(cov_mat_qp_curr, state_0_fixed_cov_qp, 2*N) ...
        %     + calc_overlap_cov_mats(cov_mat_qp_curr, state_1_fixed_cov_qp, 2*N);
       % overlap_fixed_basis_states(ii) = NaN; % Add correct code

        %-------------------------------------------%
        %  Instantaneous computational basis state overlap
        %  I'm originally did this very innefficiently. I need the localised
        %  Majorana zero modes. I should do this without rediagonalising
        %  the entire tetron.
        %  Be careful about using the same basis in calc_overlap_cov_mats
        % I think I've shown analytically that for <0|psi> and <1|psi> I
        % can use the delocalised basis for |0> and the localised basis for
        % psi and still get the same result as using the localised basis
        % for both. Abhijeet agreed with this. So I've now switched to
        % using delocalised MZMs. But be careful with the
        % comparison_determinant sign

        state_0_inst_cov_qp = state_0_fixed_cov_qp;
        state_1_inst_cov_qp = state_1_fixed_cov_qp;

        overlap_inst_basis_states(ii) = calc_overlap_cov_mats(cov_mat_qp_curr, state_0_inst_cov_qp, 2*N) ...
            + calc_overlap_cov_mats(cov_mat_qp_curr, state_1_inst_cov_qp, 2*N);

        % % Old code which explicitly uses the e-vecs with the localised MZM
        % % states. 
        % if nargout >= 15
        %     % Do the calculation in in the instantaneous QP basis. Then:
        %     state_0_inst_cov_qp = state_0_fixed_cov_qp;
        %     state_1_inst_cov_qp = state_1_fixed_cov_qp;
        % 
        %     % Express current state using instantaneous LOCALISED MZM basis
        %     [e_vecs_loc_curr, e_vals_loc_curr, majorana_zero_modes_ref] = ... 
        %         diagonalise_uncoupled_tetron_via_Kitaev_chains(H_tetron(1:2*N,1:2*N), ...
        %         H_tetron((2*N+1):end, (2*N+1):end), 'dirac_zero_modes');
        % 
        %     corr_mat_curr_loc = e_vecs_loc_curr.'*corr_mat_curr_site*conj(e_vecs_loc_curr);
        %     cov_mat_qp_loc_curr = -1i.*omega_total*(2*corr_mat_curr_loc - eye(4*N))*omega_total'; 
        % 
        %     overlap_inst_basis_states(ii) = calc_overlap_cov_mats(cov_mat_qp_loc_curr, state_0_inst_cov_qp, 2*N) ...
        %         + calc_overlap_cov_mats(cov_mat_qp_loc_curr, state_1_inst_cov_qp, 2*N);
        % 
        %     % What if I used the delocalised MZMs? This is what I
        %     % previously stored in overlap_fixed_basis_states
        %     overlap_inst_basis_states_deloc_MZM_debug(ii) = calc_overlap_cov_mats(cov_mat_qp_curr, state_0_inst_cov_qp, 2*N) ...
        %         + calc_overlap_cov_mats(cov_mat_qp_curr, state_1_inst_cov_qp, 2*N);
        % end
        % ------------------------------%
    


        % DEBUG DELETE:
        % - - -%
        % if ii == length(t_vec)
        %     a =  calc_overlap_cov_mats(cov_mat_qp_curr, state_0_inst_cov_qp, 2*N) 
        %     b = calc_overlap_cov_mats(cov_mat_qp_loc_curr, state_0_inst_cov_qp, 2*N) 
        % 
        %     c =  calc_overlap_cov_mats(cov_mat_qp_curr, state_1_inst_cov_qp, 2*N) 
        %     d = calc_overlap_cov_mats(cov_mat_qp_loc_curr, state_1_inst_cov_qp, 2*N)   
        % 
        %     e = calc_pfaffian_4x4(cov_mat_qp_curr(mzm_indices_cov_qp, mzm_indices_cov_qp))
        %     f = calc_pfaffian_4x4(cov_mat_qp_loc_curr(mzm_indices_cov_qp, mzm_indices_cov_qp))
        % end
        %- - - %


        % [~, ~, state_0_fixed_corr_qp] = ...
        %      get_tetron_init_cov_mat_general(e_vecs_init, "0");
        % [~, ~, state_1_fixed_corr_qp] = ...
        %      get_tetron_init_cov_mat_general(e_vecs_init, "1");  
        % state_0_fixed_cov_qp = convert_corr_qp_to_cov_qp(state_0_fixed_corr_qp, omega_total);
        % state_1_fixed_cov_qp = convert_corr_qp_to_cov_qp(state_1_fixed_corr_qp, omega_total);        
        % 

        % FIX THIS: - below two lines are placeholders. Calculate the 0 and
        % 1 states using the instantaneous e_vecs (can't use delocalised
        % MZMs for that I think). 
        % state_0_inst_cov_qp = state_0_fixed_cov_qp;
        % state_1_inst_cov_qp = state_1_fixed_cov_qp;        
        % 
        % overlap_inst_basis_states(ii) = calc_overlap_cov_mats(cov_mat_qp_curr, state_0_inst_cov_qp, 2*N) ...
        %     + calc_overlap_cov_mats(cov_mat_qp_curr, state_1_cov_inst_qp, 2*N);

        %Delocalised MZM code - using init basis. Don't need this in my
        %final code, I reckon. It's just a check. Could copy this in there
        %commented out however. 
        corr_mat_curr_site = conj(e_vecs_init)*corr_mat_curr*e_vecs_init.'; 
        corr_mat_curr_deloc = e_vecs_deloc_init.'*corr_mat_curr_site*conj(e_vecs_deloc_init);
        cov_mat_curr_deloc = -1i.*omega_total*(2*corr_mat_curr_deloc - eye(4*N))*omega_total'; 
        P_y1y2y3y4_v2(ii) = calc_pfaffian_4x4(cov_mat_curr_deloc(mzm_indices_cov_qp, mzm_indices_cov_qp));

        if strcmp(basis_choice, 'instant_hamil_basis')
            e_vals_mat(ii, :) = e_vals_deloc_curr(:).';
        end

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
    P_y1y2y3y4_v1 = remove_negligible_imag_parts(P_y1y2y3y4_v1);
    overlap_with_init_state = remove_negligible_imag_parts(overlap_with_init_state);
    overlap_fixed_basis_states = remove_negligible_imag_parts(overlap_fixed_basis_states);
    overlap_inst_basis_states = remove_negligible_imag_parts(overlap_inst_basis_states);
    overlap_inst_basis_states_deloc_MZM_debug = remove_negligible_imag_parts(overlap_inst_basis_states_deloc_MZM_debug);

    %Delocalised MZM code
    P_y1y2y3y4_v2 = remove_negligible_imag_parts(P_y1y2y3y4_v2);

    % % Extract final QP populations - Don't use this, you change
    % corr_mat_deloc above.
    % if strcmp(basis_choice, 'init_hamil_basis') 
    %     % Here I can extract my QP initial basis populations from:
    %     % corr_mat_curr
    %     corr_mat_final = corr_mat_curr;
    % 
    % elseif strcmp(basis_choice, 'instant_hamil_basis')
    %     % Here I can extract my instantaneous QP basis populations from
    %     % corr_mat_curr_deloc
    %     corr_mat_final = corr_mat_curr_deloc;
    % end
    % quasi_occ_final(:,1) = diag(corr_mat_final(1:N, 1:N));
    % quasi_occ_final(:,2) = diag(corr_mat_final((1:N)+2*N, (1:N)+2*N));
    % quasi_occ_final = remove_negligible_imag_parts(quasi_occ_final);
    quasi_occ_t = remove_negligible_imag_parts(quasi_occ_t); 
    quasi_occ_final(:,1)  = quasi_occ_t(1:N, end);
    quasi_occ_final(:, 2) = quasi_occ_t((1+N):end, end);
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