%-------------------------------------------------%
% Tetron Time Evolution - QP Correlation - WITH NOISE
%-------------------------------------------------%
% Local Utility function. Might be a good idea to write a more generalized 
% version of this later. 
% This function is a bit different to its predecessors. Instead of
% specifying an H_evolve, here you must specificy tetron parameters and the
% noise trajectory for the mu_vec. 
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
%
% TO DO: CHANGE e_vects_init an e_vals_init to e_vecs_av e_vals_av or
% noiseless. THey are not the initial e-vecs and e-vals, they are the that
% of the tetron without noise - or the long time averaged tetron.

function [Z_exp, X_exp, Y_exp, t_vec, first_site_occ, mu_vec] =...
    calculate_noisy_tetron_evolution_dirac_qp_corr_mat(...
    corr_mat_init, e_vecs_init, e_vals_init, tetron_params, mu_vec, delta_t, t_init, t_final)

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
 
    mu_curr = mu_vec(:,1); %2 x 1 array (top and bottom chain chemical potentials)
    H_tetron = get_tetron_BdG_Hamiltonian(mu_curr, w, delta, N, BC, 'chain_1_chain_2'); 
    H_evolve = e_vecs_init'*H_tetron*e_vecs_init; %in initial QP basis
    
    %U_delta_t = expm(1i*H_evolve.*delta_t);  

    N = length(e_vals_init)./4;
    
    % Run time evolution
    t_vec = t_init:delta_t:t_final;
    Z_exp = zeros(size(t_vec));
    X_exp = zeros(size(t_vec));
    Y_exp = zeros(size(t_vec));
    first_site_occ = zeros(size(t_vec));
    
    corr_mat_curr = corr_mat_init; % In QP basis
    cov_mat_qp_curr = convert_corr_qp_to_cov_qp(corr_mat_curr);   
    Z_exp(1) = cov_mat_qp_curr(1,1+N);
    X_exp(1) = cov_mat_qp_curr(1,1+2*N);
    Y_exp(1) = -cov_mat_qp_curr(1,1+3*N); 
    
    if nargout >= 5 % THIS SHOULD BE NARGOUT and >=5
        corr_mat_site_curr = conj(e_vecs_init)*corr_mat_curr*e_vecs_init.';
        first_site_occ(1) = corr_mat_site_curr(1,1);
    end
    
    for ii = 2:length(t_vec)
        % Get current time evolution unitary
        mu_curr = mu_vec(:,ii);
        H_tetron = get_tetron_BdG_Hamiltonian(mu_curr, w, delta, N, BC, 'chain_1_chain_2'); 
        H_evolve = e_vecs_init'*H_tetron*e_vecs_init; %in initial QP basis   
        
        U_delta_t = expm(1i*H_evolve.*delta_t);  
        % H_e = (H_evolve + H_evolve')/2;
        % [V, D_evolve] = eig(H_e);
        % D_evolve = real(D_evolve);
        % U_delta_t = V*diag(exp(1i*diag(D_evolve).*delta_t))*V';

        
        %Time evolution
        corr_mat_next = U_delta_t*corr_mat_curr*U_delta_t';
        corr_mat_curr = corr_mat_next; 
        % Calculate quantities
        cov_mat_qp_curr = convert_corr_qp_to_cov_qp(corr_mat_curr);   
        Z_exp(ii) = cov_mat_qp_curr(1,1+N); 
        X_exp(ii) = cov_mat_qp_curr(1,1+2*N);
        Y_exp(ii) = -cov_mat_qp_curr(1,1+3*N);
        
        if nargout >= 5
            corr_mat_site_curr = conj(e_vecs_init)*corr_mat_curr*e_vecs_init.';
            first_site_occ(ii) = corr_mat_site_curr(1,1);
        end    
    end
    corr_mat_final = corr_mat_curr;    
    % warning('Check minus sign in converting corr to cov mat - should be fine.');
    %Dodgely remove imaginary part to Z_exp if its negligble
    
%     if max(abs(imag(Z_exp))) < 1e-10 & ...
%             norm(imag(Z_exp))/norm(real(Z_exp)) < 1e-8
%         Z_exp = real(Z_exp);
%     end
    Z_exp = remove_negligble_imag_parts(Z_exp);
    X_exp = remove_negligble_imag_parts(X_exp);
    Y_exp = remove_negligble_imag_parts(Y_exp);
end

function cov_mat_qp = convert_corr_qp_to_cov_qp(corr_mat_qp)
    N = size(corr_mat_qp,1)./4;
    
    %Conversion matrices - inefficient recreating these all the time.
    % would be better to hold in memory.
    omega_s = 1/sqrt(2)*[eye(N),eye(N); 1i*eye(N), -1i*eye(N)]; 
    lam = [fliplr(eye(N)), zeros(N); zeros(N), eye(N)]; 
    omega_total = kron(eye(2),omega_s*lam);
    
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

function y = remove_negligble_imag_parts(x)

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