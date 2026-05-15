%-------------------------------------------------%
% Tetron Time Evolution - Site Basis Correlation
%-------------------------------------------------%
% Local Utility function. Might be a good idea to write a more generalized 
% version of this later. 
%
% Notes on inputs (not comprehensive): 
% INPUTS:
%          e_vecs - tetron eigenvectors, ordered by chain-1 chain-2 first, 
%                   and by descending eigenvalues second. 
%                   This is achieved by using diagonalise_uncoupled_tetron_via_Kitaev_chains
%                   function. Make sure you use this ordering, or this
%                   function will return an error. 
%          e_vecs_init - used to define the Pauli operators and the basis
%          for corr_mat_init. 

function [Z_exp, X_exp, Y_exp, t_vec, first_site_occ] =...
    calculate_tetron_evolution_dirac_qp_corr_mat(...
    corr_mat_init, e_vecs_init, e_vals_init, H_evolve, delta_t, t_init, t_final)

    %check_inputs(e_vecs_init, e_vals_init) - not good function.

    % I PREVIOULY USED THIS:

    %U_delta_t = expm(1i*H_evolve.*delta_t);  


    % INSTEAD OF THIS
    H_e = (H_evolve + H_evolve')/2;
    [V, D_evolve] = eig(H_e);
    D_evolve = real(D_evolve);
    U_delta_t = V*diag(exp(1i*diag(D_evolve).*delta_t))*V';


%     [V,D_evolve] = eig(H_evolve);
%     if norm(imag(D_evolve)) >1e-5*norm(real(D_evolve)) 
%        warning('Non-negligble imaginary part removed'); 
%     end
%     D_evolve = real(D_evolve);
%     U_delta_t = V*diag(exp(1i*diag(D_evolve)*delta_t))*V';    
    
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
    
    if nargin > 5
        corr_mat_site_curr = conj(e_vecs_init)*corr_mat_curr*e_vecs_init.';
        first_site_occ(1) = corr_mat_site_curr(1,1);
    end
    
    for ii = 2:length(t_vec)
        %Time evolution
        corr_mat_next = U_delta_t*corr_mat_curr*U_delta_t';
        corr_mat_curr = corr_mat_next; 
        % Calculate quantities
        cov_mat_qp_curr = convert_corr_qp_to_cov_qp(corr_mat_curr);   
        Z_exp(ii) = cov_mat_qp_curr(1,1+N); 
        X_exp(ii) = cov_mat_qp_curr(1,1+2*N);
        Y_exp(ii) = -cov_mat_qp_curr(1,1+3*N);
        
        if nargin> 5
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