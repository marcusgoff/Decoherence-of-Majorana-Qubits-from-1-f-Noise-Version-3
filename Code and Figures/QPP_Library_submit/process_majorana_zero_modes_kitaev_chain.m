%-------------------------------------------%
% Process Majorana Zero Modes Kitaev Chain
%-------------------------------------------%
% Inputs: zero_modes - 2N x 2 matrix where columns are the eigenvectors of
%         a BdG Hamiltonian with near zero energy. 
%         zero_modes_energies - 1 x 2N with corresponding energies of
%         zero_modes. 
% Outputs: majorana_zero_modes - 2N x 2 matrix where the columns are the
%      rotated eigenvectors of zero_modes in the Majorana zero mode basis. 
%          deloc_dirac_zero_mode - 2N x 2 matrix containing |d_0> and |d_0^dag>
%                             where |d_0> is the delocalised fermionic mode
%                             corresponding to the majorana_zero_modes. 

% FOR THE WRAPPER FUNCTION
% I change my mind - not going to do this - I think it's a good idea to 
% create a wrapper function to do this.
% Inputs:
%         e_vecs - eigenvectors of a BdG Hamiltonian
%         e_vals - eigenvalues of a BdG Hamiltonian
% Ouptut: e_vecs_processed - processed eigenvectors where any near zero
%         energy eigenvectors have been rotated to majorana zero modes with
%         maximimsed weight each on on side of the chain
%         e_vals_processed - e_vals adjusted to correspond to processed
%         e_vecs matrix.
%         test_neg_overlap_failed_flag - OPTIONAL output. If this is raised
%         |d_0> may have a positive energy and |d_0^dag> may have a
%         negative energy. This will often occur when you are closed to the 
%         fixed point of the Kitaev chain. You likely can work with this.
%
%
% UPDATES MADE:
%        27_06_204 - switched columns of deloc_dirac_zero_mode, so that
%        e_vecs correspond correctly to energies in zero_mode_energies. 

function [majorana_zero_modes, deloc_dirac_zero_mode, test_neg_overlap_failed_flag] = ...
    process_majorana_zero_modes_kitaev_chain(zero_modes, zero_modes_energies)
    %% Process and Check Inputs
    if size(zero_modes,2) ~= 2 || mod(size(zero_modes,1),2)
        error('zero_modes must be a 2N x 2 matrix');
    end
    if size(zero_modes_energies,2) ~= 1 && size(zero_modes_energies,1) ~= 1|| ...
            length(zero_modes_energies) ~= 2
        error('zero_modes_energy must be a 1 x 2 matrix');
    end   
    % Check whether we already have fermions or majorana fermions
    % TO DO
    
    [e_vecs, e_vals] = sort_eigenvectors_and_eigenvalues(...
        zero_modes, zero_modes_energies, 'descend');
    e_vec_pos = zero_modes(:,1); % positive energy
    e_vec_neg = zero_modes(:,2); % negative energy
    
    %% Locally Defined Variables
    norm_zero_thresh = 1e-10;
    N = length(zero_modes)./2; 
    X = [0 1; 1 0];  
    
    %% Define a delocalised majorana mode from e_vec_neg
    maj_deloc_1 = e_vec_neg + particle_hole_trans(e_vec_neg); % delocalised majorana mode
    
    norm_sq_maj_deloc = maj_deloc_1'*maj_deloc_1; 
    if norm_sq_maj_deloc < norm_zero_thresh
       %warning('Check 1 alternate passage');
       maj_deloc_1  = 1i*(e_vec_neg - particle_hole_trans(e_vec_neg)); 
       norm_sq_maj_deloc = maj_deloc_1'*maj_deloc_1; 
       if norm_sq_maj_deloc < norm_zero_thresh
          error('Need to change algorithm for creating delocalised majorana mode'); 
       end
    end
    maj_deloc_1 = maj_deloc_1./sqrt(norm_sq_maj_deloc);
    
    %% Define maj_deloc_2 as the orthogonal eigenstate in near-zero energy subspace
    % Define projector onto near-zero energy subspace
    proj_near_zero = e_vec_neg*e_vec_neg' + e_vec_pos*e_vec_pos';
    if norm(proj_near_zero^2 - proj_near_zero)./(4*N^2) > 1e-10
        error('Error in proj_near_zero calculation');
    end
    proj_maj_deloc_2 = proj_near_zero - maj_deloc_1*maj_deloc_1';
    [V,D] = eig(proj_maj_deloc_2, 'vector');
    [~, ind] = min(abs(D - 1));
    if abs(D(ind)-1) > 1e-5
        error('Cannot find maj_deloc_2 using projector method');
    end
    maj_deloc_2_guess = V(:,ind);
    % Ensure maj_deloc_2 is Hermitian (it could have come with a factor of
    % exp(1i*theta):
    maj_deloc_2 = maj_deloc_2_guess + particle_hole_trans(maj_deloc_2_guess);

    norm_sq_maj_2_deloc = maj_deloc_2'*maj_deloc_2; 
    if norm_sq_maj_2_deloc < norm_zero_thresh
       maj_deloc_2  = 1i*(maj_deloc_2_guess - particle_hole_trans(maj_deloc_2_guess)); 
       norm_sq_maj_2_deloc = maj_deloc_2'*maj_deloc_2; 
       if norm_sq_maj_2_deloc < norm_zero_thresh
          error('Need to change algorithm for creating delocalised majorana mode'); 
       end
    end
    maj_deloc_2 = maj_deloc_2./sqrt(norm_sq_maj_2_deloc);
    
    
    %% Find Majorana operators with weight on either side of chain
    % Find superpositions of maj_deloc_1 and maj_deloc_2 which minimise
    % expectation of position operator
    position_op = diag([1:N, 1:N]);
    theta_checks = linspace(-pi,pi,1e5);
    theta_min_val = 0;
    min_val = Inf;
    exp_vals = zeros(size(theta_checks)); % DEBUG
    
    for ii = 1:length(theta_checks)
       gamma_check = cos(theta_checks(ii))*maj_deloc_1 + sin(theta_checks(ii))*maj_deloc_2;
       exp_position_op =  real(gamma_check'*position_op*gamma_check);
       exp_vals(ii) = exp_position_op; % DEBUG
       if exp_position_op < min_val
           min_val = exp_position_op;
           theta_min_val = theta_checks(ii);
       end
    end
    gamma_1_temp = cos(theta_min_val)*maj_deloc_1 + sin(theta_min_val)*maj_deloc_2;
    gamma_2_temp =  cos(theta_min_val + pi/2)*maj_deloc_1 + sin(theta_min_val + pi/2)*maj_deloc_2; 
    
    % It's possible that gamma_1_temp or gamma_2_temp are off by a global 
    % factor of "-". It doesn't matter if they both have a global factor of 
    % "-", but if only one is, this corresponds to switching |d0> and |d0^dag>.
    % We could live with this, but I'd rather fix it. 
    
    dirac_zero_mode_temp = 1/sqrt(2).*[(gamma_1_temp + 1i*gamma_2_temp),...
        gamma_1_temp - 1i*gamma_2_temp]; 
    
    test_neg_overlap = abs(e_vec_neg'*dirac_zero_mode_temp);

    if test_neg_overlap(1) >= test_neg_overlap(2)
        gamma_2 = gamma_2_temp;
    else
        gamma_2 = -gamma_2_temp;
    end    
    if max(test_neg_overlap) <0.98
        test_neg_overlap_failed_flag = 1;
    else
        test_neg_overlap_failed_flag = 0;
    end
    
    gamma_1 = gamma_1_temp; 
  
    majorana_zero_modes = [gamma_1, gamma_2];       
    deloc_dirac_zero_mode = zeros(size(majorana_zero_modes));
    %deloc_dirac_zero_mode(:,1) = 1/sqrt(2)*(gamma_1 + 1i*gamma_2);
    %deloc_dirac_zero_mode(:,2) = 1/sqrt(2)*(gamma_1 - 1i*gamma_2);
    
    % CHANGE MADE ON THE 27_06_24 (I switched the columns of
    % deloc_dirac_zero_mode
    deloc_dirac_zero_mode(:,1) = 1/sqrt(2)*(gamma_1 - 1i*gamma_2);
    deloc_dirac_zero_mode(:,2) = 1/sqrt(2)*(gamma_1 + 1i*gamma_2);
      
end



%% ------------------- %%
%    Local Functions
%-----------------------%

%% Particle Hole Transformation

function vec_PHT = particle_hole_trans(vec)
    N = length(vec)/2;
    X = [0 1; 1 0]; 
    vec_PHT = kron(X, eye(N))*conj(vec);
end

%% Check Eigenvalue

% function [e_val, max_nonzero_value] = check_eigenvalue(A,vec)
%     thresh_1 = 1e-5;
%     thresh_2 = 1e-5;
%     product_term = A*vec./vec;
%     potential_e_vals = product_term(abs(vec) > thresh_1);
%     max_nonzero_value = max(abs(potential_e_vals));
%     if std(potential_e_vals) < thresh_2
%         e_val = mean(potential_e_vals);
%     else
%         e_val = NaN
%     end
%     check_near_zero_eigenvalue(A,vec);
% end
% 
% %% Check Near Zero Eigenvalues
% % Returns the largest element of abs(A*vec./vec) - ignoring near-zero
% % elements of vec.
% 
% function max_nonzero_value = check_near_zero_eigenvalue(A,vec)
%     thresh_1 = 1e-5;
%     product_term = A*vec./vec;
%     potential_e_vals = product_term(abs(vec) > thresh_1);
%     max_nonzero_value = max(abs(potential_e_vals));
% end
% 
