%-------------------------------------------%
% Process Majorana Zero Modes Kitaev Chain Wrapper - No Optimisation
%-------------------------------------------%
% This function runs a simpler and faster version of
% process_majorana_zero_modes_kitaev_chain.m where the optimsation for the
% maximally localised majoranas has been removed. 
% This means the majoranas returned here are not guarenteed to be maximally
% localised, however the MZMs will be MZMs (not a mix of majorana and dirac
% fermions). 
%
% I do however include a check to ensure the parity of < 2i gamma_1 gamma_2>
% is consistent (there are no reflections, only rotations in MZM space). 
%
% I intend to use this at each time step in my time evolution for the
% instaneneous Pauli matrices basis. 
%
% I've included process_majorana_zero_modes_kitaev_chain_no_optim as a
% local function. Feel Free to Put it in its own file. 
% ------------------------------------------------------%
% ------------------------------------------------------%
%
% DOCUMENTATION FOR ORIGINAL VERSION OF THIS FUNCTION (with optimatisation):
% This is a basic wrapper function for
% process_majorana_zero_modes_kitaev_chain. Consider the matrix V, obtained 
% from [V, D] = eig(H_BdG) - where H_BdG is a BdG Hamiltonian for a Kitaev 
% chain with at most two majorana modes. 
% The purpose of this wrapper is to replace the two 
% eigenvectors in the near-degenerate zero energy manifold with either the
% two majoranas that are localised on the left and right hand sides of the
% chain or the quasiparticle creation and annihilation operators which
% corresponds to superpositions of these two modes. 
%
% To do: add checks on the e_vals to check whether process_majorana_zero_modes_kitaev_chain
% should be called. At present this is just done blindly. 
% To do: modify e_vals out if two majoranas are returned which are not
% strictly eigenvectors (because |d_0> and |d_0^dag> are not degenerate). 
% To do: add additional input: thresh - if the absolute value of the lowest
% two eigenergies are below this value, then do nothing. 
% 
% Inputs:
%         e_vecs - eigenvectors of a BdG Hamiltonian
%         e_vals - eigenvalues of a BdG Hamiltonian
%         majorana_zero_modes_ref - 4N x 2 matrix with e_vecs corresponding to 
%               reference majorana zero modes. Output zero modes with be checked 
%               to only be an approximate rotation from these in the MZM subspace and
%               not a reflection. First vector is the positive energy
%               BDG vector and second vector has negative (CHECK) energy vector. 
%
% Ouptut: e_vecs_out - processed eigenvectors where any near zero
%         energy eigenvectors have been rotated to majorana zero modes with
%         maximimsed weight each on on side of the chain or the
%         corresponding dirac fermion modes.
%         e_vals_out - e_vals adjusted to correspond to processed. These
%         are in descending order. 
%         e_vecs matrix.
%         zero_mode_type - either majorana_zero_modes or dirac_zero_modes.
%         

% If the MZM splitting is non-zero then the outputs are not technically all
% eigenvectors - for zero_mode_type -> majorana_zero_modes

function [e_vecs_out, e_vals_out, majorana_zero_modes, comp_det] = ...
    process_majorana_zero_modes_kitaev_chain_wrapper_no_optim(e_vecs, e_vals, zero_mode_type, majorana_zero_modes_ref)
 %   if nargin > 4
 %       error('Too many input parameters');
 %   end
 %   if nargin == 4
 %       order = varargin{1}
 %   else
 %       order = 'descend';
 %   end

    if mod(size(e_vecs,2),2) || mod(size(e_vecs,1),2) || ...
         size(e_vecs,1) ~= size(e_vecs,2)
        error('e_vecs must be a 2N x 2N matrix');
    end
    if mod(size(e_vals,2),2) && mod(size(e_vals,1),2)
        error('e_vals must be a length 2N vector');
    end
    if ~strcmp(zero_mode_type, 'majorana_zero_modes') && ...
            ~strcmp(zero_mode_type, 'dirac_zero_modes') 
        error('Invalid input for zero_mode_type');
    end
    N = length(e_vals)./2;
    [e_vecs, e_vals] = sort_eigenvectors_and_eigenvalues(e_vecs, e_vals, 'descend');
    e_vecs_out = e_vecs; e_vals_out = e_vals;
    zero_modes = e_vecs(:, [N, N+1]);
    zero_modes_energies = e_vals([N,N+1]);
    
    [majorana_zero_modes, deloc_dirac_zero_mode, comp_det] = ...
    process_majorana_zero_modes_kitaev_chain_no_optim(zero_modes,...
    zero_modes_energies, majorana_zero_modes_ref);

    if strcmp(zero_mode_type, 'majorana_zero_modes')
        e_vecs_out(:,[N,N+1]) = majorana_zero_modes;
    end
    if strcmp(zero_mode_type, 'dirac_zero_modes')
        e_vecs_out(:,[N,N+1]) = deloc_dirac_zero_mode;
    end   

   % disp('Add logic for near-deg subspace threshold in process MZM wrapper');
    %disp('Also add logic for e_vals_out');
end

%------------------------------------------------------------------------%
%------------------------------------------------------------------------%
%-------------------------------------------%
% Process Majorana Zero Modes Kitaev Chain
%-------------------------------------------%
% This function runs a simpler and faster version of
% process_majorana_zero_modes_kitaev_chain.m where the optimsation for the
% maximally localised majoranas has been removed. 
% This means the majoranas returned here are not guarenteed to be maximally
% localised, however the MZMs will be MZMs (not a mix of majorana and dirac
% fermions). 
%
% I do however include a check to ensure the parity of < 2i gamma_1 gamma_2>
% is consistent (there are no reflections, only rotations in MZM space). 
%
% I intend to use this at each time step in my time evolution for the
% instaneneous Pauli matrices basis. 
%
%
% ------------------------------------------------------%
% ------------------------------------------------------%
%
% DOCUMENTATION FOR ORIGINAL VERSION OF THIS FUNCTION (with optimatisation):
%
% Inputs: zero_modes - 2N x 2 matrix where columns are the eigenvectors of
%         a BdG Hamiltonian with near zero energy. 
%         zero_modes_energies - 1 x 2N with corresponding energies of
%         zero_modes. 
% Outputs: majorana_zero_modes - 2N x 2 matrix where the columns are the
%      rotated eigenvectors of zero_modes in the Majorana zero mode basis. 
%          deloc_dirac_zero_mode - 2N x 2 matrix containing |d_0> and |d_0^dag>
%                             where |d_0> is the delocalised fermionic mode
%                             corresponding to the majorana_zero_modes. 


function [majorana_zero_modes, deloc_dirac_zero_mode, comp_det] = ...
    process_majorana_zero_modes_kitaev_chain_no_optim(zero_modes, zero_modes_energies, majorana_zero_modes_ref)
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

    %% Get signs consistent with reference gamma_1 and gamma_2
    % maj_deloc_1;
    % maj_deloc_2; 
    
    % Get reference majorana modes
    % zero_modes_ref
    %gamma_1_ref = 1/sqrt(2).*(zero_modes_ref(1) + zero_modes_ref(2));
    %gamma_2_ref = 1i/sqrt(2).*(zero_modes_ref(1) - zero_modes_ref(2));
    gamma_1_ref = majorana_zero_modes_ref(:,1);
    gamma_2_ref = majorana_zero_modes_ref(:,2);


    % Construct the comparision matrix
    comp_matrix = [maj_deloc_1'*gamma_1_ref, maj_deloc_2'*gamma_1_ref;...
                   maj_deloc_1'*gamma_2_ref, maj_deloc_2'*gamma_2_ref];
    comp_det = det(comp_matrix);

    gamma_1 = maj_deloc_1; 
    if sign(real(comp_det)) == 1
        gamma_2 = maj_deloc_2;
    else
        gamma_2 = -maj_deloc_2;
            % Reconstruct the comparison matrix
        comp_matrix_new = [gamma_1'*gamma_1_ref, gamma_2'*gamma_1_ref;...
                   gamma_1'*gamma_2_ref, gamma_2'*gamma_2_ref];
        if sign(real(det(comp_matrix_new))) == -1
            warning('Comparison matrix method incorrect.');
        end
    end

    if abs(imag(comp_det)./real(comp_det)) >1e-5
        warning('Comparison matrix method incorrect');
    end

    majorana_zero_modes = [gamma_1, gamma_2];       
    deloc_dirac_zero_mode = zeros(size(majorana_zero_modes));
    %deloc_dirac_zero_mode(:,1) = 1/sqrt(2)*(gamma_1 + 1i*gamma_2);
    %deloc_dirac_zero_mode(:,2) = 1/sqrt(2)*(gamma_1 - 1i*gamma_2);
    
    % CHANGE MADE ON THE 27_06_24 (I switched the columns of
    % deloc_dirac_zero_mode
    deloc_dirac_zero_mode(:,1) = 1/sqrt(2)*(gamma_1 - 1i*gamma_2);
    deloc_dirac_zero_mode(:,2) = 1/sqrt(2)*(gamma_1 + 1i*gamma_2);
      
end


    
            % I DON'T THINK THIS IS NECESSARY:
            % % It's possible that gamma_1_temp or gamma_2_temp are off by a global 
            % % factor of "-". It doesn't matter if they both have a global factor of 
            % % "-", but if only one is, this corresponds to switching |d0> and |d0^dag>.
            % % We could live with this, but I'd rather fix it. 
            % 
            % dirac_zero_mode_temp = 1/sqrt(2).*[(gamma_1_temp + 1i*gamma_2_temp),...
            %     gamma_1_temp - 1i*gamma_2_temp]; 
            % test_neg_overlap = abs(e_vec_neg'*dirac_zero_mode_temp);
            % if test_neg_overlap(1) >= test_neg_overlap(2)
            %     gamma_2 = gamma_2_temp;
            % else
            %     gamma_2 = -gamma_2_temp;
            % end    
            % if max(test_neg_overlap) <0.98
            %     test_neg_overlap_failed_flag = 1;
            % else
            %     test_neg_overlap_failed_flag = 0;
            % end
            % 
            % gamma_1 = gamma_1_temp; 


%% ------------------- %%
%    Local Functions
%-----------------------%

%% Particle Hole Transformation

function vec_PHT = particle_hole_trans(vec)
    N = length(vec)/2;
    X = [0 1; 1 0]; 
    vec_PHT = kron(X, eye(N))*conj(vec);
end








