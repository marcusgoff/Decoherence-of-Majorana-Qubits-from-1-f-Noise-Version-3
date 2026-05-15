%-------------------------------------------%
% Process Majorana Zero Modes Kitaev Chain Wrapper
%-------------------------------------------%
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

function [e_vecs_out, e_vals_out, majorana_zero_modes, dirac_zero_modes] = ...
    process_majorana_zero_modes_kitaev_chain_wrapper(e_vecs, e_vals, zero_mode_type)
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
    
    [majorana_zero_modes, deloc_dirac_zero_mode] = ...
    process_majorana_zero_modes_kitaev_chain(zero_modes, zero_modes_energies);
    if strcmp(zero_mode_type, 'majorana_zero_modes')
        e_vecs_out(:,[N,N+1]) = majorana_zero_modes;
    end
    if strcmp(zero_mode_type, 'dirac_zero_modes')
        e_vecs_out(:,[N,N+1]) = deloc_dirac_zero_mode;
    end   

   % disp('Add logic for near-deg subspace threshold in process MZM wrapper');
    %disp('Also add logic for e_vals_out');
end









