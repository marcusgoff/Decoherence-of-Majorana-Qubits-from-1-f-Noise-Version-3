%-------------------------------------------%
% Get Tetron Initalised Covariance Matrix
%-------------------------------------------%
%
% computes the covariance matrices for the computational basis states of the
% tetron:
%      |0> = |even>|even>
%      |1> = |odd>|odd>
%
% INPUTS:
%      e_vecs - eigenvectors of the tetron qubit in the chain_1-chain_2
%      basis. For each chain the eigenvalues must be ordered in descending
%      order. 
%      state - allowed inputs: "+", "-", "0", "1"
%      
% OUTPUTS
%     cov_site - covariance matrix in the site majorana basis
%     corr_site - correlation matrix in the site dirac fermion basis
%     corr_qp - correlation matrix in the dirac fermion quasiparticle basis
%
% Comments:
%     corr_qp is in the quasiparticle basis, and so it does not depend on
%     the e_vecs.

function [cov_site, corr_site, corr_qp] =...
    get_tetron_init_cov_mat_general(e_vecs, state)
    
   check_inputs(state); %only checks state at present.   
   
   N = size(e_vecs,1)./4;
    
   cov_0 = zeros(4*N);
   cov_1 = zeros(4*N);
     
   % Create correlation matrices in QP basis (dirac fermions).
   A = zeros(N); A(N,N) = 1; 
   B = eye(N); B(1,1) = 0; 
   
   corr_qp_KC_even = [zeros(N), zeros(N); zeros(N), eye(N)];
   corr_qp_KC_odd = [A, zeros(N); zeros(N), B];
   
   corr_qp_0 = [corr_qp_KC_even, zeros(2*N); zeros(2*N), corr_qp_KC_even];
   corr_qp_1 = [corr_qp_KC_odd, zeros(2*N); zeros(2*N), corr_qp_KC_odd];   
   
   % Create off-diagonal block matrices for correlations in +/- state:
   C = zeros(2*N);
   C(N+1, N) = 1i; C(N, N+1) = 1i;  
   gamma_plus_off_diag = [zeros(2*N), C; -C, zeros(2*N)];
   %This matrix, gamma_plus_off_diag is the off-diagonal terms (x 2) for the "+" state. The
   %"-" state simply has the negative of the off-diagonal terms in the "-"
   %state.
      
   switch state
       case "0"
          corr_qp = corr_qp_0;            
       case "1"
          corr_qp = corr_qp_1;
       case "+"
          corr_qp = 0.5*(corr_qp_0 + corr_qp_1 + gamma_plus_off_diag);
       case "-"
          corr_qp = 0.5*(corr_qp_0 + corr_qp_1 - gamma_plus_off_diag);
       case "+i"
           error('case not yet coded');
         % corr_qp = 0.5*(corr_qp_0 - corr_qp_1 + 1i*gamma_plus_off_diag);           
       case "-i"
           error('case not yet coded');
          %corr_qp = 0.5*(corr_qp_0 - corr_qp_1 - 1i*gamma_plus_off_diag); 
   end      
   
   % Convert correlation matrix to site basis (dirac fermions):  
   corr_site = conj(e_vecs)*corr_qp*e_vecs.';
    
   % Convert to covariance matrices in site basis (majorana fermions). 
   omega_s = 1/sqrt(2)*[eye(N),eye(N); 1i*eye(N), -1i*eye(N)]; 
   cov_site = -1i*kron(eye(2),omega_s)*(2*corr_site - eye(4*N))*kron(eye(2),omega_s'); 
   
   %Output both correlation and covariance matrices? or just covariance. 
  % warning('I have not checked cov_1, corr_site_1 rigorously');   
end

function check_inputs(state)
    if ~(strcmp(state, '0') || strcmp(state, '1') || strcmp(state, '+') || ...
            strcmp(state, '-') || strcmp(state, '+i') || strcmp(state, '-i'))
        error('Invalid input for parameter ''state''');
    end
end

