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


function [cov, corr_site_0, corr_site_1, corr_qp_0, corr_qp_1] =...
    get_tetron_init_cov_mat_general(e_vecs)
    % Add Input checks
    
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

   % Convert correlation matrices to site basis (dirac fermions).
   X = [0 1; 1 0];
   S = kron(X, eye(N));
   T = kron(eye(2),S)*e_vecs*kron(eye(2),S);
   corr_site_0_dep = T*corr_qp_0*T';
   corr_site_1_dep = T*corr_qp_1*T';   
   
   corr_site_0 = conj(e_vecs)*corr_qp_0*e_vecs.';
   corr_site_1 = conj(e_vecs)*corr_qp_1*e_vecs.';
   
   check = norm(corr_site_0 - corr_site_0_dep); % Delete
   
   lam = kron(eye(2),[fliplr(eye(N)), zeros(N); zeros(N), eye(N)]); % Delete
   T2 = kron(eye(2),S)*e_vecs*lam*kron(eye(2),S); % Delete
   corr_site_1_check = T2*corr_qp_1*T2';   % Delete  
   norm( corr_site_1_check + corr_site_1) % Delete
    
   % Convert to covariance matrices in site basis (majorana fermions). 
   omega_s = 1/sqrt(2)*[eye(N),eye(N); 1i*eye(N), -1i*eye(N)]; 
   cov_0 = -1i*kron(eye(2),omega_s)*(2*corr_site_0 - eye(4*N))*kron(eye(2),omega_s');
   cov_1 = -1i*kron(eye(2),omega_s)*(2*corr_site_1 - eye(4*N))*kron(eye(2),omega_s');  
   
   %Output both correlation and covariance matrices? or just covariance. 
   warning('I have not checked cov_1, corr_site_1 rigorously');   
end



