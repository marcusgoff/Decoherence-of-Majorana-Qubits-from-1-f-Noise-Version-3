%-------------------------------------------%
% Get Tetron Initalised Covariance Matrix - Check
%-------------------------------------------%
% INCOMPLETE FUNCTION - I DECIDED THAT THIS IS NOT NECESSARY
%
% Here I directly calculate the Covariance matrix in the site basis -
% primarly using this method as a check for
% convert_correlation_to_covariance_mat_check.
% Can only be used for a tetron that is initally uncoupled
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
%      

function [cov_0, cov_1, corr_site_0] = get_tetron_init_cov_mat_check(e_vecs)

    cov_0 = zeros(size(e_vecs));
    cov_1 = zeros(size(e_vecs));

    N = size(e_vecs,1)./4;
    % Here I'm going to straight up assume that e_vecs is straight up
    % uncoupled.
    U_1 = e_vecs(1:N,1:N);
    U_2 = e_vecs(N+1:end,N+1:end);
    
    
    %-----------------%
    % EVEN-EVEN STATE
    %------------------%
    % Solve cov_0 = |even>|even> which is the vacuum of all quasiparticles
    % stored in U_1 and U_2. 
    % Get exp values of quadratic terms of site-local operators (in the
    % even-even state):
    alpha_1 = U_1(1:N, N+1:end);
    beta_1  = U_1(N+1:end,N+1:end);
    
    ann_ann_1 = zeros(2*N,2*N);
    ann_cr_1 = zeros(2*N,2*N);
    cr_ann_1 = zeros(2*N,2*N);
    cr_cr_1 = zeros(2*N,2*N);
    
    
    
    for ii = 1:N
        for jj = 1:N
            for m = 1:N
            % ENTER CODE HERE.
                %corr_mat(ii,jj) = corr_mat(ii,jj) + alpha(ii,m)'*alpha(jj,m);
                %corr_mat(ii, jj+N) = corr_mat(ii, jj+N) + alpha(ii,m)'*beta(jj,m);
                %corr_mat(ii+N, jj) = corr_mat(ii+N, jj) + beta(ii,m)'*alpha(jj,m);
                %corr_mat(ii+N, jj+N) = corr_mat(ii+N, jj+N) + beta(ii,m)'*beta(jj,m);   
            end
        end
    end
    
    
                %corr_mat(ii,jj) = corr_mat(ii,jj) + alpha(ii,m)'*alpha(jj,m);
                %corr_mat(ii, jj+N) = corr_mat(ii, jj+N) + alpha(ii,m)'*beta(jj,m);
                %corr_mat(ii+N, jj) = corr_mat(ii+N, jj) + beta(ii,m)'*alpha(jj,m);
                %corr_mat(ii+N, jj+N) = corr_mat(ii+N, jj+N) + beta(ii,m)'*beta(jj,m);      
    %alpha_2 = U_2(1:N, N+1:end);
    %beta_2  = U_2(N+1:end,N+1:end);    
    
    
   
    




end