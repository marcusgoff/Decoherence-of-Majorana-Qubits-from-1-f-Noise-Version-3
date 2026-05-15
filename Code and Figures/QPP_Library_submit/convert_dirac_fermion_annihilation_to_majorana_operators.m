%---------------------------------------------------%
% Convert Dirac Fermionic Annihilation Operators 
% to Majorana Operators
%---------------------------------------------------%
% Author: Marcus Goffage
% Date: 14.04.24
%
% Input: dirac_fermion_ann - size Nx1 or 1xN cell array containing N dirac 
%        fermion annihilation operators.
% Ouput: majorana_operators - size 2*Nx1 cell array containing majorana
%        fermion operators, in format: 
%       (x_1 ... x_N, p_1 ... p_N)^T, where
%        x_j = (a_j + a_j^dag)/sqrt(2), 
%        p_j = (a_j - a_j^dag)/(1i*sqrt(2),
%        and, a_j are the annihilation operators in dirac_fermion_ann.
%

function majorana_operators =...
    convert_dirac_fermion_annihilation_to_majorana_operators(dirac_fermion_ann)
    if min(size(dirac_fermion_ann)) ~= 1
       error('Input to convert_dirac_fermion_annihilation_to_majorana_operators must be Nx1 or 1xN cell array'); 
    end
    N = length(dirac_fermion_ann);
    majorana_operators = cell(2*N,1);
    
    for ii = 1:N
        x_temp = 1/sqrt(2)*(dirac_fermion_ann{ii} + dirac_fermion_ann{ii}');
        p_temp = -1i/sqrt(2)*(dirac_fermion_ann{ii} - dirac_fermion_ann{ii}');
        majorana_operators{ii} = x_temp;
        majorana_operators{ii + N} = p_temp;
    end
end
    
    
    
    