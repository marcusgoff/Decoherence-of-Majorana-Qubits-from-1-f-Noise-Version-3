%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: Load_Constants.m
% Purpose: Load constants and parameters into the MATLAB workspace.
%
%
% Author: Marcus C. Goffage
% Date: 30-Jan-2025
% Affiliation: University of New South Wales
%
% Paper: "Decoherence in Majorana Qubits by 1/f Noise"
% Paper Authors: A. Alase^1, M. C. Goffage^2, M. C. Cassidy^2, 
%                S. N. Coppersmith^{2*}
% Affiliations:  ^1 University of Sydney
%                ^2 University of New South Wales
%                *  Corresponding Author
%
% -------------------------------------------------------------------------
% ABOUT THIS SCRIPT
% -------------------------------------------------------------------------
% Description:
%   - Loads various constants into the MATLAB workspace.
%   - Many of these parameters are listed in Table 1 of the supplementary
%     information of the paper.
%
% Requirements:
%   - MATLAB R2024 or newer
%
% Output:
%   - Variables loaded into the workspace for use in subsequent scripts.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Units
hbar_mueVs = 6.582119*10^(-10); %mu eV-s 
k_f_inverse = 40e-9; %40e-9; % m 
a = (pi/2)*k_f_inverse; % m 
delta_muev = 110; %mu eV - Superconducting gap 
S_0 = 1; %mu eV - squared value - on the DOT
mu_muev = sqrt(8)*S_0/5; %mu eV - factor of 5 is due to ratio of quantum dot
                         % capacitance to nanowire capacitance. 
                      
t_conversion = 10^(-6)*delta_muev/hbar_mueVs; % unitless time/(1 \mu s); 
