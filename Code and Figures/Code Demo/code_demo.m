%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: code_demo.m
% Purpose: Demo of the code used in our paper
%          "Decoherence in Majorana Qubits by 1/f Noise"
%
% Author: Marcus C. Goffage
% Date: 12-Aug-2025
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
%   - Demo code which runs in under 10 minutes. Calculates and plots P_qpp
%     as a function of time for L = 3 microns, and for a TLF rate of 200
%     GHz, for 3 nanosecond or untill P_{qpp} = 0.1 (whichever occurs first),
%     this corresponds to the yellow line in figure 2b in the main tex). 
%
% Requirements:
%   - MATLAB R2024 or newer
%   - Dependencies: functions in /QPP_Library directory. 
%
% Output:
%   - Saves results in current directory
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fprintf('Running Run_Figure2b \n');

%% Load Bespoke Quasiparticle Poisoning Library 
addpath('../QPP_Library_submit')
addpath('../')

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% Load Constants from Table 1 in the Supplementary Information
run load_constants.m

warning('off', 'all')
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% Declare Simulation Input Parameters

% Kitaev Chain Parameters. Note that our QPP_Library we use the Kitaev
% convention employed in our previous paper arXiv:2504.17485. This requires
% multiplying delta and w (as used in our latest paper) 
% by -1/2 and 1/2 respectively, as done below. 
delta = 1*(-1/2); % dimensionless - as used by QPP_lib
w = 3.189*(1/2);  % dimensionless - as used by QPP_lib
BC = "OBC";

init_state = '+'; % Initialise qubit in |+> state

% TLF Parameters
mu_low = 0; % dimensionless
mu_high = mu_muev./delta_muev; % dimensionless
mu_offset = 0; % dimensionless

gamma_1_ghz = 200;
gamma_2_ghz = 20;
num_gamma_points = 1; % do not change 

% Simulation parameters for each gamma point:
num_trials = 20; % number of trials at each gamma points

% Error Bar Parameters
% Number of trials in each group. Error bar lengths are taken as the
% standard error of the means of each group. 
group_size = 5; 

% Choose three chain lengths
chain_length_a = 3e-6; % metres
chain_length_b = 5e-6; % metres
chain_length_c = 10e-6; % metres


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% Convert Gamma and Chain length to Dimensionless Values
% Gamma-Sweep Variables
gamma_sim_to_ghz_conversion = (0.5*0.0598)/10; % 1/GHz
gamma_1 = gamma_1_ghz*gamma_sim_to_ghz_conversion;
gamma_2 = gamma_2_ghz*gamma_sim_to_ghz_conversion;

% Number of Lattice Sites
N_a = round(chain_length_a./a); % divide by lattice constant, a, and round.
% N_b = round(chain_length_b./a); 
% N_c = round(chain_length_c./a); 


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% Run Simulations


% - - - - - - - - - - - - - %
% Gamma 1 Simulations

% Chain Length A
fprintf('Running Simulation for Chain Length: %.0f microns \n TLF rate: %.0f\n', chain_length_a*1e6, gamma_1_ghz);
t_final_ns = 3; 
t_final = t_final_ns*t_conversion*1e-3;
[L_even_mat_a, L_odd_mat_a, X_exp_mat_a, Z_exp_mat_a, ...
          L_even_cell_a, L_even_final_mat_a, ...
          L_even_final_vec_a, L_odd_final_vec_a, gamma_vec_a, t_vec_cell_a] = ...
        run_gamma_sweep(mu_low, mu_high, mu_offset, w, delta, N_a, BC, ...
                           num_gamma_points, gamma_1, gamma_1, ...
                           t_final, num_trials, init_state);
P_qpp_a_gamma_1 = mean(L_even_cell_a{1}/2,2);
t_vec_a_gamma_1 = t_vec_cell_a{1}/t_conversion*1e3; % nanosec
fprintf('\n');
% 
% % Chain Length B
% fprintf('Running Simulation for Chain Length: %.0f microns \n TLF rate: %.0f\n', chain_length_b*1e6, gamma_1_ghz);
% t_final_ns = 1.5;
% t_final = t_final_ns*t_conversion*1e-3;
% [L_even_mat_b, L_odd_mat_b, X_exp_mat_b, Z_exp_mat_b, ...
%           L_even_cell_b, L_even_final_mat_b, ...
%           L_even_final_vec_b, L_odd_final_vec_b, gamma_vec_b, t_vec_cell_b] = ...
%         run_gamma_sweep(mu_low, mu_high, mu_offset, w, delta, N_b, BC, ...
%                            num_gamma_points, gamma_1, gamma_1, ...
%                            t_final, num_trials, init_state);
% P_qpp_b_gamma_1 = mean(L_even_cell_b{1}/2,2);
% t_vec_b_gamma_1 = t_vec_cell_b{1}/t_conversion*1e3; % nanosec
% fprintf('\n');
% 
% % Chain Length C
% fprintf('Running Simulation for Chain Length: %.0f microns \n TLF rate: %.0f\n', chain_length_c*1e6, gamma_1_ghz);
% t_final_ns = 1; 
% t_final = t_final_ns*t_conversion*1e-3;
% [L_even_mat_c, L_odd_mat_c, X_exp_mat_c, Z_exp_mat_c, ...
%           L_even_cell_c, L_even_final_mat_c, ...
%           L_even_final_vec_c, L_odd_final_vec_c, gamma_vec_c, t_vec_cell_c] = ...
%         run_gamma_sweep(mu_low, mu_high, mu_offset, w, delta, N_c, BC, ...
%                            num_gamma_points, gamma_1, gamma_1, ...
%                            t_final, num_trials, init_state);
% P_qpp_c_gamma_1 = mean(L_even_cell_c{1}/2,2);
% t_vec_c_gamma_1 = t_vec_cell_c{1}/t_conversion*1e3; % nanosec
% fprintf('\n');
% 
% % - - - - - - - - - - - - - %
% % Gamma 2 Simulations
% 
% % Chain Length A
% fprintf('Running Simulation for Chain Length: %.0f microns \n TLF rate: %.0f\n', chain_length_a*1e6, gamma_2_ghz);
% t_final_ns = 10; 
% t_final = t_final_ns*t_conversion*1e-3;
% [L_even_mat_a, L_odd_mat_a, X_exp_mat_a, Z_exp_mat_a, ...
%           L_even_cell_a, L_even_final_mat_a, ...
%           L_even_final_vec_a, L_odd_final_vec_a, gamma_vec_a, t_vec_cell_a] = ...
%         run_gamma_sweep(mu_low, mu_high, mu_offset, w, delta, N_a, BC, ...
%                            num_gamma_points, gamma_2, gamma_2, ...
%                            t_final, num_trials, init_state);
% P_qpp_a_gamma_2 = mean(L_even_cell_a{1}/2,2);
% t_vec_a_gamma_2 = t_vec_cell_a{1}/t_conversion*1e3; % nanosec
% fprintf('\n');
% 
% % Chain Length B
% fprintf('Running Simulation for Chain Length: %.0f microns \n TLF rate: %.0f\n', chain_length_b*1e6, gamma_2_ghz);
% t_final_ns = 10;
% t_final = t_final_ns*t_conversion*1e-3;
% [L_even_mat_b, L_odd_mat_b, X_exp_mat_b, Z_exp_mat_b, ...
%           L_even_cell_b, L_even_final_mat_b, ...
%           L_even_final_vec_b, L_odd_final_vec_b, gamma_vec_b, t_vec_cell_b] = ...
%         run_gamma_sweep(mu_low, mu_high, mu_offset, w, delta, N_b, BC, ...
%                            num_gamma_points, gamma_2, gamma_2, ...
%                            t_final, num_trials, init_state);
% P_qpp_b_gamma_2 = mean(L_even_cell_b{1}/2,2);
% t_vec_b_gamma_2 = t_vec_cell_b{1}/t_conversion*1e3; % nanosec
% fprintf('\n');
% 
% % Chain Length C
% fprintf('Running Simulation for Chain Length: %.0f microns \n TLF rate: %.0f\n', chain_length_c*1e6, gamma_2_ghz);
% t_final_ns = 6; 
% t_final = t_final_ns*t_conversion*1e-3;
% [L_even_mat_c, L_odd_mat_c, X_exp_mat_c, Z_exp_mat_c, ...
%           L_even_cell_c, L_even_final_mat_c, ...
%           L_even_final_vec_c, L_odd_final_vec_c, gamma_vec_c, t_vec_cell_c] = ...
%         run_gamma_sweep(mu_low, mu_high, mu_offset, w, delta, N_c, BC, ...
%                            num_gamma_points, gamma_2, gamma_2, ...
%                            t_final, num_trials, init_state);
% P_qpp_c_gamma_2 = mean(L_even_cell_c{1}/2,2);
% t_vec_c_gamma_2 = t_vec_cell_c{1}/t_conversion*1e3; % nanosec
% fprintf('\n');


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% Plot  

% === Plot styling parameters ===
marker_size = 17;
line_width = 4;
fig_border_LW = 2.3;
fig_ax_font_size = 30;
RGB = orderedcolors("gem");

col_1 = RGB(1,:);  
col_2 = RGB(2,:);  
col_3 = RGB(3,:);  

% === Create figure ===
currFig = figure();
currFig.Position = [205 258 786 503];

ax = gca;
ax.FontSize = fig_ax_font_size;
ax.LineWidth = fig_border_LW;
hold on;

xlim([0, 10]);
ax.XTick = [0, 5, 10];

% === DUMMY frequency legend handles ===
h200_dummy = plot(nan, nan, '-', 'Color', 'k', 'LineWidth', line_width);
h20_dummy = plot(nan, nan, 's-', 'Color', 'k', 'LineWidth', line_width, ...
                  'MarkerFaceColor', 'k', 'MarkerSize', marker_size);

% === 200 GHz curves ===
% Divide L_even/2 to get P_qpp per chain 
% h10um_200 = plot(t_vec_c_gamma_1, P_qpp_c_gamma_1, ...
%     'Color', col_1, 'LineWidth', line_width);
% 
% h5um_200 = plot(t_vec_b_gamma_1, P_qpp_b_gamma_1, ...
%     'Color', col_2, 'LineWidth', line_width);

h3um_200 = plot(t_vec_a_gamma_1, P_qpp_a_gamma_1, ...
    'Color', col_3, 'LineWidth', line_width);

% === 20 GHz curves w ===
col_1_light = col_1 + (1 - col_1) * 0.3;
col_2_light = col_2 + (1 - col_2) * 0.3;
col_3_light = col_3 + (1 - col_3) * 0.3;


col_1_pastel = [0.3, 0.7, 0.8]; 
col_2_pastel = [0.8, 0.5, 0.7]; 
col_3_pastel = [0.9, 0.65, 0.4];

% h10um_20 = plot(t_vec_c_gamma_2, P_qpp_c_gamma_2, ...
%     '-', 'Color', col_1_pastel, 'LineWidth', line_width);
% 
% h5um_20 = plot(t_vec_b_gamma_2, P_qpp_b_gamma_2, ...
%     'Color', col_2_pastel, 'LineWidth', line_width);
% 
% h3um_20 = plot(t_vec_a_gamma_2, P_qpp_a_gamma_2, ...
%     'Color', col_3_pastel, 'LineWidth', line_width);


% === Labels and limits ===
xlabel('Time (ns)');
ylabel('P_{QPP}');
ylim([0, 0.1]);
yticks([0, 0.05, 0.1]);

% === Axis tweaks ===
ax.TickDir = 'in';
ax.Box = 'on';
ax.XMinorTick = 'on';
ax.YMinorTick = 'off';
ax.Box = 'off';
ax.YAxis.MinorTickValues = 1:10;

hold off;


% % === Manually draw box === %
% xl = xlim(ax); yl = ylim(ax);
% hold on;
% plot([xl(1) xl(2)], [yl(1) yl(1)], 'k', 'LineWidth', ax.LineWidth) % bottom
% plot([xl(1) xl(2)], [yl(2) yl(2)], 'k', 'LineWidth', ax.LineWidth) % top
% plot([xl(1) xl(1)], [yl(1) yl(2)], 'k', 'LineWidth', ax.LineWidth) % left
% plot([xl(2) xl(2)], [yl(1) yl(2)], 'k', 'LineWidth', ax.LineWidth) % right
% hold off;


% % === Legend ===
% currLeg = legend( ...
%     [h10um_200, h5um_200, h3um_200, ...
%      h10um_20,  h5um_20,  h3um_20], ...
%     { ...
%     '$\mathcal{L} = 10\,\mu\rm{m},\;\Gamma = 200~\rm{GHz}$', ...
%     '$\mathcal{L} = 5\,\mu\rm{m},\;\;\;\Gamma = 200~\rm{GHz}$', ...
%     '$\mathcal{L} = 3\,\mu\rm{m},\;\;\;\Gamma = 200~\rm{GHz}$', ...
%     '$\mathcal{L} = 10\,\mu\rm{m},\;\Gamma = 20~\rm{GHz}$', ...
%     '$\mathcal{L} = 5\,\mu\rm{m},\;\;\;\Gamma = 20~\rm{GHz}$', ...
%     '$\mathcal{L} = 3\,\mu\rm{m},\;\;\;\Gamma = 20~\rm{GHz}$' ...
%     }, ...
%     'Interpreter', 'latex', ...
%     'NumColumns', 1, ...
%     'Location', 'southeast');

% === Legend ===
currLeg = legend( ...
     h3um_200, ...
    { ...
    '$\mathcal{L} = 3\,\mu\rm{m},\;\;\;\Gamma = 200~\rm{GHz}$', ...
    }, ...
    'Interpreter', 'latex', ...
    'NumColumns', 1, ...
    'Location', 'southeast');

currLeg.Box = 'off';
currLeg.ItemTokenSize = [30, 10];
%currLeg.Position = [0.5258 0.1639 0.3791 0.3994];
currLeg.FontSize = 20;
% First create the ylabel
ylab = ylabel('P_{QPP}');


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% Save Figure 2b

% Save as MATLAB .fig image
saveas(gcf, 'demo_figure.fig');

% Save as PNG and SVG image
saveas(gcf, 'demo_figure.png');
saveas(gcf, 'demo_figure.svg');


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% Save Data
%close all;
save('demo_figure');
