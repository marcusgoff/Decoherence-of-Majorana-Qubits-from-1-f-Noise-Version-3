%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: Run_AddFig.m
% Purpose: Generate Additional Figure, the data from which is used in S1
%          "Decoherence in Majorana Qubits by 1/f Noise"
%
% Author: Marcus C. Goffage
% Date: 23-Apr-2025
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
%   - Generate Additional Figure, the data from which is used in S1
%
% Requirements:
%   - MATLAB R2024 or newer
%   - Dependencies: functions in /QPP_Library directory. 
%
% Output:
%   - Saves Figure_add and Figure_add_inset in the current directory
%   - Saves entire worskpace in the current directory
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

fprintf('\n Running Run_AddFig \n');

%% Load Bespoke Quasiparticle Poisoning Library 
addpath('../QPP_Library_submit')
addpath('..')

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% Load Constants from Table 1 in the Supplementary Information
run load_constants.m

warning('off', 'all');
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

gamma_min_ghz = 5; % GHz 20
gamma_max_ghz = 5e2; % GHz 2e3
num_gamma_points = 11; 

% Simulation parameters for each gamma point:
num_trials = 50; % number of trials at each gamma points
t_final_ns = 1; % nanoseconds 
t_final = t_final_ns*t_conversion*1e-3;

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
%gamma_sim_to_ghz_conversion = (0.5*0.0598)/10; % 1/GHz
gamma_sim_to_ghz_conversion = 0.5*(1/(t_conversion*10^(-3))); %1/GHz
gamma_min = gamma_min_ghz*gamma_sim_to_ghz_conversion;
gamma_max = gamma_max_ghz*gamma_sim_to_ghz_conversion;

% Number of Lattice Sites
N_a = round(chain_length_a./a); % divide by lattice constant, a, and round.
N_b = round(chain_length_b./a); 
N_c = round(chain_length_c./a); 

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% Run Simulations

% Chain Length A
fprintf('Running Simulation for Chain Length: %.0f microns \n', chain_length_a*1e6);
[L_even_mat_a, L_odd_mat_a, X_exp_mat_a, Z_exp_mat_a, ...
          L_even_cell_a, L_even_final_mat_a, ...
          L_even_final_vec_a, L_odd_final_vec_a, gamma_vec_a] = ...
        run_gamma_sweep(mu_low, mu_high, mu_offset, w, delta, N_a, BC, ...
                           num_gamma_points, gamma_min, gamma_max, ...
                           t_final, num_trials, init_state);
P_qpp_final_a = L_even_final_vec_a/2; %/2 to get P_qpp per chain in tetron
P_qpp_final_mat_a = L_even_final_mat_a/2;
fprintf('\n');
save('figure_add_fig_just_3_micron');

% Chain Length B
fprintf('Running Simulation for Chain Length: %.0f microns \n', chain_length_b*1e6);
[L_even_mat_b, L_odd_mat_b, X_exp_mat_b, Z_exp_mat_b, ...
          L_even_cell_b, L_even_final_mat_b, ...
          L_even_final_vec_b, L_odd_final_vec_b, gamma_vec_b] = ...
        run_gamma_sweep(mu_low, mu_high, mu_offset, w, delta, N_b, BC, ...
                           num_gamma_points, gamma_min, gamma_max, ...
                           t_final, num_trials, init_state);
P_qpp_final_b = L_even_final_vec_b/2; %/2 to get P_qpp per chain in tetron
P_qpp_final_mat_b = L_even_final_mat_b/2; 
fprintf('\n');

% Chain Length C
fprintf('Running Simulation for Chain Length: %.0f microns \n', chain_length_c*1e6);
[L_even_mat_c, L_odd_mat_c, X_exp_mat_c, Z_exp_mat_c, ...
          L_even_cell_c, L_even_final_mat_c, ...
          L_even_final_vec_c, L_odd_final_vec_c, gamma_vec_c] = ...
        run_gamma_sweep(mu_low, mu_high, mu_offset, w, delta, N_c, BC, ...
                           num_gamma_points, gamma_min, gamma_max, ...
                           t_final, num_trials, init_state);
P_qpp_final_c = L_even_final_vec_c/2;  %/2 to get P_qpp per chain in tetron
P_qpp_final_mat_c = L_even_final_mat_c/2;
fprintf('\n');

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%%
% Convert dimensionless gamma_vec (1/2 switching rate in normalized units) to
% gamma_vec_ghz (switching rate in GHz)
gamma_vec_ghz = 2*10^(-9)*gamma_vec_a*delta_muev/hbar_mueVs;

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% Get error Bars
[~, err_bars_a] = get_gamma_sweep_error_bars_func(P_qpp_final_mat_a, group_size);
[~, err_bars_b] = get_gamma_sweep_error_bars_func(P_qpp_final_mat_b, group_size);
[~, err_bars_c] = get_gamma_sweep_error_bars_func(P_qpp_final_mat_c, group_size);

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% Plot Figure Add Fig

% === Plot styling ===
marker_size = 16;
line_width = 4;
fig_ax_font_size = 30;
fig_border_LW = 2.3;
RGB = orderedcolors("gem");

col_1 = RGB(1,:);  
col_2 = RGB(2,:); 
col_3 = RGB(3,:);  

% === Create figure ===
currFig = figure();
currFig.Position = [205   341   823   505];
ax = gca;

% === Plot Curves ===
h_a = errorbar(gamma_vec_ghz, P_qpp_final_a, err_bars_a, '-.o', ...
   'Color', col_3, 'MarkerSize', marker_size, 'LineWidth', line_width*0.7, ...
   'MarkerFaceColor', col_3, 'MarkerEdgeColor', 'w'); hold on;
h_b = errorbar(gamma_vec_ghz, P_qpp_final_b, err_bars_b, '-.o', ...
   'Color', col_2, 'MarkerSize', marker_size, 'LineWidth', line_width*0.7, ...
   'MarkerFaceColor', col_2, 'MarkerEdgeColor', 'w'); hold on;
h_c = errorbar(gamma_vec_ghz, P_qpp_final_c, err_bars_c, '-.o', ...
   'Color', col_1, 'MarkerSize', marker_size, 'LineWidth', line_width*0.7, ...
   'MarkerFaceColor', col_1, 'MarkerEdgeColor', 'w'); hold on;

% === Log scale axes ===
set(ax, 'XScale', 'log');


% === Axis settings ===
ax.FontSize = fig_ax_font_size;
ax.LineWidth = fig_border_LW;
ax.TickDir = 'in';
ax.Box = 'off';
ax.XMinorTick = 'on';
ax.YMinorTick = 'off';

xlim([gamma_min_ghz*0.9, gamma_max_ghz*1.1]);


% ==== Labels ====
xlabel('Two level fluctuator switching rate, $\Gamma$ (GHz)', 'Interpreter', 'latex');
ylabel('P_{QPP, 1ns}');


hLeg = legend([h_a, h_b, h_c], ...
    {sprintf('$\\mathcal{L} = %i ~\\mu\\rm{m}$', 1e6*chain_length_a), ...
     sprintf('$\\mathcal{L} = %i ~\\mu\\rm{m}$', 1e6*chain_length_b), ...
     sprintf('$\\mathcal{L} = %i ~\\mu\\rm{m}$', 1e6*chain_length_c)}, ...
    'Interpreter', 'latex', ...
    'Location', 'northwest', ...
    'Box', 'off');

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% Save Figure Add Fig

% Save as MATLAB .fig image
saveas(gcf, 'figure_add_fig.fig');

% Save as PNG image
saveas(gcf, 'figure_add_fig.png');
saveas(gcf, 'figure_add_fig.svg');

%-------------------------------------------------------------------------
%% -----------------------------------------------------------------------
% Plot Figure Add Fig Inset
% ------------------------------------------------------------------------
% Control Parameters
gamma_indices_plot = [1, 6]; % 1 x 2 matrix
0
%% 

% Extract P_QPP versus L Data
length_vec = [chain_length_a, chain_length_b, chain_length_c]*1e6;
P_vers_L_data = [P_qpp_final_a(gamma_indices_plot), ...
    P_qpp_final_b(gamma_indices_plot), P_qpp_final_c(gamma_indices_plot)];

error_bars_vers_L = [err_bars_a(gamma_indices_plot), ...
    err_bars_b(gamma_indices_plot), err_bars_c(gamma_indices_plot)];

% Get Linear Fits 
% Fit to y = mx
length_vec_ext = linspace(0, 15, 20); 

warning('off', 'all'); % supress warning of raning start point
ft_a = fittype('m * x', 'independent', 'x', 'coefficients', 'm');
fit_a = fit(length_vec(:), P_vers_L_data(1,:).', ft_a);
fit_y_a = fit_a.m*length_vec_ext;

ft_b = fittype('m * x', 'independent', 'x', 'coefficients', 'm');
fit_b = fit(length_vec(:), P_vers_L_data(2,:).', ft_b);
fit_y_b = fit_b.m*length_vec_ext;

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% Plot Figure Add Fig Inset

% === Plot styling === %
marker_size_inset = 17*1.4;
line_width_inset = 4;
fig_ax_font_size_inset = 30*1.6;
col_1_light = col_1 + (1 - col_1) * 0.3;
col_2_light = col_2 + (1 - col_2) * 0.3;
col_3_light = col_3 + (1 - col_3) * 0.3;

% === Create figure ===
currFig = figure();
currFig.Position = [205 352 746 409];
ax = gca;
hold on;


% === Plot Linear Fits === %
h_a_fit = plot(length_vec_ext, fit_y_a, ...
    'Color', col_1_light, 'MarkerSize', marker_size_inset, 'LineWidth', line_width_inset, ...
    'MarkerFaceColor', col_1_light, 'MarkerEdgeColor', 'w');

h_b_fit = plot(length_vec_ext, fit_y_b, ...
    'Color', col_2_light, 'MarkerSize', marker_size_inset, 'LineWidth', line_width_inset, ...
    'MarkerFaceColor', col_2_light, 'MarkerEdgeColor', 'w');

% === Plot error bars === %
h_a = errorbar(length_vec, P_vers_L_data(1,:), error_bars_vers_L(1,:), 'o', ...
    'Color', col_1, 'MarkerSize', marker_size_inset, 'LineWidth', line_width_inset, ...
    'MarkerFaceColor', col_1, 'MarkerEdgeColor', 'w');

hold on;
h_b = errorbar(length_vec, P_vers_L_data(2,:), error_bars_vers_L(2,:), 'o', ...
    'Color', col_2, 'MarkerSize', marker_size_inset, 'LineWidth', line_width_inset, ...
    'MarkerFaceColor', col_2, 'MarkerEdgeColor', 'w');

% === Axis settings ===
ax.FontSize = fig_ax_font_size_inset;
ax.LineWidth = 2.3;
ax.TickDir = 'in';
ax.Box = 'off';
ax.XMinorTick = 'off';
axis tight
%ax.YMinorTick = 'off';
ax.XTick = [0, 5, 10];

xlim([0, 11]);
ylim([0, 0.06]);
xlabel('$\mathcal{L} (\mu \rm m)$', 'interpreter', 'latex');
ylabel('P_{QPP, 5ns}');

% % Set major Y ticks, excluding 0.05 and 0.15
ax.YTick = [0, 0.05];
% Enable minor ticks
ax.YMinorTick = 'on';
% Manually specify the positions of minor ticks
ax.YAxis.MinorTickValues = [0.025];

% === Manually draw box === %
xl = xlim(ax); yl = ylim(ax);
hold on;
plot([xl(1) xl(2)], [yl(1) yl(1)], 'k', 'LineWidth', ax.LineWidth) % bottom
plot([xl(1) xl(2)], [yl(2) yl(2)], 'k', 'LineWidth', ax.LineWidth) % top
plot([xl(1) xl(1)], [yl(1) yl(2)], 'k', 'LineWidth', ax.LineWidth) % left
plot([xl(2) xl(2)], [yl(1) yl(2)], 'k', 'LineWidth', ax.LineWidth) % right
hold off;

text(5.2, 0.05, ...
    sprintf('$\\Gamma = %.0f\\,\\mathrm{GHz}$', gamma_vec_ghz(gamma_indices_plot(1))), ...
    'Color', col_2, 'Interpreter', 'latex', ...
    'FontSize', 0.9 * fig_ax_font_size_inset, ...
    'HorizontalAlignment', 'center');

text(7.5, 0.01, ...
    sprintf('$\\Gamma = %.0f\\,\\mathrm{GHz}$', gamma_vec_ghz(gamma_indices_plot(2))), ...
    'Color', col_1, 'Interpreter', 'latex', ...
    'FontSize', 0.8 * fig_ax_font_size_inset, ...
    'HorizontalAlignment', 'center');

hold off;

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% Plot Figure Add Fig Inset

% === Plot styling ===
marker_size_inset = 17 * 1.4;
line_width_inset       = 4;
fig_ax_font_size_inset = 30 * 1.6;

col_1_light = col_1 + (1 - col_1) * 0.3;
col_2_light = col_2 + (1 - col_2) * 0.3;
col_3_light = col_3 + (1 - col_3) * 0.3;

% === Create figure ===
currFig          = figure();
currFig.Position = [205, 352, 746, 409];
ax               = gca;
hold on;

% === Plot linear fits ===
plot(length_vec_ext, fit_y_a, 'Color', col_1_light, ...
    'MarkerSize', marker_size_inset, 'LineWidth', line_width_inset, ...
    'MarkerFaceColor', col_1_light, 'MarkerEdgeColor', 'w');

plot(length_vec_ext, fit_y_b, 'Color', col_2_light, ...
    'MarkerSize', marker_size_inset, 'LineWidth', line_width_inset, ...
    'MarkerFaceColor', col_2_light, 'MarkerEdgeColor', 'w');

% === Plot error bars ===
errorbar(length_vec, P_vers_L_data(1,:), error_bars_vers_L(1,:), 'o', ...
    'Color', col_1, 'MarkerSize', marker_size_inset, 'LineWidth', line_width_inset, ...
    'MarkerFaceColor', col_1, 'MarkerEdgeColor', 'w');

errorbar(length_vec, P_vers_L_data(2,:), error_bars_vers_L(2,:), 'o', ...
    'Color', col_2, 'MarkerSize', marker_size_inset, 'LineWidth', line_width_inset, ...
    'MarkerFaceColor', col_2, 'MarkerEdgeColor', 'w');

% === Axis settings ===
ax.FontSize      = fig_ax_font_size_inset;
ax.LineWidth     = 2.3;
ax.TickDir       = 'in';
ax.Box           = 'off';
ax.XMinorTick    = 'off';
ax.XTick         = [0, 5, 10];
ax.YTick         = [0, 0.05];
ax.YMinorTick    = 'on';
ax.YAxis.MinorTickValues = [0.01, 0.02, 0.03, 0.04, 0.06];

xlim([0, 11]);
ylim([0, 0.055]);
xlabel('$\mathcal{L} (\mu \mathrm{m})$', 'Interpreter', 'latex');
ylabel('P_{QPP, 5ns}');

% === Manually draw box ===
xl = xlim(ax);
yl = ylim(ax);
plot([xl(1) xl(2)], [yl(1) yl(1)], 'k', 'LineWidth', ax.LineWidth); % bottom
plot([xl(1) xl(2)], [yl(2) yl(2)], 'k', 'LineWidth', ax.LineWidth); % top
plot([xl(1) xl(1)], [yl(1) yl(2)], 'k', 'LineWidth', ax.LineWidth); % left
plot([xl(2) xl(2)], [yl(1) yl(2)], 'k', 'LineWidth', ax.LineWidth); % right

% === Add annotations ===
text(8.2, 0.014, ...
    sprintf('$\\Gamma = %.0f\\,\\mathrm{GHz}$', gamma_vec_ghz(gamma_indices_plot(1))), ...
    'Color', col_1, 'Interpreter', 'latex', ...
    'FontSize', 0.8 * fig_ax_font_size_inset, ...
    'HorizontalAlignment', 'center');

text(6.5, 0.04, ...
    sprintf('$\\Gamma = %.0f\\,\\mathrm{GHz}$', gamma_vec_ghz(gamma_indices_plot(2))), ...
    'Color', col_2, 'Interpreter', 'latex', ...
    'FontSize', 0.8 * fig_ax_font_size_inset, ...
    'HorizontalAlignment', 'center');

hold off;

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% Save Figure Add Fig

% Save as MATLAB .fig image
saveas(gcf, '../results/figure_add_fig_inset.fig');

% Save as PNG and SVG image
saveas(gcf, '../results/figure_add_fig_inset.png');
saveas(gcf, '../results/figure_add_fig_inset.svg');


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% Save Data
%close all;
save('figure_add_fig_data');

%% === Export plotted data for Figure Add Fig to CSV ===

% Ensure output directory exists
output_dir = "../results";

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

filename = fullfile(output_dir, "Figure_add_fig.csv");

% === Build data matrix ===

N_gamma = numel(gamma_vec_ghz);

data_mat = zeros(N_gamma, 7);

data_mat(:,1) = gamma_vec_ghz(:);

data_mat(:,2) = P_qpp_final_a(:);
data_mat(:,3) = err_bars_a(:);

data_mat(:,4) = P_qpp_final_b(:);
data_mat(:,5) = err_bars_b(:);

data_mat(:,6) = P_qpp_final_c(:);
data_mat(:,7) = err_bars_c(:);

% === Column headers ===

headers = { ...
    'Gamma_GHz', ...
    sprintf('P_QPP_L_%i_um', round(1e6*chain_length_a)), ...
    sprintf('Err_L_%i_um', round(1e6*chain_length_a)), ...
    sprintf('P_QPP_L_%i_um', round(1e6*chain_length_b)), ...
    sprintf('Err_L_%i_um', round(1e6*chain_length_b)), ...
    sprintf('P_QPP_L_%i_um', round(1e6*chain_length_c)), ...
    sprintf('Err_L_%i_um', round(1e6*chain_length_c)) ...
    };

% === Write table ===

T = array2table(data_mat, 'VariableNames', headers);

writetable(T, filename);

disp("CSV file written to:");
disp(filename);