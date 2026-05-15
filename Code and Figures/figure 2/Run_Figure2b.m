%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: Run_Figure2b.m
% Purpose: Generates Data for Figure 2b for the Paper
%          "Decoherence in Majorana Qubits by 1/f Noise"
%
% Author: Marcus C. Goffage
% Date: 23-Apr-2026
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
%   - Generates Figure 2b from our paper.
%
% Requirements:
%   - MATLAB R2024 or newer
%   - Dependencies: functions in /QPP_Library directory. 
%
% Output:
%   - Saves entire worskpace in the current directory
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

%%
fprintf('\n Running Run_FigureS1 \n');

%% Load Bespoke Quasiparticle Poisoning Library 
addpath('../QPP_Library_submit')
addpath('../')

warning('off', 'all');
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% Load Constants from Table 1 in the Supplementary Information
run load_constants.m
delta_0_muev = delta_muev;

% Run at the higher value of delta_mu:
%mu_muev = mu_muev*5; 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% Declare Simulation Input Parameters

% Take from run_script_delta_search.m and go from there. 

% Kitaev Chain Parameters. Note that our QPP_Library we use the Kitaev
% convention employed in our previous paper arXiv:2504.17485. This requires
% multiplying delta and w (as used in our latest paper) 
% by -1/2 and 1/2 respectively, as done below. 
delta_0 = 1*(-1/2); % dimensionless - as used by QPP_lib
w = 3.189*(1/2);  % dimensionless - as used by QPP_lib
BC = "OBC";
delta_0_GHz = 10^(-9)*delta_0_muev/(2*pi*hbar_mueVs);

init_state = '+'; % Initialise qubit in |+> state

% TLF Parameters
mu_low = 0; % dimensionless
mu_high = mu_muev./delta_muev; % dimensinless
mu_offset = 0; % dimensionless

gamma_min_ghz = 1e-1; % GHz 
gamma_max_ghz = 5e2;  % GHz  % 5e2
num_gamma_points = 20; % previously 36

% Simulation parameters for each gamma point:
num_trials = 50;  % number of trials at each gamma points - previously 50
t_final_ns = 30; % nanoseconds
t_final = t_final_ns*t_conversion*1e-3;

% Error Bar Parameters
% Number of trials in each group. Error bar lengths are taken as the
% standard error of the means of each group. 
group_size = 5; 

% Choose chain length
chain_length = 3e-6; % metres
N = round(chain_length/a);

% Choose delta-factors
delta_factors = [0.1, 1/3, 1, 3, 9];

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% Convert Gamma and Chain length to Dimensionless Values
% Gamma-Sweep Variables
%gamma_sim_to_ghz_conversion = (0.5*0.0598)/10; % 1/GHz
gamma_sim_to_ghz_conversion = 0.5*(1/(t_conversion*10^(-3))); %1/GHz
gamma_min = gamma_min_ghz*gamma_sim_to_ghz_conversion;
gamma_max = gamma_max_ghz*gamma_sim_to_ghz_conversion;


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% Run Simulations

% Preallocate storage
P_qpp_final      = cell(1, numel(delta_factors));
P_qpp_final_mat  = cell(1, numel(delta_factors));

for i = 1:length(delta_factors)
    fprintf('Running Simulation for Delta: %.0f GHz \n', delta_factors(i) * delta_0_GHz);

    [~, ~, ~, ~, ...
              ~, L_even_final_mat, ...
              L_even_final_vec, ~, gamma_vec] = ...
            run_gamma_sweep(mu_low, mu_high, mu_offset, w, ...
                            delta_factors(i)*delta_0, N, BC, ...
                            num_gamma_points, gamma_min, gamma_max, ...
                            t_final, num_trials, init_state);

    % Store results
    P_qpp_final{i}     = L_even_final_vec / 2; % /2 to get P_qpp per chain in tetron
    P_qpp_final_mat{i} = L_even_final_mat/2;

    fprintf('\n');
end


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%%
% Convert dimensionless gamma_vec (1/2 switching rate in normalized units) to
% gamma_vec_ghz (switching rate in GHz)
gamma_vec_ghz = 2*10^(-9)*gamma_vec*delta_muev/hbar_mueVs;

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% Get error Bars
err_bars = cell(1, numel(P_qpp_final_mat));

for i = 1:numel(P_qpp_final_mat)
    [~, err_bars{i}] = get_gamma_sweep_error_bars_func(P_qpp_final_mat{i}, group_size);
end


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% Plot Figure S1

% === Plot styling ===
marker_size = 16;
line_width = 4;
fig_ax_font_size = 30;
fig_border_LW = 2.3;
RGB = orderedcolors("gem");


col_mat = RGB([4, 3, 2, 1, 5], :);
blend_factor = 0.3;
col_mat_light = col_mat + (1 - col_mat) * blend_factor;

% === Create figure ===
currFig = figure();
currFig.Position = [205   341   823   505];
ax = gca;

% === Plot Curves ===
h = cell(1, numel(P_qpp_final));

t_final_s = t_final_ns*1e-9;
% Plot R_qpp = P_qpp/t_final_s *10^(-6) [MHz]
for i = 1:numel(P_qpp_final)
    y_plot = P_qpp_final{i}*1e-6/(t_final_s);
    h{i} = errorbar(gamma_vec_ghz, ...
                    y_plot/max(y_plot), ...
                    err_bars{i}/max(y_plot), ...
                    'o', ...
                    'Color',           col_mat(i, :), ...
                    'MarkerSize',      marker_size, ...
                    'LineWidth',       line_width * 0.7, ...
                    'MarkerFaceColor', col_mat(i, :), ...
                    'MarkerEdgeColor', 'w');
    hold on;
end

% === Log scale axes ===
set(ax, 'XScale', 'log');
set(ax, 'YScale', 'log');


% === Axis settings ===
ax.FontSize = fig_ax_font_size;
ax.LineWidth = fig_border_LW;
ax.TickDir = 'in';
ax.Box = 'off';
ax.XMinorTick = 'off';
ax.YMinorTick = 'off';

xlim([gamma_min_ghz*0.9, gamma_max_ghz*1.1]);
%ylim([0.08, 1.2]);
%xlim([0.1,1000]);

% ==== Labels ====
xlabel('Two level fluctuator switching rate, $\Gamma$ (GHz)', 'Interpreter', 'latex');
%ylabel('P_{QPP, 5ns}');
%ylabel(sprintf('P_{QPP, %ins}', t_final_ns));
ylabel('R_{QPP}/R_{QPP}^{max}');


% Build legend labels dynamically
legend_labels = cell(1, numel(delta_factors));
for i = 1:numel(delta_factors)
    legend_labels{i} = sprintf('$\\Delta = %.2f ~ \\mu \\rm{eV}$', ...
                               delta_factors(i) * delta_muev);
end

% Create legend from all handles
hLeg = legend([h{:}], legend_labels, ...
    'Interpreter', 'latex', ...
    'Location', 'south', ...
    'Box', 'off');


%%  Plot Fermi Golden Rule
hold on;
delta_vec = delta_factors*delta_0_muev;


delta_mu = (mu_high - mu_low);
delta_mu_muev = delta_mu*delta_0_muev;
gamma_vec_hz = gamma_vec_ghz*(1e9);

fermi_velocity = (2*w*delta_muev)*a/hbar_mueVs; % w*a/hbar
% Factor of 2 is simply because "w" here is defined with factor of 1/2
% compared to convention in paper (code uses conventions of our earlier
% PRB paper). 

for ii = 1:length(delta_vec)
    delta_curr_muev = delta_vec(ii);
    
    % fix fermi velocity
    fgr = chain_length*delta_mu_muev.^2./(16*hbar_mueVs*fermi_velocity*delta_curr_muev) ...
        *gamma_vec_hz.*(1 + hbar_mueVs^2*gamma_vec_hz.^2/delta_curr_muev.^2).^(-3/2);

    y_plot = fgr*1e-6;
    fgr_no_fact = gamma_vec_hz.*(1 + hbar_mueVs^2.*gamma_vec_hz.^2./delta_curr_muev).^(-3/2);
    plot(gamma_vec_ghz, y_plot/max(y_plot), 'color', col_mat(ii, :), 'LineWidth', 3, 'LineStyle', '--'); 
    % For P_qpp, plot: fgr*t_final_ns*1e-9
end

% Create legend from all handles
hLeg = legend([h{:}], legend_labels, ...
    'Interpreter', 'latex', ...
    'Location', 'southeast', ...
    'Box', 'off');

fermi_velocity = (2*w*delta_muev)*a/hbar_mueVs; % w*a/hbar
R_max = (chain_length.*delta_mu_muev.^2)./(8*3^(3/2)*hbar_mueVs.^2*fermi_velocity);
P_qpp_max = R_max*t_final_ns*1e-9;



hLeg = legend([h{:}], legend_labels, ...
    'Interpreter', 'latex', ...
    'Location', 'south', ...
    'Box', 'off');

ax = gca;

ylims = ax.YLim;
ax.YLim = [ylims(1)*0.7, ylims(2)];
ax.XTick = 10.^(-10:10);

ax.YTick = 10.^(-10:10);

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% Save Figure 

% Save as MATLAB .fig image
saveas(gcf, 'figure_2b_data.fig');

% Save as PNG image
saveas(gcf, 'figure_2b_data.png');
saveas(gcf, 'figure_2b_data.svg');



% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% Save Data
close all;
save('figure_2b_data');

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% === Export plotted data to CSV ===

% Number of gamma points and curves
N_gamma = numel(gamma_vec_ghz);
N_curves = numel(P_qpp_final);

% Preallocate data matrix
% Columns:
% 1: Gamma
% then for each curve:
%   R_QPP
%   Error
% then FGR curves
N_cols = 1 + 2*N_curves + N_curves;

data_mat = zeros(N_gamma, N_cols);

% ---- Column headers ----
headers = strings(1, N_cols);

col = 1;

headers(col) = "Gamma_GHz";
data_mat(:, col) = gamma_vec_ghz(:);
col = col + 1;

t_final_s = t_final_ns * 1e-9;

% === Add simulation data ===

for i = 1:N_curves

    % R_QPP
    R_qpp = real(P_qpp_final{i}) * 1e-6 / t_final_s;

    headers(col) = sprintf("R_QPP_Delta_%.1f_GHz", ...
        delta_factors(i) * delta_0_GHz);

    data_mat(:, col) = R_qpp(:);
    col = col + 1;

    % Error bars
    headers(col) = sprintf("Err_Delta_%.1f_GHz", ...
        delta_factors(i) * delta_0_GHz);

    data_mat(:, col) = err_bars{i}(:);
    col = col + 1;

end

% === Add Fermi Golden Rule curves ===

delta_vec = delta_factors * delta_0_muev;

delta_mu = (mu_high - mu_low);
delta_mu_muev = delta_mu * delta_0_muev;

gamma_vec_hz = gamma_vec_ghz * 1e9;

fermi_velocity = (2*w*delta_muev)*a/hbar_mueVs;

for ii = 1:N_curves

    delta_curr_muev = delta_vec(ii);

    fgr = chain_length * delta_mu_muev.^2 ./ ...
        (16 * hbar_mueVs * fermi_velocity * delta_curr_muev) ...
        .* gamma_vec_hz ...
        .* (1 + hbar_mueVs^2 * gamma_vec_hz.^2 / delta_curr_muev.^2).^(-3/2);

    headers(col) = sprintf("FGR_Delta_%.1f_GHz", ...
        delta_factors(ii) * delta_0_GHz);

    data_mat(:, col) = fgr(:) * 1e-6;

    col = col + 1;

end

% === Write to CSV in Current Directory ===

output_dir = ".";

% Create directory if it doesn't exist (recommended)
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

filename = fullfile(output_dir, "figure_2b_data.csv");

T = array2table(data_mat, 'VariableNames', headers);

writetable(T, filename);

disp("CSV file written to:");
disp(filename);

