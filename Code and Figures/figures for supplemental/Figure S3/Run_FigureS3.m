%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: Run_FigureS3.m
% Purpose: Generate Figure S3 from the paper 
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
%   - Generates Figure S3 from our paper.
%
% Requirements:
%   - MATLAB R2024 or newer
%   - Dependencies: functions in /QPP_Library directory. 
%
% Output:
%   - Saves Figure_S3 in ../results
%   - Saves entire worskpace in the current directory
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
%%
fprintf('\n Running Run_FigureS3 \n');

%% Load Bespoke Quasiparticle Poisoning Library 
addpath('../../QPP_Library_submit')
addpath('../..')

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% Load Constants from Table 1 in the Supplementary Information
run load_constants.m

warning('off', 'all');
%%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

delta = -0.5; %(0.033/2); %0.033
%w = 0; %0.5*3.32*10^3/delta_muev; %1
%N = 122; %122;%253;    
BC = "OBC";
N = 159; %159, 48

init_state = '+';

mu = 0; % also the low value (tetron initialised at H(mu_init))

mu_offset = 0;
 
num_w_points = 18;
% Note: w_vec is unitless 
w_vec = linspace(150, 500, num_w_points)/(2*delta_muev)*delta_muev/(110); %scaling factor to get in regime which worked for delta = 110 mueV.
%w_paper = 350.8/(2*delta_muev);
w_special = [197.52, 350.8, 446.9]/(2*delta_muev)*delta_muev/(110);
%w_vec = w_vec([1:2,4:end]);
w_vec = sort([w_vec, w_special]);


%% Call Constants

k_f_inverse = 40e-9; %40e-9; % m
a = (pi/2)*k_f_inverse; % m 

%% Diagonalise the Hamiltonians
mzm_cell = {};
for ii = 1:size(w_vec,2)
    fprintf('%i ', ii);
    if ~mod(ii, 10)
        fprintf('\n');
    end
    % [e_vecs_1, e_vals_1, e_vecs_2, e_vals_2, e_vecs_2_init_b, e_vals_2_init_b, ...
    % majorana_zero_modes_ref] = ...
    %     diagonalise_TLS_Hamiltonians(mu_low, mu_high, mu_offset, w_vec(ii), delta, N, BC);
    % mzm_cell{end+1} = majorana_zero_modes_ref;
    H_tetron = get_tetron_BdG_Hamiltonian(mu, w_vec(ii), delta, N, BC, 'chain_1_chain_2');
    [e_vecs_1, e_vals_1, majorana_zero_modes_ref] = ... 
        diagonalise_uncoupled_tetron_via_Kitaev_chains(H_tetron(1:2*N,1:2*N), ...
        H_tetron((2*N+1):end, (2*N+1):end), 'dirac_zero_modes');  
    mzm_cell{end+1} = majorana_zero_modes_ref;
end


%% Extra Plots
mzm_length_vec = zeros(length(w_vec), 1);

for ii = 1:length(w_vec)

    lat_sites_odd = (1:2:floor(N/2)).';
    x_odd = (lat_sites_odd-1)*a*1e6 ; % \mu m 
    lat_sites_all = (1:1:(N)).';
    x = (lat_sites_all-1)*a*1e6; % \mu m  
    % 
    mzm_temp = mzm_cell{ii};
    mzm_temp_up = abs(mzm_temp(lat_sites_odd,1));
    mzm_temp_up_l = log(mzm_temp_up);
    ind_fit = mzm_temp_up_l > -25;
    lat_sites_fit = lat_sites_odd(ind_fit);
    x_fit = (lat_sites_fit-1)*a*1e6; %\mu m
    
    figure();
    plot(x, abs(mzm_temp(1:N,1)), 'o'); hold on;
    plot(x_odd, abs(mzm_temp_up), 'x');
    title(sprintf('w = %.3g', w_vec(ii)));
   % Do the fitting here

    figure();
    plot(x_odd, log(mzm_temp_up), 'o')
    P = polyfit(x_fit, log(mzm_temp_up(ind_fit)), 1);
    hold on;
    plot(x_odd, P(1)*x_odd + P(2), 'r');
    mzm_length_vec(ii) = -1./P(1);
    title(sprintf('w = %.3g = %.3g \\mu eV', w_vec(ii), w_vec(ii)*2*delta_muev));

end

mzm_length_vec_nm = mzm_length_vec*1e3;

%%


w_vec_mueV = w_vec*2*delta_muev;
%%
figure();
plot(w_vec_mueV, mzm_length_vec_nm, '-.o');
xlabel('w_{Alicea} [\mu eV]'); ylabel('MZM Length [nm]');
hold on;
%plot(w_vec_mueV([3,10]), mzm_length_vec_nm([3,10]), '*');


fit_loc= polyfit(w_vec_mueV, mzm_length_vec_nm, 1);
%hold on;
%plot(w_vec_mueV, w_vec_mueV*fit_loc(1) + fit_loc(2), 'r-');

% w for 250nm coherence length
w_max_mueV = (250 - fit_loc(2))./fit_loc(1)
w_min_mueV = (100 - fit_loc(2))./fit_loc(1)

w_max = w_max_mueV/(2*delta_muev)
w_min = w_min_mueV/(2*delta_muev)


%%
% - - - - - - 
% PLOT SCRIPT - MZM_FINDER_PLOTS_FOR_FIGURE
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


w_vec_mueV = 2*delta_muev*w_vec;

%% Extra Point - what I actually solved For

% w_mueV_solve = 350.79;
% mzm_length_nm_solve = 121.8017; %

%%

figure();
plot(w_vec_mueV, mzm_length_vec_nm, '-.o');
xlabel('w_{Alicea} [\mu eV]'); ylabel('MZM Length [nm]');
hold on;

%%
% === Styling parameters ===
marker_size = 17;
line_width = 4;
fig_ax_font_size = 30;
fig_border_LW = 2.3;

% === Colors ===
RGB = orderedcolors("gem");
col_1 = RGB(1,:);  % or choose your preferred color
col_1_light = col_1 + (1 - col_1) * 0.3;

col_1 = RGB(1,:);  % 3 μm
col_2 = RGB(2,:);  % 5 μm
col_3 = RGB(3,:)*0.8;  % 10 μm
 

% === Create figure ===
currFig = figure();
currFig.Position = [205 258 786 503];

ax = gca;
ax.FontSize = fig_ax_font_size;
ax.LineWidth = fig_border_LW;
hold on;

xlim([0, 150]);
ylim([0, 300]);

% === Fit === %
loc_fit = polyfit(w_vec_mueV, mzm_length_vec_nm,1);
w_fit = linspace(0, 800, 10); 
y_fit = loc_fit(1)*w_fit + loc_fit(2);

% === Plot data ===
h_fit_mzm = plot(w_fit, y_fit, '-', ...
    'Color', col_1_light, ...
    'MarkerSize', marker_size, ...
    'LineWidth', 0.8*line_width, ...
    'MarkerFaceColor', col_1, ...
    'MarkerEdgeColor', 'w');
hold on;

chosen_ind =  [4, 18];
solve_ind = 12;
%hide_ind = chosen_ind + 1;
hide_ind = []; %[3, 13, 17, 19];

w_mueV_solve = w_vec_mueV(solve_ind);
mzm_length_nm_solve = mzm_length_vec_nm(solve_ind); 

w_chosen = w_vec_mueV(chosen_ind);
mzm_length_chosen = mzm_length_vec_nm(chosen_ind);

w_vec_mueV_reduced = w_vec_mueV;
mzm_length_vec_nm_reduced = mzm_length_vec_nm;
w_vec_mueV_reduced([chosen_ind, hide_ind, solve_ind]) = [];
mzm_length_vec_nm_reduced([chosen_ind, hide_ind, solve_ind]) = [];

h_mzm = plot(w_vec_mueV_reduced(:), mzm_length_vec_nm_reduced(:), 'o', ...
    'Color', col_1, ...
    'MarkerSize', marker_size*0.8, ...
    'LineWidth', line_width, ...
    'MarkerFaceColor', col_1, ...
    'MarkerEdgeColor', 'w');
hold on; 
h_mzm = plot(w_chosen(:), mzm_length_chosen(:), '*', ...
    'Color', col_3, ...
    'MarkerSize', marker_size, ...
    'LineWidth', line_width, ...
    'MarkerFaceColor', 'w', ...
    'MarkerEdgeColor', col_3);

h_chosen = plot(w_mueV_solve, mzm_length_nm_solve, '*', ...
    'Color', col_3, ...
    'MarkerSize', marker_size, ...
    'LineWidth', line_width, ...
    'MarkerFaceColor', 'w', ...
    'MarkerEdgeColor', col_2);


% Optional: highlight specific points
% plot(w_vec_mueV([3,10]), mzm_length_vec_nm([3,10]), '*', ...
%     'Color', col_1, 'MarkerSize', marker_size * 1.2, 'LineWidth', line_width);

% === Axis labels ===
xlabel('Kitaev chain hopping strength, $w$ ($\mu$eV)', 'Interpreter', 'latex');
ylabel({'MZM Localization Length,', '$\zeta$ (nm)'}, 'Interpreter', 'latex');

% === Ticks and limits ===
ax.TickDir = 'in';
ax.Box = 'off';  % we'll draw box manually
ax.XMinorTick = 'off';
ax.YMinorTick = 'off';
% ax.YAxis.MinorTickValues = 1:10; % optional if needed

% === Draw custom box ===
xl = xlim(ax);
yl = ylim(ax);
plot([xl(1) xl(2)], [yl(2) yl(2)], 'k', 'LineWidth', ax.LineWidth) % top
plot([xl(2) xl(2)], [yl(1) yl(2)], 'k', 'LineWidth', ax.LineWidth) % right


% First special point
text(3, 110, ...
    ['$\begin{array}{l}' ...
     'w=' num2str(w_chosen(1), '%.0f') '\,\mu \rm{eV} \\ ' ...
     '\zeta=' num2str(mzm_length_chosen(1), '%.0f') '\,\rm{nm}' ...
     '\end{array}$'], ...
    'FontSize', fig_ax_font_size*0.8 , ...
    'Interpreter', 'latex', ...
    'Color', col_3);

% Second special point
text(55, 260, ...
    ['$\begin{array}{l}' ...
     'w=' num2str(w_chosen(2), '%.0f') '\,\mu \rm{eV} \\ ' ...
     '\zeta=' num2str(mzm_length_chosen(2), '%.0f') '\,\rm{nm}' ...
     '\end{array}$'], ...
    'FontSize', fig_ax_font_size*0.8 , ...
    'Interpreter', 'latex', ...
    'Color', col_3);


% Third special point
text(340/4.685, 160, ...
    ['$\begin{array}{l}' ...
     'w=' num2str(w_mueV_solve, '%.1f') '\,\mu \rm{eV} \\ ' ...
     '\zeta=' num2str(mzm_length_nm_solve, '%.0f') '\,\rm{nm}' ...
     '\end{array}$'], ...
    'FontSize', fig_ax_font_size*0.8 , ...
    'Interpreter', 'latex', ...
    'Color', col_2);


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% Save Figure S3

% Save as MATLAB .fig image
saveas(gcf, 'figure_S3.fig');

% Save as PNG and SVG image
saveas(gcf, 'figure_S3.png');
saveas(gcf, 'figure_S3.svg');


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%% Save Data
%close all;
save('figure_S3_data');


%% Export data for Figure S3 to CSV

output_dir = ".";
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

filename = fullfile(output_dir, "Figure_S3.csv");

% --- Fit line data (exact plotted line) ---
loc_fit = polyfit(w_vec_mueV, mzm_length_vec_nm, 1);
w_fit = linspace(0, 800, 10);
y_fit = loc_fit(1)*w_fit + loc_fit(2);

% --- Special indices from plotting script ---
chosen_ind = [4, 18];
solve_ind = 12;
hide_ind = [];

% --- Reproduce reduced vectors exactly as plotted ---
w_vec_mueV_reduced = w_vec_mueV;
mzm_length_vec_nm_reduced = mzm_length_vec_nm;

w_vec_mueV_reduced([chosen_ind, hide_ind, solve_ind]) = [];
mzm_length_vec_nm_reduced([chosen_ind, hide_ind, solve_ind]) = [];

% --- Chosen and solve points ---
w_chosen = w_vec_mueV(chosen_ind);
mzm_length_chosen = mzm_length_vec_nm(chosen_ind);

w_mueV_solve = w_vec_mueV(solve_ind);
mzm_length_nm_solve = mzm_length_vec_nm(solve_ind);

% --- Build table ---
T = table( ...
    w_vec_mueV(:), ...
    mzm_length_vec_nm(:), ...
    'VariableNames', { ...
        'w_u_eV', ...
        'MZM_length_nm' ...
    });

% --- Add fit line as separate columns (padded with NaN to match length) ---
N = height(T);

fit_w_col = nan(N,1);
fit_y_col = nan(N,1);

fit_w_col(1:length(w_fit)) = w_fit(:);
fit_y_col(1:length(y_fit)) = y_fit(:);

T.w_fit_u_eV = fit_w_col;
T.fit_length_nm = fit_y_col;

% --- Add chosen points ---
chosen_w_col = nan(N,1);
chosen_len_col = nan(N,1);

chosen_w_col(1:length(w_chosen)) = w_chosen(:);
chosen_len_col(1:length(mzm_length_chosen)) = mzm_length_chosen(:);

T.w_chosen_u_eV = chosen_w_col;
T.length_chosen_nm = chosen_len_col;

% --- Add solve point ---
solve_w_col = nan(N,1);
solve_len_col = nan(N,1);

solve_w_col(1) = w_mueV_solve;
solve_len_col(1) = mzm_length_nm_solve;

T.w_solve_u_eV = solve_w_col;
T.length_solve_nm = solve_len_col;

% --- Write CSV ---
writetable(T, filename);

disp("CSV file written:");
disp(filename);