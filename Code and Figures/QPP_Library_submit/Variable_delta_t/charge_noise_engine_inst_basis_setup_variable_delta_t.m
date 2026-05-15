%--------------------------------------%
% Charge Noise Engine - General - INIT CONDITION
%--------------------------------------%
% DO NOT USE THIS FOR ACTUAL SIMULATIONS. ONLY USED FOR CHECKS IN SETUP.
% ^ I think this statement got overwritten.
%
% General version of the charge noise engine for different inbuilt noise
% types OR with you own imported noise trajectory. 
%
%
% INPUTS:
%   noise_type: A string denoting the noise type. Allowed inputs:
%               'gaussian', 'two_state_exp', or 'custom'
%   noise_parameters: An array or matrix, depending on the the entry on the
%               noise_type. For: 
%               Gaussian: [t_noise_factor, fluctuation_amplitude, characteristic_noise_corr_freq]
%                       two_state_exp: [sigma, tau, fluctuation_amplitude]
%               custom: matrix of size "num_time_steps x 2*num_trials", whose
%                       columns contain independent noise trajectories. The
%                       first num_trial columns are for the top kitaev
%                       chain and the following num_trial columns are for
%                       the bottom Kitaev chain. 
%
%   time_evolution_method: Choose method for obtaining the unitary U_t generating
%                          time evolution at each time step, from the noisy
%                          Hamiltonian, H_t, at that time. 
%               'time_ev_diagonalise' - Diagonalise H_t to construct U_t.
%               'time_ev_exp' - Use matlab function "expm" to obtain U_t
%               'time_ev_1st_order' - Use a 1st order Taylor series
%                              expansion to obtain U_t.
%               'time_ev_Magnus' - use the Magnus expansion - TO BE IMPLEMENTED.  
%
%   save_output = 'save_output' or 'do_not_save_output'
%
%   plot_results = 'plot_all_results', 'plot_X_exp_only' or 'do_not_plot_results'
%
%   noise_trajectories = 1st element of varargin. 2*num_trial x
%       num_time_steps_noise matrix, where each row is a noise realisation. 
%       num_time_steps_noise = (t_final-t_init)/delta_t. First num_trial rows
%       are for the top Kitaev chain, the next num_trial rows are for the 
%       bottom Kitaev chain. If this input is not provided, then a new noise 
%       trajectory is generated. 
%
%   init_condition - Select whether at t=0 the system is initialised 
%                       according to the clean Hamiltonian or the first time 
%                       step of the noise-trajectory. 
%                       Allowed values are:
%                       'initial_state_no_noise', 'initial_state_noise'.
%
%   basis_choice - Select whether to define the Pauli operators based of the 
%                  Hamiltonian at t = 0 or the instanteneous Hamiltonian
%                  H(t). Allowed values are:
%                  'init_hamil_basis', 'instant_hamil_basis'.
%
% OUTPUTS
%   e_vals_mat - this contains the eigenvalues only for the first trajectory. 
%                I don't want to pass all trajectory's eigenvalues out. 
%                The (i, j)'th element corresponds to the j'th eigenvalue
%                at the t_i time step. 
%
% To do - change above plot parameter to distinguish, the full plots of
% individual trajectories and everythign, or just the X_exp mean results
% and fit.
%
% INIT Eigenvectors are both the initial eigenvectors and the "averaged" Hamiltonian
% % eigenvectors - ie. with no noise.
% 
    
function [omega_0_from_fit, T2_from_fit, omega_0_guess, X_exp_mean,...
         P_y1y2y3y4_v1_mean, P_y1y2y3y4_v2_mean, overlap_with_init_state_mean, overlap_fixed_basis_states_mean, ...
         t_vec, top_power_mat, bot_power_mat, mzm_parity_decay, ...
         comparison_determinants_cell, gap_mean, e_vals_mat, quasi_occ_t, ...
         overlap_inst_basis_states_mean, Z_exp_mean] = ...
    charge_noise_engine_inst_basis_setup_variable_delta_t(mu_mean, mu_offset, w, delta, N, BC, t_init, t_final,...
     delta_t, noise_type, noise_parameters, num_trials, time_evolution_method,...
     save_output, plot_results, init_condition, basis_choice, varargin)  

    check_valid_inputs(save_output, plot_results);
  

     %% Varargin Processing (added April 10th - Late Addition)
     %Current optional variables (in order): 
     % [init_state] - allowed values: "0", "1", "+", "-"
     % Feel free to extend varargin to more inputs
     if nargin == 17 % Default parameters
         init_state = "+"; 
     elseif nargin == 18
        init_state = varargin{1}; % No need to check if valid here, is checked in get_tetron_init_cov_mat_general. 
     elseif nargin > 19
             error("Too many input parameters");
     end
     %% Process Input Parameters
     mu = mu_mean + mu_offset*[+1, -1];
    
     %% Get Initial Hamiltonian - if defined by clean Hamiltonian
    if strcmp(init_condition, 'initial_state_no_noise')
        H_tetron_init = get_tetron_BdG_Hamiltonian(mu, w, delta, N, BC, 'chain_1_chain_2');
        [e_vecs_init, e_vals_init, majorana_zero_modes_ref] = ... 
            diagonalise_uncoupled_tetron_via_Kitaev_chains(H_tetron_init(1:2*N,1:2*N), ...
            H_tetron_init((2*N+1):end, (2*N+1):end), 'dirac_zero_modes');
    
        %[~, ~, corr_qp_init] = ...
        %     get_tetron_init_cov_mat_general(e_vecs_init, "+");  
        [~, ~, corr_qp_init] = ...
            get_tetron_init_cov_mat_general(e_vecs_init, init_state); 

    end

    %% Run Simulations
    t_vec = t_init:delta_t:t_final;
    Z_exp_mat = zeros(num_trials, length(t_vec));
    X_exp_mat = zeros(num_trials, length(t_vec));
    Y_exp_mat = zeros(num_trials, length(t_vec));
    P_y1y2y3y4_v1_mat = zeros(num_trials, length(t_vec));
    P_y1y2y3y4_v2_mat = zeros(num_trials, length(t_vec));    
    mu_mat_top = zeros(num_trials, length(t_vec));
    mu_mat_bot = zeros(num_trials, length(t_vec));
    top_power_mat = zeros(num_trials,1);
    bot_power_mat = zeros(num_trials,1);
    comparison_determinants_mat_top = zeros(num_trials, length(t_vec));
    comparison_determinants_mat_bot = zeros(num_trials, length(t_vec));

    %warning('comparison_determinant_mat only calcd for one traj');

    overlap_with_init_state_mat = zeros(num_trials, length(t_vec));
    overlap_fixed_basis_states_mat = zeros(num_trials, length(t_vec));
    overlap_inst_basis_states_mat = zeros(num_trials, length(t_vec));

    tStart = tic;

    if num_trials > 1
        fprintf('Trial Number: \n');
    end
    for trial_no = 1:num_trials
        if trial_no > 1
            fprintf('%i ', trial_no);
            if mod(trial_no, 10) == 0
                fprintf('\n');
            end
        end
        % % Get Noise Trajectories for top and bottom chains
        % [gaussian_noise_top, t_noise_top, white_noise_top] =...
        %     generate_gaussian_noise(t_init, t_final, delta_t, ...,
        %     t_noise_factor, fluctuation_amplitude, characteristic_noise_corr_freq);
        % [gaussian_noise_bot, t_noise_bot, white_noise_bot] =...
        %     generate_gaussian_noise(t_init, t_final, delta_t, ...,
        %     t_noise_factor, fluctuation_amplitude, characteristic_noise_corr_freq);
        % 
        % % Grab a centre segment of the noise trajectories.
        % t_ind_start = floor(length(t_noise_top)/4);
        % t_ind_final = t_ind_start + length(t_vec) - 1;
        %gaussian_noise_top(t_ind_start:t_ind_final);
       
        if ~strcmp(noise_type, 'custom')
            [noise_top, noise_bot, fluctuation_amplitude, top_power, bot_power] = ...
                process_noise_generation(noise_type, noise_parameters, ...
            t_init, t_final, delta_t);
        else %Custom noise
            % num_time_steps x 2*num_trials
            noise_top = noise_parameters(:, trial_no).';
            noise_bot = noise_parameters(:, trial_no + num_trials).';
            top_power = sqrt(sum(abs(noise_top).^2).*(t_vec(2) - t_vec(1)));
            bot_power = sqrt(sum(abs(noise_bot).^2).*(t_vec(2) - t_vec(1)));
            fluctuation_amplitude = sqrt(top_power);             
        end

        top_power_mat(trial_no) = top_power;
        bot_power_mat(trial_no) = bot_power;

        % Store the top/bot noise trajectories in a matrix
        mu_vec = zeros(2, length(t_vec));
        mu_vec(1,:) = mu_mean + mu_offset +  ...
            noise_top;
        
        mu_vec(2,:) = mu_mean - mu_offset +  ...
            noise_bot;
       
        mu_mat_top(trial_no, :) = mu_vec(1,:);
        mu_mat_bot(trial_no, :) = mu_vec(2,:);
  
        % Get initial Hamiltonian - if defining by first step of noise
        % trajectory
        if strcmp(init_condition, 'initial_state_noise')
            H_tetron_init = get_tetron_BdG_Hamiltonian(mu_vec(:,1), w, delta, N, BC, 'chain_1_chain_2');
            [e_vecs_init, e_vals_init, majorana_zero_modes_ref] = ... 
                diagonalise_uncoupled_tetron_via_Kitaev_chains(H_tetron_init(1:2*N,1:2*N), ...
                H_tetron_init((2*N+1):end, (2*N+1):end), 'dirac_zero_modes');
        
            % [~, ~, corr_qp_init] = ...
            %      get_tetron_init_cov_mat_general(e_vecs_init, "+");   
            [~, ~, corr_qp_init] = ...
                 get_tetron_init_cov_mat_general(e_vecs_init, init_state);     
        end
        gap_mean = sort(abs(e_vals_init)); 
        gap_mean = gap_mean(5); %omit zero modes

        % Do time evolution
        tetron_params = {w, delta, N, BC};
    
        % Dodgey handling of the method I used to calculate
        % overlap_inst_basis_states
        if nargout >= 17
            [Z_exp, X_exp, Y_exp, P_y1y2y3y4_v1, P_y1y2y3y4_v2, overlap_with_init_state, ...
                overlap_fixed_basis_states, ~,~,~, comparison_determinants_mat, ...
                e_vals_mat_curr, ~, quasi_occ_t, ...
                overlap_inst_basis_states] =...
                calc_noisy_tetron_ev_inst_basis_setup(...
                corr_qp_init, e_vecs_init, e_vals_init, tetron_params,...
                mu_vec, delta_t, t_init, t_final, time_evolution_method, ...
                basis_choice, majorana_zero_modes_ref, mu);
        else
            [Z_exp, X_exp, Y_exp, P_y1y2y3y4_v1, P_y1y2y3y4_v2, overlap_with_init_state, ...
                overlap_fixed_basis_states, ~,~,~, comparison_determinants_mat, ...
                e_vals_mat_curr, ~, quasi_occ_t] =...
                calc_noisy_tetron_ev_inst_basis_setup(...
                corr_qp_init, e_vecs_init, e_vals_init, tetron_params,...
                mu_vec, delta_t, t_init, t_final, time_evolution_method, ...
                basis_choice, majorana_zero_modes_ref, mu);

                overlap_inst_basis_states = NaN*ones(size(overlap_with_init_state));
        end

        if trial_no ==1
            e_vals_mat = e_vals_mat_curr; % Only save e_vals for first trial
        end

        comparison_determinants_mat_top(trial_no, :) = comparison_determinants_mat(1, :);
        comparison_determinants_mat_bot(trial_no, :) = comparison_determinants_mat(2, :);

        Z_exp_mat(trial_no, :) = Z_exp;
        X_exp_mat(trial_no, :) = X_exp;
        Y_exp_mat(trial_no, :) = Y_exp;
        P_y1y2y3y4_v1_mat(trial_no, :) = P_y1y2y3y4_v1;
        P_y1y2y3y4_v2_mat(trial_no, :) = P_y1y2y3y4_v2;        
        overlap_with_init_state_mat(trial_no, :) = overlap_with_init_state;
        overlap_fixed_basis_states_mat(trial_no, :) = overlap_fixed_basis_states;
        overlap_inst_basis_states_mat(trial_no, :) = overlap_inst_basis_states;        
       
    end

    comparison_determinants_cell = {comparison_determinants_mat_top, ...
        comparison_determinants_mat_bot};

    fprintf('\n');
    run_time = toc(tStart);

    % Clean sim:
    H_evolve = H_tetron_init; H_evolve_qp = e_vecs_init'*H_tetron_init*e_vecs_init; % correct ordering

    [Z_exp_clean, X_exp_clean, Y_exp_clean] = calculate_tetron_evolution_dirac_qp_corr_mat(...
         corr_qp_init, e_vecs_init, e_vals_init, H_evolve_qp, delta_t, t_init, t_final);
    
    %% Do Fit on X_exp Data
    % Fit a damped cosine to the 

    X_exp_mean = mean(X_exp_mat,1); 
    Z_exp_mean = mean(Z_exp_mat,1); 
    P_y1y2y3y4_v1_mean = mean(P_y1y2y3y4_v1_mat,1);
    P_y1y2y3y4_v2_mean = mean(P_y1y2y3y4_v2_mat,1);  
    overlap_with_init_state_mean = mean(overlap_with_init_state_mat,1);
    overlap_fixed_basis_states_mean = mean(overlap_fixed_basis_states_mat,1);
    overlap_inst_basis_states_mean = mean(overlap_inst_basis_states_mat,1);

    if norm(imag(X_exp_mean)) > 1e-10
        warning('non-negligble imaginary componenent in X_exp_mean');
    end
    X_exp_mean = real(X_exp_mean); 

    cosine_damped = 'exp(-a*x)*cos(b*x)';

    %Initial guess for a and b in equation above
    start_point = [fluctuation_amplitude*1e-2, 2*min(abs(e_vals_init))];
    %weights = linspace(1,0,length(t_vec)).^2;
    weights = ones(size(t_vec));
    try
        X_exp_fit = fit(t_vec(:), X_exp_mean(:), cosine_damped, 'Start', start_point,...
                        'Weights', weights);
        best_fit = exp(-X_exp_fit.a.*t_vec).*cos(X_exp_fit.b.*t_vec);
        omega_0_from_fit = X_exp_fit.b; 
        T2_from_fit = 1./X_exp_fit.a;
    catch
        warning('Fit failed');
        T2_from_fit = NaN;
        omega_0_from_fit = NaN;
        best_fit = zeros(size(t_vec));
    end

    if strcmp(plot_results, 'plot_X_exp_only') ||...
            strcmp(plot_results, 'plot_all_results') 
        figure();
        plot(t_vec, X_exp_mean,'LineWidth', 1);
        xlabel('t');ylabel('<X> = <-2i \gamma_1 \gamma_3>');
        grid on;
        hold on; plot(t_vec, best_fit, '-.', 'LineWidth', 1);
        if T2_from_fit >=0 
            plot([T2_from_fit,T2_from_fit], [-1,1], '-.', 'LineWidth', 1);
        else 
            warning('Negative T2 time calculated');
        end
        hold off;
        legend('Numerics', 'fit', 'T2 time');
    end
    
    %% Calculate Decay time of Total MZM Parity

    % Constrained polynomial of log(parity);
    linear_symb = '-a*x';
    %Initial guess for a and b in equation above
    start_point = [-1/1000];
    %weights = linspace(1,0,length(t_vec)).^2;
    weights = ones(size(t_vec));
    log_mzm_parity = log(P_y1y2y3y4_v1_mean);
    parity_leakage_fit = fit(t_vec(:), log_mzm_parity(:), linear_symb, 'Start', start_point,...
                    'Weights', weights);

    mzm_parity_decay = 1./parity_leakage_fit.a; 
    

    %% Plot results
% 
    if strcmp(plot_results, 'plot_all_results')
        figure(); subplot(3,2,1);
        plot(t_vec, real(Z_exp_mat), '-.'); hold on;
        plot(t_vec, Z_exp_clean, 'k', 'LineWidth', 1.5);
        xlabel('t');ylabel('<Z> = <-2i \gamma_1 \gamma_2>');
        grid on;
        legend('clean');

        subplot(3,2,3);
        plot(t_vec, real(X_exp_mat), '-.'); hold on;
        plot(t_vec, X_exp_clean, 'k', 'LineWidth', 1.5); 
        xlabel('t');ylabel('<X> = <-2i \gamma_1 \gamma_3>');
        grid on;
        legend('clean');

        subplot(3,2,5); 
        plot(t_vec, real(Y_exp_mat), '-.'); hold on;
        plot(t_vec, Y_exp_clean, 'k', 'LineWidth', 1.5);
        xlabel('t'); ylabel('<Y> = <2i \gamma_1 \gamma_4>');
        grid on;
        shg;
        legend('clean');

        subplot(3,2,2); 
        plot(t_vec, mean(real(Z_exp_mat),1))
        xlabel('t');ylabel('<Z> = <-2i \gamma_1 \gamma_2>');
        grid on;
        legend('Mean');

        subplot(3,2,4);
        plot(t_vec, mean(real(X_exp_mat),1));
        xlabel('t');ylabel('<X> = <-2i \gamma_1 \gamma_3>');
        grid on;
        legend('Mean');

        subplot(3,2,6);
        plot(t_vec, mean(real(Y_exp_mat),1));
        xlabel('t'); ylabel('<Y> = <2i \gamma_1 \gamma_4>');
        grid on;
        shg;
        legend('Mean');

        figure(); subplot(1,2,1);
        plot(t_vec, mu_mat_top); xlabel('t'); ylabel('\mu_{top}');
        grid on;
        subplot(1,2,2);
        plot(t_vec, mu_mat_bot); xlabel('t'); ylabel('\mu_{bot}');
        grid on;

        figure(); 
        plot(t_vec, P_y1y2y3y4_v1_mean);
        xlabel('t'); ylabel('<P_{total}> = <-4 \gamma_1 \gamma_2 \gamma_3 \gamma_4');
        hold on;
        plot(t_vec, exp(-t_vec./mzm_parity_decay));
    end

     %Debug - but is outputted
    a_temp = sort(e_vals_init, 'ascend', 'ComparisonMethod','abs');
    omega_0_guess = a_temp(1).*2;

    %% Save Output

    % Format a nice string which gives all the parameters and saves them in
    % this cell. 
    if strcmp(save_output, 'save_output')
        noise_type = 'Gaussian';
        curr_time = datetime(now,'ConvertFrom','datenum');

        % Add Run-time to info_str
        info_str = sprintf(['Tetron Charge Noise Results\nTime: %s \n', ...
            'Run time: %g seconds\n',...
            'Noise Type: %s \n\n',...
            'N = %i \nw = %g \ndelta = %g \nBC = %s',...
            '\nmu_mean = %g \nmu_offset = %g', ...
            '\n\nfluctuation_amplitude = %g, \ncharacteristic_noise_corr_freq = %g', ...
            '\nt_noise_factor = %g', ...
            '\nnum_trials = %g', ...
            '\n\nt_init = %g, \nt_final = %g \ndelta_t = %g', ...
            '\n\nCell array output elements, in order are:\n',...
            't_vec, Z_exp_mat, X_exp_mat, Y_exp_mat\n',...
            'Z_exp_clean, X_exp_clean, Y_exp_clean, mu_mat_top, mu_mat_bot.'],...
            curr_time, run_time, noise_type, N, w, delta, BC, mu_mean, mu_offset,...
            fluctuation_amplitude, characteristic_noise_corr_freq, t_noise_factor, ...
            num_trials, t_init, t_final, delta_t );

        output = {info_str, t_vec, Z_exp_mat, X_exp_mat, Y_exp_mat, ...
            Z_exp_clean, X_exp_clean, Y_exp_clean, mu_mat_top, mu_mat_bot};

        save('./Results/trajectories_run_generic', "output");

        %% Save with next available file name count
        % Should put a catch if this file doesn't exist:
        load('traj_run_count');
        traj_run_count = traj_run_count +1;
        % traj_run_count = 2;
        save('traj_run_count', 'traj_run_count');
        save_file_name = sprintf('./Results/trajectories_run_%i', traj_run_count);

        save(save_file_name, "output");

        traj_run_count = 1;
        save('traj_run_count', 'traj_run_count');
    end
end

%% Function - Check Valid Inputs
function check_valid_inputs(save_output, plot_results)

    if ~strcmp(save_output, 'save_output') & ...
            ~strcmp(save_output, 'do_not_save_output')
        error('Invalid Input for save_output');
    end
    
    if ~strcmp(plot_results, 'plot_all_results') & ...
           ~strcmp(plot_results, 'plot_X_exp_only') &...
            ~strcmp(plot_results, 'do_not_plot_results') 
        error('Invalid Input for plot_results');
    end

end

%% Function - Process Noise Generation

function [noise_top, noise_bot, fluctuation_amplitude, top_power, bot_power] = ...
    process_noise_generation(noise_type, noise_parameters, ...
    t_init, t_final, delta_t)

    t_vec = t_init:delta_t:t_final;
    dt = t_vec(2) - t_vec(1);         
    total_time = max(t_vec) - min(t_vec);

    switch noise_type
        case 'gaussian'
            %[t_noise_factor, fluctuation_amplitude, characteristic_noise_corr_freq]
            t_noise_factor = noise_parameters(1);
            fluctuation_amplitude = noise_parameters(2);
            characteristic_noise_corr_freq = noise_parameters(3); 
    
            [gaussian_noise_top, t_noise_top] =...
                generate_gaussian_noise(t_init, t_final, delta_t, ...,
                t_noise_factor, fluctuation_amplitude, characteristic_noise_corr_freq);
            [gaussian_noise_bot, t_noise_bot] =...
                generate_gaussian_noise(t_init, t_final, delta_t, ...,
                t_noise_factor, fluctuation_amplitude, characteristic_noise_corr_freq);
    
            % Grab a centre segment of the noise trajectories.
            t_ind_start = floor(length(t_noise_top)/4);
            t_ind_final = t_ind_start + length(t_vec) - 1;
            noise_top = gaussian_noise_top(t_ind_start:t_ind_final);  
            noise_bot = gaussian_noise_bot(t_ind_start:t_ind_final);  

        case 'two_state_exp'
            %[sigma, tau, signal_power]
            sigma = noise_parameters(1);
            tau = noise_parameters(2);
            signal_power = noise_parameters(3);

            [noise_top, t_noise_top] = ...
                generate_two_parameter_noise(sigma, tau, t_init, t_final, delta_t, signal_power, 'do_not_remove_dc');
            [noise_bot, t_noise_bot] = ...
                 generate_two_parameter_noise(sigma, tau, t_init, t_final, delta_t, signal_power, 'do_not_remove_dc');
            warning('DC not reomved');
            %noise_top = fluctuation_amplitude.*gaussian_noise_top;
            %noise_bot = fluctuation_amplitude.*gaussian_noise_bot;
            fluctuation_amplitude = sqrt(signal_power);

        case 'two_state_exp_no_dc'
            %[sigma, tau, signal_power]
            sigma = noise_parameters(1);
            tau = noise_parameters(2);
            signal_power = noise_parameters(3);

            [noise_top, t_noise_top] = ...
                generate_two_parameter_noise(sigma, tau, t_init, t_final, delta_t, signal_power, 'remove_dc');
            [noise_bot, t_noise_bot] = ...
                 generate_two_parameter_noise(sigma, tau, t_init, t_final, delta_t, signal_power, 'remove_dc');
            fluctuation_amplitude = sqrt(signal_power);
            %gaussian_noise_top = gaussian_noise_top - 0.5; % remove dc component
            %gaussian_noise_bot = gaussian_noise_bot - 0.5; 
            %noise_top = fluctuation_amplitude.*gaussian_noise_top;
            %noise_bot = fluctuation_amplitude.*gaussian_noise_bot;  

            % %DEBUG
            % t_vec = t_init:delta_t:t_final;
            % dt = t_vec(2) - t_vec(1);            
            % power_2state = sum(abs(noise_top).^2).*dt./max(t_vec)

        case 'lorentzian_white_filtered'      
            % [sigma, tau, signal_power, t_noise_factor]
            sigma = noise_parameters(1);
            tau = noise_parameters(2);
            signal_power = noise_parameters(3);
            t_noise_factor = noise_parameters(4);

            [noise_top, t_noise_top] = ...
                 generate_filtered_lorentzian_noise(t_init, t_final, delta_t, t_noise_factor,...
                sigma, tau, signal_power);
            [noise_bot, t_noise_top] = ...
                 generate_filtered_lorentzian_noise(t_init, t_final, delta_t, t_noise_factor,...
                sigma, tau, signal_power);
            fluctuation_amplitude = sqrt(signal_power);

            % Grab a centre segment of the noise trajectories.
            %t_vec = t_init:delta_t:t_final;
            t_ind_start = floor(length(t_noise_top)/4);
            t_ind_final = t_ind_start + length(t_vec) - 1;
            noise_top = noise_top(t_ind_start:t_ind_final);  
            noise_bot = noise_bot(t_ind_start:t_ind_final);  

        case 'lorentzian_white_filtered_cut_off'      
            % [sigma, tau, signal_power, t_noise_factor, lor_cut_off_freq]
            sigma = noise_parameters(1);
            tau = noise_parameters(2);
            signal_power = noise_parameters(3);
            t_noise_factor = noise_parameters(4);
            lor_cut_off_freq = noise_parameters(5);

            [noise_top, t_noise_top] = ...
                 generate_filtered_lorentzian_noise_cut_off(t_init, t_final, delta_t, t_noise_factor,...
                sigma, tau, signal_power, lor_cut_off_freq);
            [noise_bot, t_noise_top] = ...
                 generate_filtered_lorentzian_noise_cut_off(t_init, t_final, delta_t, t_noise_factor,...
                sigma, tau, signal_power, lor_cut_off_freq);
            fluctuation_amplitude = sqrt(signal_power);

            % Grab a centre segment of the noise trajectories.
            %t_vec = t_init:delta_t:t_final;
            t_ind_start = floor(length(t_noise_top)/4);
            t_ind_final = t_ind_start + length(t_vec) - 1;
            noise_top = noise_top(t_ind_start:t_ind_final);  
            noise_bot = noise_bot(t_ind_start:t_ind_final);  

        case 'combined_gaussian_lorentzian' 
        % [sigma, tau, signal_power_lor, signal_power_gauss, ...
        % noise_corr_freq_gauss, signal_power_total, t_noise_factor, ...
        % cut_off_freq]
        sigma = noise_parameters(1);
        tau = noise_parameters(2);
        signal_power_lor = noise_parameters(3);
        signal_power_gauss = noise_parameters(4);
        noise_corr_freq_gauss = noise_parameters(5);
        signal_power_total = noise_parameters(6);
        t_noise_factor = noise_parameters(7);
        cut_off_freq = noise_parameters(8);

            [noise_top, t_noise_top] = ...
                generate_gauss_lor_tail_noise(t_init, t_final, delta_t, t_noise_factor,...
                sigma, tau, signal_power_lor, signal_power_gauss, noise_corr_freq_gauss, ...
                signal_power_total, cut_off_freq);

            [noise_bot, t_noise_bot] = ...
                generate_gauss_lor_tail_noise(t_init, t_final, delta_t, t_noise_factor,...
                sigma, tau, signal_power_lor, signal_power_gauss, noise_corr_freq_gauss, ...
                signal_power_total, cut_off_freq);

        fluctuation_amplitude = sqrt(signal_power_total);

        case 'custom'
            error('Inneficient to pass all custom trajectories in here.');
        otherwise
            error('Invalid input');
    end

    % Calc Power
    top_power = calc_signal_power(noise_top, dt, total_time);    
    bot_power = calc_signal_power(noise_bot, dt, total_time);  
            
end

function p = calc_signal_power(s, dt, T)
    p = sum(abs(s).^2).*dt./T;

end



        % % Get Noise Trajectories for top and bottom chains
        % [gaussian_noise_top, t_noise_top, white_noise_top] =...
        %     generate_gaussian_noise(t_init, t_final, delta_t, ...,
        %     t_noise_factor, fluctuation_amplitude, characteristic_noise_corr_freq);
        % [gaussian_noise_bot, t_noise_bot, white_noise_bot] =...
        %     generate_gaussian_noise(t_init, t_final, delta_t, ...,
        %     t_noise_factor, fluctuation_amplitude, characteristic_noise_corr_freq);
        % 
        % % Grab a centre segment of the noise trajectories.
        % t_ind_start = floor(length(t_noise_top)/4);
        % t_ind_final = t_ind_start + length(t_vec) - 1;
        %gaussian_noise_top(t_ind_start:t_ind_final);

    