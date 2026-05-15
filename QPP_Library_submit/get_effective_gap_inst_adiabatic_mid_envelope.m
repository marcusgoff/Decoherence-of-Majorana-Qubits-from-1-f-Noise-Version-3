%-----------------------%
% Get Effective Gap     %
%-----------------------%
% Returns the effective gap up to some global scaling factor. This is used 
% for instantaneous-adiabatic regime calculations. 
%
% Here I used the "middle-envelope" in the leakage, which is half way
% between the top and bottom envelope of the leakage oscillations. This is
% particularly important for low_mu simulations. 
%
% Note that here, unlike in my earlier code, I use x_1 =
% ramp_rate/gap_1_min^2. This is the reprocal of what I used in my earlier
% code. 
%
% This function
% fits a straight line with gradient "-2" to the adiabatic regime and then
% finds the effective gaps which causes all the fits to overlap. Note that 
% multiplying the vector gap_eff by any constant will just translate all
% curves along the x-axis by a constant. To deal with this redundancy, this 
% code assumes that "g_1" for n_chosen mu_final is "correct" and aligns all
% fits to this curve's fit. 
%
% I will need to find a better way to deal with the freedom of choice of
% this scaling factor. 
%
% Note that this function uses for the x-ais 
%
% INPUTS:
%       transition_rate - L_p, the parity leakage


function gap_eff = get_effective_gap_inst_adiabatic_mid_envelope(n_chosen,...
    adiabatic_regime_up_bound, ramp_height_vec, ramp_rate_mat, gap_1_min_vec, ... 
    transition_rate, plot_calcs)

    %log10_x_mat = log10(gap_1_min_vec.^2./ramp_rate_mat); % Use gap_1 to define the adiabatic regime
    x_mat = ramp_rate_mat./gap_1_min_vec.^2;
    y_mat = transition_rate;

    %log10_x_mat = log10(x_mat); % Use gap_1 to define the adiabatic regime
    %log10_y_mat = log10(y_mat);
    
    % Use 'fit'
    fit_eqn = '2*x + a'; %'-2*x + a';
    init_guess = 0.5;
    log_y_intercepts = zeros(length(ramp_height_vec),1);
    fit_cell = {};
    
    x_cell = {};
    min_cell = {};
    max_cell = {};
    top_env_cell = {};
    bot_env_cell = {};
    mid_env_cell = {};

   if plot_calcs
        fig1 = figure(); 
        loglog((x_mat), log10(y_mat));
        xlabel('v/g_1^2');
        ylabel('Leakage');
        axis tight;
    end


    for ii = 1:length(ramp_height_vec)
    
        log_x_ind_ii = log10(x_mat(:, ii)) <= adiabatic_regime_up_bound;
        
        x = x_mat(log_x_ind_ii, ii);
        y = y_mat(log_x_ind_ii, ii);
        [y_mid, top_envelope, bot_envelope, min_mat, max_mat] =  fit_top_and_bottom(x, y);


        %f1 = fit(log10_x_mat(log_x_ind_ii, ii), log10_y_mat(log_x_ind_ii, ii), fit_eqn, 'Start', init_guess);
        f1 = fit(log10(x) , log10(y_mid), fit_eqn, 'Start', init_guess);        
        log_y_intercepts(ii) = f1.a;
        fit_cell{end+1} = f1;

        x_cell{end+1} = x; min_cell{end+1} = min_mat; max_cell{end+1} = max_mat;
        top_env_cell{end+1} = top_envelope; bot_env_cell{end+1} = bot_envelope;
        mid_env_cell{end+1} = y_mid; % Consider outputting all of these
    
        if plot_calcs
            title_str = sprintf('\\Delta mu = %.3g', ramp_height_vec(ii));

            figure(); 
            subplot(1,2,1);
            plot(x, y); hold on;
            plot(min_mat(:,1), min_mat(:,2), 'o'); 
            plot(max_mat(:,1), max_mat(:,2), 'x'); 
            plot(x, y_mid, '-', 'LineWidth', 1.2); plot(x, top_envelope, '-.');
            plot(x, bot_envelope, '-.');  
            loglog(x, 10.^(f1.a).*x.^2, 'r-.', 'LineWidth', 1.2);
            xlabel('v/g_1^2');
            title(title_str); 
     
            subplot(1,2,2);            
            loglog(x, y); hold on;
            loglog(min_mat(:,1), min_mat(:,2), 'o'); 
            loglog(max_mat(:,1), max_mat(:,2), 'x'); 
            loglog(x, y_mid, '-', 'LineWidth', 1.2); plot(x, top_envelope, '-.'); 
            loglog(x, bot_envelope, '-.');   
            loglog(x, 10.^(f1.a).*x.^2, 'r-.', 'LineWidth', 1.2);
            xlabel('v/g_1^2');
            title(title_str);         
        end
    end  
    % TO DO:
    %
    % Why is my constrained fit not lining up with my y_mid line?

    %Note -2 is the hard-coded slope
    delta_X =  (log_y_intercepts(n_chosen) - log_y_intercepts)./(2); %log_x_inter - log_x_inter(n_chosen);
    %gap_eff = gap_1_min_vec.*sqrt(10.^(-delta_X.'));
    gap_eff = gap_1_min_vec(n_chosen).*10.^(delta_X.'/2);


    if plot_calcs
        figure();
        subplot(1,2,1);
        ii = 1; plot(log10(x_cell{ii}./gap_eff(ii).^2), log10(mid_env_cell{ii})); hold on;
        ii = 2; plot(log10(x_cell{ii}./gap_eff(ii).^2), log10(mid_env_cell{ii})); hold on;
        ii = 3; plot(log10(x_cell{ii}./gap_eff(ii).^2), log10(mid_env_cell{ii})); hold on;
        title('Translated Mid Envelopes');
        xlabel('v/g_{eff}^2');

        subplot(1,2,2);
        ii = 1;  log_x = log10(x_cell{ii});
        loglog(x_cell{ii}./gap_eff(ii).^2, 10.^(2*(log_x) + log_y_intercepts(ii))); hold on;

        ii = 2;  log_x = log10(x_cell{ii});
        loglog(x_cell{ii}./gap_eff(ii).^2, 10.^(2*(log_x) + log_y_intercepts(ii)), '-.'); hold on;

        ii = 3;  log_x = log10(x_cell{ii});
        loglog(x_cell{ii}./gap_eff(ii).^2, 10.^(2*(log_x) + log_y_intercepts(ii)), ':'); hold on;  
        title('Translated Fits');
        xlabel('v/g_{eff}^2');
    end

end



%% Local Functions

function [y_mid, top_envelope, bot_envelope, min_mat, max_mat] =  fit_top_and_bottom(x, y)

    min_ind = islocalmin(y);
    max_ind = islocalmax(y);
    
    y_max = y(max_ind);
    x_max = x(max_ind);
    max_mat = [x_max(:), y_max(:)];
    
    y_min = y(min_ind);
    x_min = x(min_ind); 
    min_mat = [x_min(:), y_min(:)];
    
    fit_top = polyfit(log10(x_max), log10(y_max), 1);
    fit_bot = polyfit(log10(x_min), log10(y_min), 1);
    
    top_envelope = (10^fit_top(2).*x.^fit_top(1));
    bot_envelope = (10^fit_bot(2).*x.^fit_bot(1));
    
    y_mid = (bot_envelope + top_envelope)./2;

end



        % 
        % figure();
        % ii = 1;  log_x = log10(x_cell{ii});
        % plot(log_x, 2*(log_x + delta_X(ii)) + log_y_intercepts(ii)); hold on;
        % ii = 2;  log_x = log10(x_cell{ii});
        % plot(log_x, 2*(log_x + delta_X(ii)) + log_y_intercepts(ii), '-.'); hold on;
        % ii = 3;  log_x = log10(x_cell{ii});
        % plot(log_x, 2*(log_x + delta_X(ii)) + log_y_intercepts(ii), ":"); 

        % figure();
        % ii = 1;  log_x = log10(x_cell{ii});
        % plot(log_x - delta_X(ii), 2*(log_x) + log_y_intercepts(ii)); hold on;
        % 
        % ii = 2;  log_x = log10(x_cell{ii});
        % plot(log_x - delta_X(ii), 2*(log_x) + log_y_intercepts(ii), '-.'); hold on;
        % 
        % ii = 3;  log_x = log10(x_cell{ii});
        % plot(log_x - delta_X(ii), 2*(log_x) + log_y_intercepts(ii), ":");    
