%-----------------------%
% Get Effective Gap     %
%-----------------------%
% Returns the effective gap up to some global scaling factor. This is used 
% for instantaneous-adiabatic regime calculations. This function
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


function gap_eff = get_effective_gap_inst_adiabatic(n_chosen, adiabatic_regime_low_bound, ramp_height_vec, ramp_rate_mat, gap_1_min_vec, transition_rate)

    log10_x_mat = log10(gap_1_min_vec.^2./ramp_rate_mat); % Use gap_1 to define the adiabatic regime
    log10_y_mat = log10(transition_rate);
    
    % Use 'fit'
    fit_eqn = '-2*x + a';
    init_guess = -0.2;
    log_y_intercepts = zeros(length(ramp_height_vec),1);
    fit_cell = {};
    
    for ii = 1:length(ramp_height_vec)
    
        log_x_ind_ii = log10_x_mat(:, ii) >= adiabatic_regime_low_bound;
        
        f1 = fit(log10_x_mat(log_x_ind_ii, ii), log10_y_mat(log_x_ind_ii, ii), fit_eqn, 'Start', init_guess);
        log_y_intercepts(ii) = f1.a;
        fit_cell{end+1} = f1;
    
    end  
    
    
    %Note -2 is the hard-coded slope
    delta_X =  (log_y_intercepts(n_chosen) - log_y_intercepts)./(-2); %log_x_inter - log_x_inter(n_chosen);
    gap_eff = gap_1_min_vec.*sqrt(10.^(-delta_X.'));

end

