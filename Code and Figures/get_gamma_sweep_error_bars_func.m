
function [error_bar_boot, error_bar_sem] = get_gamma_sweep_error_bars_func(N_qp_chain_final_mat, group_size)


    num_trials = size(N_qp_chain_final_mat, 2);
    num_gamma_points = size(N_qp_chain_final_mat, 1);
    num_groups = floor(num_trials/group_size);
    error_bar_sem = zeros(num_gamma_points, 1);
    error_bar_boot = zeros(num_gamma_points, 1);
    num_resamp = 10000;
    for jj = 1:num_gamma_points
    
        means_vec = zeros(num_groups, 1);
    
        for ii = 1:num_groups
            curr_start =group_size*(ii-1)+1;
            curr_group = N_qp_chain_final_mat(jj, curr_start:(curr_start + group_size - 1));
            means_vec(ii) = mean(curr_group);
        end
        error_bar_sem(jj) = std(means_vec)/sqrt(num_groups);
        error_bar_boot(jj) = bootstrap_resampling(means_vec, num_resamp);
    end
end


%% Local Functions

function s_boot = bootstrap_resampling(sample_data, num_resamp)

    boot_means = zeros(num_resamp,1);
    num_data = length(sample_data); 
    for ii = 1:num_resamp
        resample = sample_data(randi(num_data, num_data, 1));
        boot_means(ii) = mean(resample);    
    end
    s_boot = std(boot_means); 

end
