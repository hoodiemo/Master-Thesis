% This file calculates percentiles of the samples after the burn in phase

median_per_row = median(param_chain, 2);
percentile_5 = prctile(param_chain, 5, 2);
percentile_95 = prctile(param_chain, 95, 2);
median_log_posterior = median(logposterior_chain, 1);

save('burned_in_results.mat', 'median_per_row', 'percentile_95', 'percentile_5', 'median_log_posterior');