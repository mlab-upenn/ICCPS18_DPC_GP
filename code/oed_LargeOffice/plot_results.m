
sim_samples = 96*[3, 4, 7, 10, 14, 21];
kernel = 'sqexp'; % choose 'truong' or 'sqexp'

% * ae     ... the mean absolute error
% * se     ... the mean squared error 
% * rmse   ... the root mean squared error 
% * lpd    ... the log-predicted loss
% * mrse   ... the mean relative square error 
% * smse   ... the standardized mean squared error
% * msll   ... the mean standardized log loss
metric = 'ae'; 

figure; hold on;
for ids = 1:numel(sim_samples)

    if ids ==2
        continue;
    end
    load(fullfile(['results_' kernel], 'doe_LargeOffice_MV_2ramped_3input_28day_20171003_1817.mat'));
    h1 = plot(ids, results(ids).test_results.(metric), '-db');
    
    load(fullfile(['results_' kernel], 'doe_LargeOffice_IG_2ramped_3input_21day_20171003_1804.mat'));
    h2 = plot(ids, results(ids).test_results.(metric), '-or');
    
    load(fullfile(['results_' kernel], 'random_LargeOffice_uniform_2ramped_3input_28day_20171003_1922.mat'));
    h3 = plot(ids, results(ids).test_results.(metric), '-*k');
    
    load(fullfile(['results_' kernel], 'random_LargeOffice_prbs_2ramped_3input_28day_20171004_2034.mat'));
    h4 = plot(ids, results(ids).test_results.(metric), '-*g');
    
end

legend([h1, h2, h3, h4], 'MV', 'IG', 'uniform', 'prbs')