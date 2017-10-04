

sim_samples = 96*[3, 4, 7, 10, 14, 21];

figure; hold on;
for idsim = 1:numel(sim_samples)
    load 'results/random_LargeOffice_uniform_2ramped_3input_28day_20171003_1922.mat';
    plot(idsim, sqrt(results(idsim).test_results.ae), '-db')
    
    load 'results/doe_LargeOffice_IG_2ramped_3input_21day_20171003_1804.mat';
    plot(idsim, sqrt(results(idsim).test_results.ae), '-or')
    
    load('results/random_LargeOffice_prbs_2ramped_3input_28day_20171003_1922.mat');
    plot(idsim, sqrt(results(idsim).test_results.ae), '-*k')
    
end
