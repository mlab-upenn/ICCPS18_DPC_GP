function plot_results
plot_for_postfixes({'3days_final', '7days_final', '14days_final'}, 'SE');
plot_for_postfixes({'3days_truongkernel', '7days_truongkernel', '14days_truongkernel'}, 'Truong''s kernel');
end

function plot_for_postfixes(file_postfixes, txt)

file_prefixes = {'doe_sampling_noreset_IG_3input_ahead00',... % 'doe_sampling_noreset_IG_3input_ahead00_prior',...
    'random_sampling_uniform_3input_ahead00',...
    'random_sampling_prbs_3input_ahead00'};

result_rmse = nan(numel(file_postfixes), numel(file_prefixes));
result_msll = result_rmse;

for k = 1:numel(file_postfixes)
    for p = 1:numel(file_prefixes)
        filename = [file_prefixes{p} '_' file_postfixes{k}];
        result = load(filename, 'validation_result');
        result_rmse(k, p) = result.validation_result.loss_rmse;
        result_msll(k, p) = result.validation_result.loss_msll;
    end
end

clear result

figure
bar(result_rmse);
title([txt ' - RMSE']);
legend('OED', 'uniform', 'prbs'); % 'OED-prior', 

figure
bar(result_msll);
title([txt ' - MSLL']);
legend('OED', 'uniform', 'prbs'); %'OED-prior', 
end