function plot_oed
%plot_for_postfixes({'3days_final', '7days_final', '14days_final'}, '', 'SE');
%plot_for_postfixes({'3days_truongkernel', '7days_truongkernel', '14days_truongkernel'}, '', 'Truong''s kernel');
plot_for_postfixes({'3days_truongkernel', '4days_truongkernel', '7days_truongkernel', '14days_truongkernel'}, '2ramped_', 'RAMPED - Truong''s kernel');
end

function plot_for_postfixes(file_postfixes, file_midfix, txt)

prefix = '../training/feature_selection/';
file_prefixes = {[prefix 'doe_sampling_noreset_IG_' file_midfix '3input_ahead00'],... % 'doe_sampling_noreset_IG_3input_ahead00_prior',...
%     [prefix 'doe_sampling_noreset_IG_' file_midfix '3input_final_ahead00'],... 
    [prefix 'random_sampling_uniform_' file_midfix '3input_ahead00']}; %,...
    %['random_sampling_prbs_' file_midfix '3input_ahead00']};

result_rmse = nan(numel(file_postfixes), numel(file_prefixes));
result_msll = result_rmse;
result_ae = result_rmse;
result_lpd = result_rmse;

for k = 1:numel(file_postfixes)
    for p = 1:numel(file_prefixes)
        filename = [file_prefixes{p} '_' file_postfixes{k}];
        result = load(filename, 'validation_result');
        result_rmse(k, p) = result.validation_result.loss_rmse;
        result_msll(k, p) = result.validation_result.loss_msll;
        result_ae(k, p) = result.validation_result.loss_ae;
        result_lpd(k, p) = result.validation_result.loss_lpd;
    end
end

clear result

result_ae = result_ae/1e3;
result_rmse = result_rmse/1e3;

figure; grid off; hold on;
linwidth = 1;
markersize = 6;
h1 = plot(1:numel(result_rmse(:,1)), result_rmse(:,1), '-ok', 'LineWidth', linwidth, 'MarkerSize', markersize);
h2 = plot(1:numel(result_ae(:,1)), result_ae(:,1), '-dk', 'LineWidth', linwidth, 'MarkerSize', markersize);
h3 = plot(1:numel(result_rmse(:,2)), result_rmse(:,2), '--ok', 'LineWidth', linwidth, 'MarkerSize', markersize);
h4 = plot(1:numel(result_ae(:,2)), result_ae(:,2), '--dk', 'LineWidth', linwidth, 'MarkerSize', markersize);

hleg = legend([h1, h2, h3, h4], 'OED: RMSE', 'OED: AE', 'Random: RMSE', 'Random: AE');
set(hleg, 'box', 'off', 'Location', 'NorthEast');
ax = gca;
ax.XTick = [1, 2, 3, 4];
ax.XTickLabel = {'3 days','4 days','7 days','15 days'}';
xlim([0.75, 4.25])
ylim([5, 50])
ylabel('error [kW]')

matlab2tikz('width', '\fwidth', 'height', '\hwidth', '../../paper/src/figures/oed-acc.tex');

% figure
% bar(result_rmse);
% title([txt ' - RMSE']);
% legend('OED', 'OED-final', 'uniform'); % 'OED-prior', , 'prbs'
% 
% figure
% bar(result_msll);
% title([txt ' - MSLL']);
% legend('OED', 'OED-final', 'uniform'); %'OED-prior', , 'prbs'
% 
% figure
% bar(result_ae);
% title([txt ' - AE']);
% legend('OED', 'OED-final', 'uniform'); %'OED-prior', , 'prbs'
% 
% figure
% bar(result_lpd);
% title([txt ' - LPD']);
% legend('OED', 'OED-final', 'uniform'); %'OED-prior', , 'prbs'
end