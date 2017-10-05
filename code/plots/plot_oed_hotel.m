
% figure settings
lw = 1;
ms = 6;

n_days = [3, 4, 7, 10, 14];
sim_samples = 96*n_days;
prefix = '../oedLargeHotel/';
kernel = 'sqexp'; % choose 'truong' or 'sqexp'
metrics = {'rmse', 'msll'};
labels = {'RMSE [kW]', 'MSLL'};

figure;
hs = cell(1,2);

for idm = 1:numel(metrics)
    
    metric = metrics{idm};
    
    iter = 1;
    mv = nan(1,numel(sim_samples));
    ig = nan(1,numel(sim_samples));
    unif = nan(1,numel(sim_samples));
    prbs = nan(1,numel(sim_samples));
    
    for ids = 1:numel(sim_samples)
        
        if ids ==2
            continue;
        end
        load(fullfile([prefix 'results_' kernel], 'doe_LargeHotel_MV_2ramped_3input_14day_20171004_0313.mat'));
        mv(iter) = results(ids).test_results.(metric);
        
        load(fullfile([prefix 'results_' kernel], 'doe_LargeHotel_IG_2ramped_3input_14day_20171003_2220.mat'));
        ig(iter) = results(ids).test_results.(metric);
        
        load(fullfile([prefix 'results_' kernel], 'random_LargeHotel_uniform_2ramped_3input_30day_20171004_2030.mat'));
        unif(iter) = results(ids).test_results.(metric);
        
        load(fullfile([prefix 'results_' kernel], 'random_LargeHotel_prbs_2ramped_3input_30day_20171004_2025.mat'));
        prbs(iter) = results(ids).test_results.(metric);
        
        iter = iter + 1;
        
    end
    
    if find(strcmp(metric, {'ae', 'rmse'}))
        mv = mv/1e3;
        ig = ig/1e3;
        unif = unif/1e3;
        prbs = prbs/1e3;
    end
    
    hs{idm} = subplot(1,2,idm); hold on; grid off; box on;
    h1 = plot(1:numel(sim_samples), ig, '-d', 'LineWidth', lw, 'MarkerSize', ms);
    h2 = plot(1:numel(sim_samples), mv, '-d', 'LineWidth', lw, 'MarkerSize', ms);
    h3 = plot(1:numel(sim_samples), unif, '--o', 'LineWidth', lw, 'MarkerSize', ms, 'Color',[0.4940    0.1840    0.5560]);
    h4 = plot(1:numel(sim_samples), prbs, '--o', 'LineWidth', lw, 'MarkerSize', ms, 'Color', [0.4660    0.6740    0.1880]);
    yy = ylabel(labels{idm});
    xlabel('no of days');
    ax = gca;
    ax.XTick = [1, 2, 3, 4];
    ax.XTickLabel = {'3','7','10','14'}';
    xlim([0.75, 4.25])
%     ylim([5, 50])
    
end

hleg = legend([h1, h2, h3, h4], 'IG', 'MV', 'Uniform', 'PRBS');
set(hleg, 'box', 'off', 'Location', 'NorthOutside', 'Orientation', 'Horizontal');
gcf;
hs{1}.Position(4) = hs{2}.Position(4);
matlab2tikz('width', '\fwidth', 'height', '\hwidth', '../../paper/src/figures/casestudy-oed-hotel.tex',...
    'extraaxisoptions',['xlabel style={font=\footnotesize},'...
                       'ylabel style={font=\footnotesize},',...
                       'legend style={font=\footnotesize},',...
                       'ticklabel style={font=\footnotesize},'...
                       'ylabel shift = -5 pt,'...
                       'xlabel shift = -5 pt,'...
                       'legend style={at={(0,1)}}, align=center,']);
