
% figure settings
lw = 1;
ms = 6;

n_days = [3, 4, 7, 10, 14, 21];
sim_samples = 96*n_days;
prefix = '../oedLargeOffice/';
kernel = 'sqexp'; % choose 'truong' or 'sqexp'
metrics = {'rmse', 'smse'};
labels = {'RMSE [kW]', '1-SMSE [%]'};

figure;
title('Office')
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
        
        if idm==1
            load(fullfile([prefix 'results_' kernel], 'doe_LargeOffice_MV_2ramped_3input_28day_20171003_1817.mat'));
            mv(iter) = results(ids).test_results.(metric);
            
            load(fullfile([prefix 'results_' kernel], 'doe_LargeOffice_IG_2ramped_3input_21day_20171003_1804.mat'));
            ig(iter) = results(ids).test_results.(metric);
            
            load(fullfile([prefix 'results_' kernel], 'random_LargeOffice_uniform_2ramped_3input_28day_20171003_1922.mat'));
            unif(iter) = results(ids).test_results.(metric);
            
            load(fullfile([prefix 'results_' kernel], 'random_LargeOffice_prbs_2ramped_3input_28day_20171004_2034.mat'));
            prbs(iter) = results(ids).test_results.(metric);
        else
            load(fullfile([prefix 'results_' kernel], 'doe_LargeOffice_MV_2ramped_3input_28day_20171003_1817.mat'));
            mv(iter) = 100*(1-results(ids).test_results.(metric));
            
            load(fullfile([prefix 'results_' kernel], 'doe_LargeOffice_IG_2ramped_3input_21day_20171003_1804.mat'));
            ig(iter) = 100*(1-results(ids).test_results.(metric));
            
            load(fullfile([prefix 'results_' kernel], 'random_LargeOffice_uniform_2ramped_3input_28day_20171003_1922.mat'));
            unif(iter) = 100*(1-results(ids).test_results.(metric));
            
            load(fullfile([prefix 'results_' kernel], 'random_LargeOffice_prbs_2ramped_3input_28day_20171004_2034.mat'));
            prbs(iter) = 100*(1-results(ids).test_results.(metric));
        end
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
    ylabel(labels{idm});
    xlabel('no of hrs');
    ax = gca;
    ax.XTick = [1, 2, 3, 4, 5];
%     ax.XTickLabel = {'3','7','10','14'}';    
    ax.XTickLabel = {'72','168','240','336', '504'}';
    xlim([0.75, 5.25])
    hold off
    title('OFFICE')
    if idm==2
        ax = gca;
        ax.YTick = [50, 60, 70, 80, 90];
    end
end

hleg = legend([h1, h2, h3, h4], 'IG', 'MV', 'Uniform', 'PRBS');
set(hleg, 'box', 'off', 'Location', 'NorthOutside', 'Orientation', 'Horizontal');
gcf;
hs{1}.Position(4) = hs{2}.Position(4);
matlab2tikz('width', '\fwidth', 'height', '\hwidth', '../../paper/src/figures/casestudy-oed-office.tex',...
    'extraaxisoptions',['xlabel style={font=\footnotesize},'...
                       'ylabel style={font=\footnotesize},',...
                       'legend style={font=\footnotesize},',...
                       'ticklabel style={font=\footnotesize},'...
                       'ylabel shift = -5 pt,'...
                       'xlabel shift = -5 pt,'...
                       'title style={font=\normalsize},'...
                       'title style={yshift=1.5ex},'...
                       'legend style={at={(0,1)}}, align=center']);

