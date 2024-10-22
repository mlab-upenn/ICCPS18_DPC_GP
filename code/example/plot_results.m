hold('off');

fill([x_star(:,1); flipud(x_star(:,1))], ...
     [       f_star_mean + 2 * sqrt(f_star_variance); ...
      flipud(f_star_mean - 2 * sqrt(f_star_variance))], ...
     [0.9, 0.9, 1], ...
     'edgecolor', 'none');

hold('on');

plot(x_star(:,1), y_star, 'r.');
plot(x(:,1), y, 'k+');
plot(x_star(:,1), f_star_mean, '-', ...
     'color', [0, 0, 0.8]);

axis([-5, 5, -4, 6]);
set(gca, 'tickdir', 'out', ...
         'box',     'off');
