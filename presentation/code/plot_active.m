
figure; hold on;
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 7 4];
xfill = [x_star; flipdim(x_star,1)]; 
yfill = [mugp+2*stdvgp; flipdim(mugp-2*stdvgp,1)];
h1 = fill(xfill, yfill, [7 7 7]/8, 'EdgeColor', [7 7 7]/8);
alpha(0.75);
h2 = plot(x_star(:,1), mugp, 'k--', 'LineWidth', 2);
h4 = plot(x_star(:,1), y_star, 'b.', 'MarkerSize', 6);
h3 = plot(x(:,1), y, 'rp', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
hleg = legend([h1, h2, h3, h4], '\mu \pm 2\sigma', '\mu', 'selected data', 'given data', ...
        'Location', 'NorthOutside', 'Orientation', 'horizontal');
set(hleg, 'box', 'off');
axis([-5, 5, -3, 7])
xlabel('input x')
ylabel('output y')