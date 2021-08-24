% Plot pdf for Weibull(a = 4 (scale), b = 2 (shape)) random variable

x = 0:0.01:15;
xpdf = wblpdf(x, 4, 2);

figure
plot(-x, xpdf, 'k-', 'LineWidth', 1.5);
box off

set(gca, 'FontSize', 14, 'LineWidth', 2)

xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold')
ylabel('$f(x)$', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold')
title('-Weibull($a=4$, $b=2$) pdf', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold')

%%
% Plot pdf for Weibull(a = 1.5 (scale), b = 2 (shape)) random variable

x = 0:0.01:15;
xpdf = wblpdf(x, 1.5, 2);

figure
plot(-x, xpdf, 'k-', 'LineWidth', 1.5);
box off

set(gca, 'FontSize', 14, 'LineWidth', 2)

xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold')
ylabel('$f(x)$', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold')
title('-Weibull($a=1.5$, $b=2$) pdf', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold')

%%
% Plot pdf for Weibull(a = 1 (scale), b = 2 (shape)) random variable

x = 0:0.01:15;
xpdf = wblpdf(x, 1, 2);

figure
plot(-x, xpdf, 'k-', 'LineWidth', 1.5);
box off

set(gca, 'FontSize', 14, 'LineWidth', 2)

xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold')
ylabel('$f(x)$', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold')
title('-Weibull($a=1.5$, $b=2$) pdf', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold')


%%
% Plot pdf for Chi-Squared(4) random variable

y = 0:0.01:10;
ypdf = chi2pdf(y, 4);

figure
plot(y, ypdf, 'k-', 'LineWidth', 1.5);
box off

set(gca, 'FontSize', 14, 'LineWidth', 2)

xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold')
ylabel('$f(x)$', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold')
title('$\chi^2$($\nu=4$) pdf', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold')

%%
% Plot pdf for Lognormal(mu = 0, sigma = 1.5) random variable

z = 0:0.01:10;
zpdf = lognpdf(z, 0, 1.5);

figure
plot(-z, zpdf, 'k-', 'LineWidth', 1.5);
box off

set(gca, 'FontSize', 14, 'LineWidth', 2)

xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold')
ylabel('$f(x)$', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold')
title('Lognormal($\mu=0$, $\sigma = 1.5$) pdf', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold')
