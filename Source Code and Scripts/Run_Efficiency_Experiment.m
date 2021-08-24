%% PGS Experiments

% Run PGS stopping and EA allocation
[total_samples_bound_PGS_EA, total_samples_exact_PGS_EA] = StoppingEfficiency('PGS', 'EA', 'unknown');
frac_savings_PGS_EA = (total_samples_bound_PGS_EA - total_samples_exact_PGS_EA)./total_samples_bound_PGS_EA;

% Run PGS stopping and OCBA allocation
[total_samples_bound_PGS_OCBA, total_samples_exact_PGS_OCBA] = StoppingEfficiency('PGS', 'OCBA-PGS', 'unknown');
frac_savings_PGS_OCBA = (total_samples_bound_PGS_OCBA - total_samples_exact_PGS_OCBA)./total_samples_bound_PGS_OCBA;

% Run PGS stopping and EA allocation
[total_samples_bound_PGS_TS, total_samples_exact_PGS_TS] = StoppingEfficiency('PGS', 'TS', 'unknown');
frac_savings_PGS_TS = (total_samples_bound_PGS_TS - total_samples_exact_PGS_TS)./total_samples_bound_PGS_TS;

%%

figure
hold on
stairs(sort(frac_savings_PGS_EA), (1:length(frac_savings_PGS_EA))/length(frac_savings_PGS_EA), 'm-', 'LineWidth', 1.5);
stairs(sort(frac_savings_PGS_OCBA), (1:length(frac_savings_PGS_OCBA))/length(frac_savings_PGS_OCBA), 'b-', 'LineWidth', 1.5);
stairs(sort(frac_savings_PGS_TS), (1:length(frac_savings_PGS_TS))/length(frac_savings_PGS_TS), 'r-', 'LineWidth', 1.5);
hold off
set(gca, 'FontSize', 14, 'LineWidth', 2)

legend('EA', 'OCBA', 'TS', 'Location', 'southeast')
legend boxoff

xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold')
ylabel('P(Fraction Savings $\leq x$)', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold')

%%

% Run EOC stopping and EA allocation
[total_samples_bound_EOC_EA, total_samples_exact_EOC_EA] = StoppingEfficiency('EOC', 'EA', 'unknown');
frac_savings_EOC_EA = (total_samples_bound_EOC_EA - total_samples_exact_EOC_EA)./total_samples_bound_EOC_EA;

% Run EOC stopping and OCBA allocation
[total_samples_bound_EOC_OCBA, total_samples_exact_EOC_OCBA] = StoppingEfficiency('EOC', 'OCBA-EOC', 'unknown');
frac_savings_EOC_OCBA = (total_samples_bound_EOC_OCBA - total_samples_exact_EOC_OCBA)./total_samples_bound_EOC_OCBA;

% Run EOC stopping and EA allocation
[total_samples_bound_EOC_TS, total_samples_exact_EOC_TS] = StoppingEfficiency('EOC', 'TS', 'unknown');
frac_savings_EOC_TS = (total_samples_bound_EOC_TS - total_samples_exact_EOC_TS)./total_samples_bound_EOC_TS;

%%

figure
hold on
stairs(sort(frac_savings_EOC_EA), (1:length(frac_savings_EOC_EA))/length(frac_savings_EOC_EA), 'm-', 'LineWidth', 1.5);
stairs(sort(frac_savings_EOC_OCBA), (1:length(frac_savings_EOC_OCBA))/length(frac_savings_EOC_OCBA), 'b-', 'LineWidth', 1.5);
stairs(sort(frac_savings_EOC_TS), (1:length(frac_savings_EOC_TS))/length(frac_savings_EOC_TS), 'r-', 'LineWidth', 1.5);
hold off
set(gca, 'FontSize', 14, 'LineWidth', 2)

legend('EA', 'OCBA', 'TS', 'Location', 'southeast')
legend boxoff

xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold')
ylabel('P(Fraction Savings $\leq x$)', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold')