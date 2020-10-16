% Plot the Bonferroni and Slepian bounds for posterior PGS

clear
clc

% Confidence parameters to test
alpha_list = [0.01, 0.05, 0.1]; % 1-alpha = 0.90, 0.95, and 0.99

% Number of alternatives to test
k_list = [2:2:10, 20:10:100, 200:100:1000];

% Initialize
Bonf_pPGS = zeros(length(alpha_list), length(k_list));
Slep_pPGS = zeros(length(alpha_list), length(k_list));

% Test different values of alpha
for j = 1:length(alpha_list)
    
    alpha = alpha_list(j);
    
    % Test different values of k
    for i = 1:length(k_list)
        
        k = k_list(i);

        fprintf('Testing alpha = %.2f and k = %d.\n', alpha, k)

        % Find slippage configuration of posterior means for which the
        % posterior PGS of Alternative (k) is exactly alpha.
        % Posterior variances are assumed to all be one.
        h = calcBechhoferh(k, alpha);
        
        % Calculate posterior PGS Bonferroni bound for the best alternative
        Bonf_pPGS(j,i) = 1 - (k-1)*normcdf(-h);
        
        % Calculate posteiror PGS Slepian bound for the best alternative
        Slep_pPGS(j,i) = normcdf(h)^(k-1);
        
    end
end

%%

% figure
% hold on
% h1 = plot(k_list, Bonf_pPGS(1,:), 'b:', 'LineWidth', 1.5);
% plot(k_list, Bonf_pPGS(2,:), 'b:', 'LineWidth', 1.5);
% plot(k_list, Bonf_pPGS(3,:), 'b:', 'LineWidth', 1.5);
% h2 = plot(k_list, Slep_pPGS(1,:), 'r-', 'LineWidth', 1.5);
% plot(k_list, Slep_pPGS(2,:), 'r-', 'LineWidth', 1.5);
% plot(k_list, Slep_pPGS(3,:), 'r-', 'LineWidth', 1.5);
% hold off

figure
h1 = semilogx(k_list, Bonf_pPGS(1,:), 'b:', 'LineWidth', 1.5);
hold on
semilogx(k_list, Bonf_pPGS(2,:), 'b:', 'LineWidth', 1.5);
semilogx(k_list, Bonf_pPGS(3,:), 'b:', 'LineWidth', 1.5);
h2 = semilogx(k_list, Slep_pPGS(1,:), 'r-', 'LineWidth', 1.5);
semilogx(k_list, Slep_pPGS(2,:), 'r-', 'LineWidth', 1.5);
semilogx(k_list, Slep_pPGS(3,:), 'r-', 'LineWidth', 1.5);
hold off
box off

set(gca, 'FontSize', 14, 'LineWidth', 2)

legend([h1,h2],{' pPGS^{Bonf}_{(k)}', ' pPGS^{Slep}_{(k)}'}, 'Location', 'southwest')
legend boxoff

%title('pPGS Bounds in Slippage Configuration', 'FontSize', 16, 'FontWeight', 'bold')
xlabel('Number of Alternatives (k)', 'FontSize', 16, 'FontWeight', 'bold')
ylabel('pPGS^{Bonf}_{(k)} and pPGS^{Slep}_{(k)}', 'FontSize', 16, 'FontWeight', 'bold')