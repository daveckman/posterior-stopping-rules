% Plot the Bonferroni and Slepian bounds for posterior EOC
% Assumes known variances --> posterior is MVN distribution

clear
clc

% EOC parameters to test
beta_list = [0.05, 0.1, 0.25];

% Number of alternatives to test
k_list = [2:2:10, 20:10:100, 200:100:1000];

% Initialize
Bonf_pEOC = zeros(length(beta_list), length(k_list));
Slep_pEOC = zeros(length(beta_list), length(k_list));

DeltaSC = 0; % Initial value for root search

% Test different values of beta
for j = 1:length(beta_list)
    
    beta = beta_list(j);
    
    % Test different values of k
    for i = 1:length(k_list)
        
        k = k_list(i);
        
        fprintf('Testing beta = %.2f and k = %d.\n', beta, k)
        
        % These functions calculate the expected value of the maximum of (k-1)
        % N(0, 1) rvs and a N(Delta, 1) rv - i.e., a slippage configuration.
        prodfun = @(z, Delta) z.*normcdf(z).^(k-2).*((k-1)*normcdf(z-Delta).*normpdf(z) + normcdf(z).*normpdf(z-Delta));
        Emaxintegral = @(Delta) integral(@(z)prodfun(z, Delta), -Inf, Inf);

        % Find the value of Delta that describes a slippage configuration of posterior means
        % for which the posterior EOC of Alternative (k) is exactly beta.
        % Posterior variances are assumed to all be one.

        % Formula pEOC = Emax - Delta 
        pEOCrootfun = @(Delta) Emaxintegral(Delta) - Delta - beta;

        % Root search starting near the old value of DeltaSC
        DeltaSC = fzero(pEOCrootfun, DeltaSC); 

        % Monte Carlo check to see that the integral calculation is working.
        % A = mvnrnd([zeros(1,k-1), DeltaSC], eye(k), 10000);
        % mean(max(A,[],2)) - DeltaSC - beta; % Should be close to zero

        % Calculate posterior EOC Bonferroni bound for the best alternative
        Bonf_pEOC(j,i) = (k-1)*sqrt(2)*(normpdf(DeltaSC/sqrt(2)) - (DeltaSC/sqrt(2))*normcdf(-DeltaSC/sqrt(2)));
        % Bonf_pEOC(j,i) = CalcPostQty('EOC-Bonf', k, [zeros(1, k-1), DeltaSC], ones(1, k), zeros(1, k), 'known', 0); % Check

        % Calculate posterior EOC Slepian bound for the best alternative
        SlepPBS = @(delta) 1 - normcdf((delta + DeltaSC)/sqrt(2)).^(k-1);
        Slep_pEOC(j,i) = integral(@(delta)SlepPBS(delta), 0, Inf);      
        % Slep_pEOC(j,i) = CalcPostQty('EOC-Slep', k, [zeros(1, k-1), DeltaSC], ones(1, k), zeros(1, k), 'known', 0); % Check

    end
end

%%

% figure
% hold on
% h1 = plot(k_list, Bonf_pEOC(1,:), 'b:', 'LineWidth', 1.5);
% plot(k_list, Bonf_pEOC(2,:), 'b:', 'LineWidth', 1.5);
% plot(k_list, Bonf_pEOC(3,:), 'b:', 'LineWidth', 1.5);
% h2 = plot(k_list, Slep_pEOC(1,:), 'r-', 'LineWidth', 1.5);
% plot(k_list, Slep_pEOC(2,:), 'r-', 'LineWidth', 1.5);
% plot(k_list, Slep_pEOC(3,:), 'r-', 'LineWidth', 1.5);
% hold off

figure
h1 = semilogx(k_list, Bonf_pEOC(1,:), 'b:', 'LineWidth', 1.5);
hold on
semilogx(k_list, Bonf_pEOC(2,:), 'b:', 'LineWidth', 1.5);
semilogx(k_list, Bonf_pEOC(3,:), 'b:', 'LineWidth', 1.5);
h2 = semilogx(k_list, Slep_pEOC(1,:), 'r-', 'LineWidth', 1.5);
semilogx(k_list, Slep_pEOC(2,:), 'r-', 'LineWidth', 1.5);
semilogx(k_list, Slep_pEOC(3,:), 'r-', 'LineWidth', 1.5);
hold off
box off

set(gca, 'FontSize', 14, 'LineWidth', 2)

legend([h1,h2],{' pEOC^{Bonf}_{(k)}', ' pEOC^{Slep}_{(k)}'}, 'Location', 'northwest')
legend boxoff

%title('pEOC Bounds in Slippage Configuration', 'FontSize', 16, 'FontWeight', 'bold')
xlabel('Number of Alternatives (k)', 'FontSize', 16, 'FontWeight', 'bold')
ylabel('pEOC^{Bonf}_{(k)} and pEOC^{Slep}_{(k)}', 'FontSize', 16, 'FontWeight', 'bold')