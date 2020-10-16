%%

% Run EOC stopping and OCBA allocation
[total_samples_bound, total_samples_exact, total_samples_MC, check_times_exact, check_times_bound, check_times_MC] = MonteCarloEfficiency('EOC', 'TS', 'unknown');

%%

figure
hold on
stairs([0, sort(total_samples_bound)], (0:length(total_samples_bound))/length(total_samples_bound), 'm-', 'LineWidth', 1.5);
stairs([0, sort(total_samples_exact)], (0:length(total_samples_exact))/length(total_samples_exact), 'b-', 'LineWidth', 1.5);
stairs([0, sort(total_samples_MC)], (0:length(total_samples_MC))/length(total_samples_MC), 'r-', 'LineWidth', 1.5);
hold off
set(gca, 'FontSize', 14, 'LineWidth', 2)

legend('Bound', 'Exact', 'MC + Exact', 'Location', 'southeast')
legend boxoff

ylim([0,1])

xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold')
ylabel('P(Sample Size $\leq x$)', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold')

%%

figure
stairs([0, sort(check_times_bound)], (0:length(check_times_bound))/length(check_times_bound), 'm-', 'LineWidth', 1.5);
hold on
stairs([0, sort(check_times_exact)], (0:length(check_times_exact))/length(check_times_exact), 'b-', 'LineWidth', 1.5);
stairs([0, sort(check_times_MC)], (0:length(check_times_MC))/length(check_times_MC), 'r-', 'LineWidth', 1.5);
hold off
set(gca, 'FontSize', 14, 'LineWidth', 2)

legend('Bound', 'Exact', 'MC + Exact', 'Location', 'southeast')
legend boxoff

ylim([0,1])

xlabel('$x$ (seconds)', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold')
ylabel('P(Computation Time $\leq x$)', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold')

%%
% %%
% tic;
% W_draw_matrix = zeros(R, k);
% %    for r = 1:R
% %        W_draw_matrix(r,:) = normrnd(sample_means, sqrt((sigma.^2)./n));
% %    end
%    for r = 1:R
%         W_draw_matrix(r,:) = sample_means + sqrt(sample_vars./n).*trnd(n-1);
%     end
% pEOC_MC = mean(max(W_draw_matrix, [], 2) - W_draw_matrix(:,k_brack));
% toc;
% 
% tic;
% W_draw_matrix = zeros(R, k);
% %    for r = 1:R
% %        W_draw_matrix(r,:) = normrnd(sample_means, sqrt((sigma.^2)./n));
% %    end
%    for i = 1:k
%         W_draw_matrix(:, i) = sample_means(i) + sqrt(sample_vars(i)/n(i)).*trnd(n(i)-1, [R, 1]);
%     end
% pEOC_MC = mean(max(W_draw_matrix, [], 2) - W_draw_matrix(:,k_brack));
% toc;

% tic;
% W_draw_matrix = mvnrnd(sample_means, diag(sigma.^2./n), R);
% pEOC_MC = mean(max(W_draw_matrix, [], 2) - W_draw_matrix(:,k_brack));
% toc;