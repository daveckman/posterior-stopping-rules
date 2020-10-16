% Produce plots for Monte Carlo experiments

post_quantity = 'EOC'; %'EOC'
rpi_mode = 2;
M = 100;

% Load data
load(['MonteCarlo_TS_',post_quantity,'_Exact1_rpi_mode=',num2str(rpi_mode),'_M=',num2str(M),'.mat'])
%load(['MonteCarlo_TS_',post_quantity,'_Exact2_rpi_mode=',num2str(rpi_mode),'_M=',num2str(M),'.mat'])

%% Sample size ecdf

figure
hold on
stairs([0, sort(total_samples_Bonf_bound)], (0:length(total_samples_Bonf_bound))/length(total_samples_Bonf_bound), 'g-', 'LineWidth', 2);
stairs([0, sort(total_samples_Slep_bound)], (0:length(total_samples_Slep_bound))/length(total_samples_Slep_bound), 'm-', 'LineWidth', 2);
stairs([0, sort(total_samples_exact)], (0:length(total_samples_exact))/length(total_samples_exact), 'b-', 'LineWidth', 2);
stairs([0, sort(total_samples_MC)], (0:length(total_samples_MC))/length(total_samples_MC), 'r-', 'LineWidth', 2);
hold off
set(gca, 'FontSize', 14, 'LineWidth', 2)

legend('$N_b^{\mathrm{Bonf}}$', '$N_b^{\mathrm{Slep}}$', '$N_e$', '$N_{mc}$', 'Location', 'southeast', 'Interpreter', 'latex')
legend boxoff

ylim([0,1])

xlabel('Sample Size ($n$)', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold')
ylabel('$$P(N \leq n)$', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold')

% Report errors for empirical cdfs
sample_size_struct = {total_samples_Bonf_bound, total_samples_Slep_bound, total_samples_exact, total_samples_MC};

disc = linspace(0, max([total_samples_Bonf_bound, total_samples_Slep_bound, total_samples_exact, total_samples_MC]), 500); % discretize total sample sizes (s)
stderrors = zeros(size(sample_size_struct,2), length(disc));
for struct_index = 1:size(sample_size_struct,2)
    sample_size = sample_size_struct{struct_index};
    for col_index = 1:length(disc)
        n = disc(col_index);
        proportions = (sample_size < n);
        stderrors(struct_index, col_index) = 1.96*std(proportions)/sqrt(M);
    end
    fprintf('Max half width (1.96*std error) = %f.\n', max(stderrors(struct_index,:)))
end

%% Compuational time ecdf

figure
hold on
stairs([0, sort(check_times_Bonf_bound)], (0:length(check_times_Bonf_bound))/length(check_times_Bonf_bound), 'g-', 'LineWidth', 2);
stairs([0, sort(check_times_Slep_bound)], (0:length(check_times_Slep_bound))/length(check_times_Slep_bound), 'm-', 'LineWidth', 2);
stairs([0, sort(check_times_exact)], (0:length(check_times_exact))/length(check_times_exact), 'b-', 'LineWidth', 2);
stairs([0, sort(check_times_MC)], (0:length(check_times_MC))/length(check_times_MC), 'r-', 'LineWidth', 2);
hold off
set(gca, 'FontSize', 14, 'LineWidth', 2)

legend('$T_b^{\mathrm{Bonf}}$', '$T_b^{\mathrm{Slep}}$', '$T_e$', '$T_{mc}$', 'Location', 'southeast', 'Interpreter', 'latex')
legend boxoff

box off

ylim([0,1])

xlabel('Computation Time ($t$) (seconds)', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold')
ylabel('$P(T \leq t)$', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold')

% Report errors for empirical cdfs
comp_time_struct = {check_times_Bonf_bound, check_times_Slep_bound, check_times_exact, check_times_MC};

disc = linspace(0, max([check_times_Bonf_bound, check_times_Slep_bound, check_times_exact, check_times_MC]), 500); % discretize total sample sizes (s)
stderrors = zeros(size(comp_time_struct,2), length(disc));
for struct_index = 1:size(comp_time_struct,2)
    comp_time = comp_time_struct{struct_index};
    for col_index = 1:length(disc)
        t = disc(col_index);
        proportions = (comp_time < t);
        stderrors(struct_index, col_index) = 1.96*std(proportions)/sqrt(M);
    end
    fprintf('Max half width (1.96*std error) = %f.\n', max(stderrors(struct_index,:)))
end

%% Plot histogram of pEOC at time of pure MC stopping

figure
hold on
%histogram(pure_MC_stop_post_qty, 10)
%bins = [0.4:0.025:0.55];
%bin_centers = [0.4:0.025:0.525] + 0.0125;
bins = [0.4:0.02:0.54];
bin_centers = [0.41:0.02:0.53];
%bins = [0.4:0.01:0.53];
%bin_centers = [0.405:0.01:0.525];
[bin_counts, ~] = histcounts(pure_MC_stop_post_qty, bins);
bar(bin_centers, bin_counts, 'LineWidth', 1)
line([0.5, 0.5], [0, 40], 'Color', 'black', 'LineStyle', '--', 'LineWidth', 2);
hold off
set(gca, 'FontSize', 14, 'LineWidth', 2)

xlabel('$\mathrm{pPEOC}_{d}$', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold')
ylabel('Frequency')

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