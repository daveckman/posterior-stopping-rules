clear
clc

% Experimental setup
post_obj = 'PGS'; % string ('PGS' or 'EOC')
known_var = 'unknown'; % unknown variances
k = 100; % number of systems
n0 = 5; % common initial sample size
M = 100; % number of macroreplications
Q = 50; % number of splits

% Experiments return matrices of size M rows x Q cols
n_obs = M*Q; % common number of observations

%% Run experiments

% Run stopping and EA allocation
[total_samples_slep_bound_EA, total_samples_bonf_bound_EA, total_samples_exact_EA] = StoppingEfficiencyDoubleSplit(post_obj, 'EA', known_var, k, n0, M, Q);
frac_savings_slep_EA = (total_samples_slep_bound_EA - total_samples_exact_EA)./total_samples_slep_bound_EA;
frac_savings_bonf_EA = (total_samples_bonf_bound_EA - total_samples_exact_EA)./total_samples_bonf_bound_EA;

% Run stopping and OCBA allocation
[total_samples_slep_bound_OCBA, total_samples_bonf_bound_OCBA, total_samples_exact_OCBA] = StoppingEfficiencyDoubleSplit(post_obj, ['OCBA-',post_obj], known_var, k, n0, M, Q);
frac_savings_slep_OCBA = (total_samples_slep_bound_OCBA - total_samples_exact_OCBA)./total_samples_slep_bound_OCBA;
frac_savings_bonf_OCBA = (total_samples_bonf_bound_OCBA - total_samples_exact_OCBA)./total_samples_bonf_bound_OCBA;

% Run stopping and EA allocation
[total_samples_slep_bound_TS, total_samples_bonf_bound_TS, total_samples_exact_TS] = StoppingEfficiencyDoubleSplit(post_obj, 'TS', known_var, k, n0, M, Q);
frac_savings_slep_TS = (total_samples_slep_bound_TS - total_samples_exact_TS)./total_samples_slep_bound_TS;
frac_savings_bonf_TS = (total_samples_bonf_bound_TS - total_samples_exact_TS)./total_samples_bonf_bound_TS;

%% Save data to workspace
%save('..\Fraction Savings - Figures and Workspaces\PGS_runs.mat')
save([post_obj,'_runs_070921.mat'])

%% Report average fraction savings with error estimates

% For error estimates, take mean over each row, then std over the iid means
avg_frac_savings_slep_EA = mean(mean(frac_savings_slep_EA));
se_frac_savings_slep_EA = 1.96*std(mean(frac_savings_slep_EA, 2))/sqrt(M); 
avg_frac_savings_bonf_EA = mean(mean(frac_savings_bonf_EA));
se_frac_savings_bonf_EA = 1.96*std(mean(frac_savings_bonf_EA, 2))/sqrt(M); 
avg_frac_savings_slep_OCBA = mean(mean(frac_savings_slep_OCBA));
se_frac_savings_slep_OCBA = 1.96*std(mean(frac_savings_slep_OCBA, 2))/sqrt(M);
avg_frac_savings_bonf_OCBA = mean(mean(frac_savings_bonf_OCBA));
se_frac_savings_bonf_OCBA = 1.96*std(mean(frac_savings_bonf_OCBA, 2))/sqrt(M); 
avg_frac_savings_slep_TS = mean(mean(frac_savings_slep_TS));
se_frac_savings_slep_TS = 1.96*std(mean(frac_savings_slep_TS, 2))/sqrt(M);
avg_frac_savings_bonf_TS = mean(mean(frac_savings_bonf_TS));
se_frac_savings_bonf_TS = 1.96*std(mean(frac_savings_bonf_TS, 2))/sqrt(M);

% Print table to screen
fprintf(['Bound \t p',post_obj,'_bonf \t\t p',post_obj,'_slep\n'])
fprintf('EA \t\t %.3f+/-%.3f \t %.2f+/-%.3f\n', avg_frac_savings_bonf_EA, se_frac_savings_bonf_EA, avg_frac_savings_slep_EA, se_frac_savings_slep_EA)
fprintf('OCBA \t %.3f+/-%.3f \t %.2f+/-%.3f\n', avg_frac_savings_bonf_OCBA, se_frac_savings_bonf_OCBA, avg_frac_savings_slep_OCBA, se_frac_savings_slep_OCBA)
fprintf('TS \t\t %.3f+/-%.3f \t %.2f+/-%.3f\n', avg_frac_savings_bonf_TS, se_frac_savings_bonf_TS, avg_frac_savings_slep_TS, se_frac_savings_slep_TS)

%% Plot empirical cdfs

figure
hold on
stairs(sort(reshape(frac_savings_bonf_EA,1,n_obs)), (1:n_obs)/n_obs, 'm:', 'LineWidth', 1.5);
stairs(sort(reshape(frac_savings_slep_EA,1,n_obs)), (1:n_obs)/n_obs, 'm-', 'LineWidth', 1.5);
stairs(sort(reshape(frac_savings_bonf_OCBA,1,n_obs)), (1:n_obs)/n_obs, 'b:', 'LineWidth', 1.5);
stairs(sort(reshape(frac_savings_slep_OCBA,1,n_obs)), (1:n_obs)/n_obs, 'b-', 'LineWidth', 1.5);
stairs(sort(reshape(frac_savings_bonf_TS,1,n_obs)), (1:n_obs)/n_obs, 'r:', 'LineWidth', 1.5);
stairs(sort(reshape(frac_savings_slep_TS,1,n_obs)), (1:n_obs)/n_obs, 'r-', 'LineWidth', 1.5);

xlim([0,1])
ylim([0,1])

hold off
set(gca, 'FontSize', 14, 'LineWidth', 2)

legend(['EA - ',post_obj,' Bonf'], ['EA - ',post_obj,' Slep'], ['OCBA - ',post_obj,' Bonf'], ['OCBA - ',post_obj,' Slep'], ['TS - ',post_obj,' Bonf'], ['TS - ',post_obj,' Slep'], 'Location', 'southeast')
legend boxoff

xlabel('Fractional Savings ($s$)', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold')
ylabel('$P(S \leq s$)', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold')

%rgb2gray(RGB) % Make image grayscale

%% Report errors for empirical cdfs
frac_savings_struct = {frac_savings_bonf_EA, frac_savings_slep_EA, frac_savings_bonf_OCBA, frac_savings_slep_OCBA, frac_savings_bonf_TS,frac_savings_slep_TS};

disc = 0:0.01:1; % discretize fraction savings (s)
stderrors = zeros(size(frac_savings_struct,2), length(disc));
for struct_index = 1:size(frac_savings_struct,2)
    frac_savings = frac_savings_struct{struct_index};
    for col_index = 1:length(disc)
        s = disc(col_index);
        proportions = mean((frac_savings < s), 2);
        stderrors(struct_index, col_index) = 1.96*std(proportions)/sqrt(M);
    end
    fprintf('Max std error = %f.\n', max(stderrors(struct_index,:)))
end