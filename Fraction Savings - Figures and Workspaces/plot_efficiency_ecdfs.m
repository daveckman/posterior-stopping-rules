% Plot ECDF of fractional savings

post_quantity = 'PGS'; %'EOC'
rpi_mode = 3;
M = 100;
Q = 50;

% Load data
load([post_quantity,'_rpi_mode=',num2str(rpi_mode),'_M=',num2str(M),'_Q=',num2str(Q),'.mat'])

%%

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
fprintf('EA \t\t %.2f+/-%.2f \t %.2f+/-%.2f\n', avg_frac_savings_bonf_EA, se_frac_savings_bonf_EA, avg_frac_savings_slep_EA, se_frac_savings_slep_EA)
fprintf('OCBA \t %.2f+/-%.2f \t %.2f+/-%.2f\n', avg_frac_savings_bonf_OCBA, se_frac_savings_bonf_OCBA, avg_frac_savings_slep_OCBA, se_frac_savings_slep_OCBA)
fprintf('TS \t\t %.2f+/-%.2f \t %.2f+/-%.2f\n', avg_frac_savings_bonf_TS, se_frac_savings_bonf_TS, avg_frac_savings_slep_TS, se_frac_savings_slep_TS)

%% Plot empirical cdfs

figure
hold on
stairs(sort(reshape(frac_savings_bonf_EA,1,n_obs)), (1:n_obs)/n_obs, 'color', '#1b0ef8', 'linestyle', '--', 'LineWidth', 2);
stairs(sort(reshape(frac_savings_slep_EA,1,n_obs)), (1:n_obs)/n_obs, 'color', '#1b0ef8', 'linestyle', '-', 'LineWidth', 2);
stairs(sort(reshape(frac_savings_bonf_OCBA,1,n_obs)), (1:n_obs)/n_obs, 'color', '#ae1fe5', 'linestyle', '--', 'LineWidth', 2);
stairs(sort(reshape(frac_savings_slep_OCBA,1,n_obs)), (1:n_obs)/n_obs, 'color', '#ae1fe5', 'linestyle', '-', 'LineWidth', 2);
stairs(sort(reshape(frac_savings_bonf_TS,1,n_obs)), (1:n_obs)/n_obs, 'color', '#f67c80', 'linestyle', '--', 'LineWidth', 2);
stairs(sort(reshape(frac_savings_slep_TS,1,n_obs)), (1:n_obs)/n_obs, 'color', '#f67c80', 'linestyle', '-', 'LineWidth', 2);

xlim([0,1])
ylim([0,1])

hold off
set(gca, 'FontSize', 14, 'LineWidth', 2)

legend(['EA: ',post_obj,' Bonf'], ['EA: ',post_obj,' Slep'], ['OCBA: ',post_obj,' Bonf'], ['OCBA: ',post_obj,' Slep'], ['TS: ',post_obj,' Bonf'], ['TS: ',post_obj,' Slep'], 'Location', 'southeast')
legend boxoff

xlabel('Fractional Savings ($s$)', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold')
ylabel('$P(S \leq s)$', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold')

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
    fprintf('Max half width (1.96*std error) = %f.\n', max(stderrors(struct_index,:)))
end