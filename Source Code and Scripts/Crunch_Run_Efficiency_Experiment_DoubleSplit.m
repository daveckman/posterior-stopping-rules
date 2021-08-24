function Crunch_Run_Efficiency_Experiment_DoubleSplit(M, Q, pqi_mode, rpi_mode)

fprintf('M = %d and Q = %d and pqi_mode = %d and rpi_mode = %d.\n', M, Q, pqi_mode, rpi_mode)

cluster = parcluster;
%pool_obj = parpool(cluster);
maxNumCompThreads(8);

% Experimental setup

%post_obj = 'PGS'; % string ('PGS' or 'EOC')
known_var = 'unknown'; % unknown variances
k = 100; % number of systems
n0 = 5; % common initial sample size
%M = 100; % number of macroreplications
%Q = 50; % number of splits
%rpi_mode = 1 or 2 % Weibull distribution

% Set posterior quantity of interest
if pqi_mode == 1
    post_obj = 'PGS';
elseif pqi_mode == 2
    post_obj = 'EOC';
end
        
% Experiments return matrices of size M rows x Q cols
n_obs = M*Q; % common number of observations

%% Run experiments

% Run stopping and EA allocation
[total_samples_slep_bound_EA, total_samples_bonf_bound_EA, total_samples_exact_EA] = CrunchStoppingEfficiencyDoubleSplit(cluster, rpi_mode, post_obj, 'EA', known_var, k, n0, M, Q);
frac_savings_slep_EA = (total_samples_slep_bound_EA - total_samples_exact_EA)./total_samples_slep_bound_EA;
frac_savings_bonf_EA = (total_samples_bonf_bound_EA - total_samples_exact_EA)./total_samples_bonf_bound_EA;

% Run stopping and OCBA allocation
[total_samples_slep_bound_OCBA, total_samples_bonf_bound_OCBA, total_samples_exact_OCBA] = CrunchStoppingEfficiencyDoubleSplit(cluster, rpi_mode, post_obj, ['OCBA-',post_obj], known_var, k, n0, M, Q);
frac_savings_slep_OCBA = (total_samples_slep_bound_OCBA - total_samples_exact_OCBA)./total_samples_slep_bound_OCBA;
frac_savings_bonf_OCBA = (total_samples_bonf_bound_OCBA - total_samples_exact_OCBA)./total_samples_bonf_bound_OCBA;

% Run stopping and EA allocation
[total_samples_slep_bound_TS, total_samples_bonf_bound_TS, total_samples_exact_TS] = CrunchStoppingEfficiencyDoubleSplit(cluster, rpi_mode, post_obj, 'TS', known_var, k, n0, M, Q);
frac_savings_slep_TS = (total_samples_slep_bound_TS - total_samples_exact_TS)./total_samples_slep_bound_TS;
frac_savings_bonf_TS = (total_samples_bonf_bound_TS - total_samples_exact_TS)./total_samples_bonf_bound_TS;

%% Save data to workspace
%save('../Fraction Savings - Figures and Workspaces\PGS_runs.mat')
save(['data/',post_obj,'_rpi_mode=',num2str(rpi_mode),'_M=',num2str(M),'_Q=',num2str(Q),'.mat'])