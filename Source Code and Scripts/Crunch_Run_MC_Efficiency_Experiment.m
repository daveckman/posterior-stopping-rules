function Crunch_Run_MC_Efficiency_Experiment(M, pqi_mode, rpi_mode)

fprintf('M = %d and pqi_mode = %d and rpi_mode = %d.\n', M, pqi_mode, rpi_mode)

cluster = parcluster;
%pool_obj = parpool(cluster);
maxNumCompThreads(8);

% Experimental setup

%post_obj = 'PGS'; % string ('PGS' or 'EOC')
known_var = 'unknown'; % unknown variances
alloc_rule = 'TS';
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

%% Run EOC stopping and OCBA allocation
[total_samples_Bonf_bound, total_samples_Slep_bound, total_samples_exact, total_samples_MC, total_samples_pure_MC, check_times_Bonf_bound, check_times_Slep_bound, check_times_exact, check_times_MC, check_times_pure_MC, pure_MC_stop_post_qty] = CrunchMonteCarloEfficiency(cluster, rpi_mode, post_obj, alloc_rule, known_var, k, n0, M);

%% Save data to workspace
save(['data/MonteCarlo_',alloc_rule,'_',post_obj,'_Exact1_rpi_mode=',num2str(rpi_mode),'_M=',num2str(M),'.mat'])

