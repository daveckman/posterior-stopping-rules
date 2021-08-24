function Crunch_Run_Timing(M, pqi_mode, var_mode)

fprintf('num_obs = %d and pqi_mode = %d and var_mode = %d.\n', M, pqi_mode, var_mode)

% M = total number of macroreplications for each k

%cluster = parcluster;
%pool_obj = parpool(cluster);
%maxNumCompThreads(8);

% Report the time needed to compute the exact posterior PGS and exact 
% posterior EOC as a function of the number of alternatives, k.

% Set parameters
delta = 1; % the good selection parameter... shouldn't affect the runtimes

if pqi_mode == 1 % pPGS
    pqi_string = 'PGS-Exact';
    slep_string = 'PGS-Slep';
    bonf_string = 'PGS-Bonf';
    k_vector = [10, 100, 1000, 10000, 100000];
elseif pqi_mode == 2 % pEOC single integral
    pqi_string = 'EOC-Exact1';
    slep_string = 'EOC-Slep';
    bonf_string = 'EOC-Bonf';
    k_vector = [10, 100, 1000];
elseif pqi_mode == 3 % pEOC double integral
    pqi_string = 'EOC-Exact2';
    slep_string = 'EOC-Slep';
    bonf_string = 'EOC-Bonf';
    k_vector = [10, 100, 1000];
end

if var_mode == 1 % known variances
    var_string = 'known';
    df = 0;
elseif var_mode == 2 % unknown variances
    var_string = 'unknown';
    df = 9;
end

% Initialize matrix for output
times = zeros(length(k_vector), M);
pqtys = zeros(length(k_vector), M);
pqtys_slep = zeros(length(k_vector), M);
pqtys_bonf = zeros(length(k_vector), M);

% Reset random number stream
rng default

% Warm up the tic toc function
tic;
toc;

for k_index = 1:length(k_vector)
    k = k_vector(k_index);

    % Generate all random posterior instances (RPI) up front
    % mu_i ~ N(0, 25*delta^2) and sigma_i^2 ~ ChiSquared(4)
    %post_means_matrix = normrnd(0, 5*delta, [M, k]);
    %post_vars_matrix = chi2rnd(4, [M, k]);
    
    fprintf('Running replications for k = %d.\n',k);
    
    %parfor (m = 1:M, cluster)
    for m = 1:M
        if mod(m,100) == 0
            fprintf('\tTesting RPI %d of %d\n', m, M);
        end
        
        % Extract random problem instance
        %post_means = post_means_matrix(m,:);
        %post_vars = post_vars_matrix(m,:);
        post_means = normrnd(0, 5*delta, [1, k]);
        post_vars = chi2rnd(4, [1, k]);
        
        % Calculate the exact posterior quantity of the Alternative 1
        [~, best] = max(post_means);
        
        % Assuming known variances (no degrees of freedom) ... shouldn't affect the runtimes
        % Record the computational time
        
        if m == 1 % warmup integral function
            CalcPostQty(pqi_string, best, post_means, post_vars, df*ones(1,k), var_string, delta);
        end
        
        tic;
        pqtys(k_index, m) = CalcPostQty(pqi_string, best, post_means, post_vars, df*ones(1,k), var_string, delta);
        times(k_index, m) = toc;
        pqtys_slep(k_index, m) = CalcPostQty(slep_string, best, post_means, post_vars, df*ones(1,k), var_string, delta);
        pqtys_bonf(k_index, m) = CalcPostQty(bonf_string, best, post_means, post_vars, df*ones(1,k), var_string, delta);
    end
   
end

% Record summary statistics
mean_time = mean(times,2);
z_alpha = 1.96;
halfwidth = z_alpha*std(times,0,2)/sqrt(M);
lb = mean_time - z_alpha*std(times,0,2)/sqrt(M);
ub = mean_time + z_alpha*std(times,0,2)/sqrt(M);
q_10pct = quantile(times, 0.1, 2);
q_90pct = quantile(times, 0.9, 2);

for k_index = 1:length(k_vector)
    fprintf('Number of Alternatives \t %d \n', k_vector(k_index))
    fprintf('Average Time (seconds) \t %.3f +/- %.3f \n', mean_time(k_index), halfwidth(k_index))
    fprintf('0.10 Quantile \t %.3f \n', q_10pct(k_index))
    fprintf('0.90 Quantile \t %.3f \n', q_90pct(k_index)) 
end

save(['times_errors_',pqi_string,'_',var_string,'_M=',num2str(M),'.mat'])