clear
clc

% Report the time needed to compute the exact posterior PGS and exact 
% posterior EOC as a function of the number of alternatives, k.

% Set parameters
k_vector = [10, 100, 1000, 10000]; %, 100000]; % number of systems
%k_vector = [10, 100, 1000]; % number of systems
delta = 1; % the good selection parameter... shouldn't affect the runtimes
num_obs = 100; % total number of macroreplications for each k
%num_obs = 1000; % total number of macroreplications for each k


% Initialize matrix for output
times = zeros(length(k_vector), num_obs);

% Warm up the tic toc function
tic;
toc;

for k_index = 1:length(k_vector)
    k = k_vector(k_index);

    fprintf('Running replications for k = %d.\n',k);
    
    for m = 1:num_obs
        if mod(m,100) == 0
            fprintf('\tGenerating RPI %d of %d\n', m, num_obs);
        end
        
        % Generate a random posterior instance (RPI) where
        % mu_i ~ N(0, 25*delta^2) and sigma_i^2 ~ ChiSquared(4)
        post_means = normrnd(0, 5*delta, 1, k);
        post_vars = chi2rnd(4, 1, k);
        % Assuming known variances... shouldn't affect the runtimes
        
        % Calculate the exact posterior quantity of the Alternative 1
        % Record the computational time
        tic;
        post_qty = CalcPostQty('PGS-Exact', 1, post_means, post_vars, 0, 'known', delta); 
        %post_qty = CalcPostQtyPlusError('PGS-Exact', 1, post_means, post_vars, 0, 'known', delta); 
        times(k_index, m) = toc;

    end
   
end

save('times_errors_ppgs_exact.mat')
%save('times_errors_peoc_exact_1int.mat')
%save('times_errors_peoc_exact_2int.mat')

% Record summary statistics
mean_time = mean(times,2);
z_alpha = 1.96;
halfwidth = z_alpha*std(times,0,2)/sqrt(num_obs);
lb = mean_time - z_alpha*std(times,0,2)/sqrt(num_obs);
ub = mean_time + z_alpha*std(times,0,2)/sqrt(num_obs);

fprintf('Number of Alternatives \t %d \t %d \t %d \t %d \t %d \n', k_vector(1), k_vector(2), k_vector(3), k_vector(4), k_vector(5))
fprintf('Average Time (seconds) \t %.3f +/- %.3f \t %.3f +/- %.3f \t %.3f +/- %.3f \t %.3f +/- %.3f \t %.3f +/- %.3f \n', mean_time(1), halfwidth(1), mean_time(2), halfwidth(2), mean_time(3), halfwidth(3), mean_time(4), halfwidth(4), mean_time(5), halfwidth(5))