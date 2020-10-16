function [total_samples_bound, total_samples_exact, total_samples_MC, check_times_exact, check_times_bound, check_times_MC] = MonteCarloEfficiency(post_obj, alloc_rule, known_var)
% Run different allocation procedures (EA, TS, OCBA) for different stopping
% rule (PGS and EOC). Compare runtimes of procedures and time spent checking
% stopping rules for (1) Bound, (2) Exact, (3) MC + Exact.

% INPUTS
% -------------------------------------------------------------------------
% post_obj      : string specifying the posterior quantity of interest
%                   - 'PGS'         : Posterior PGS
%                   - 'EOC'         : Posterior EOC
% alloc_rule    : string specifying the allocation rule
%                   - 'EA'          : Equal allocation
%                   - 'OCBA-PGS'    : OCBA with PGS Slepian
%                   - 'OCBA-EOC'    : OCBA with EOC Bonferroni
%                   - 'TS'          : Thompson sampling
% known_var     : string specifying whether the true variances are known
%                   ('known' if known and 'unknown' if unknown)

% MAIN
% -------------------------------------------------------------------------

% Set parameters
k = 50; % number of systems
n0 = 10; % the initial sample size

delta = 1; % the good selection parameter
% will be unused if EOC is selected

if strcmp(post_obj, 'PGS') == 1
    threshold = 0.9; % 1 - alpha = confidence
    minmax = 1; % do NOT flip the direction of inequalities/maxima
elseif strcmp(post_obj, 'EOC') == 1    
    threshold = 0.5; % beta = EOC tolerance
    minmax = -1; % flip the direction of inequalities/maxima
end

% Run the procedure many times
M = 100; % number of macroreplications

% Monte Carlo parameters
R = 10000;

% Track the overall number of samples taken in each macroreplication
total_samples_bound = zeros(1,M); % ... using bound
total_samples_exact = zeros(1,M); % ... using numerical integration
total_samples_MC = zeros(1,M); % ... using Monte Carlo

% Warm up the tic toc function
tic;
toc;

% Track the overall number of samples taken in each macroreplication
check_times_bound = zeros(1,M); % ... using bound
check_times_exact = zeros(1,M); % ... using numerical integration
check_times_MC = zeros(1,M); % ... using Monte Carlo

% counter = 0;
% iter = 0;
% eoc = 0;
% eoc_bound = 0;

for m = 1:M % loop over macroreplications
    
    fprintf('Running macroreplication %d ...\n',m);

    stopping_flags = 0; % flag indicating how many methods have terminated
    bound_flag = 0;
    exact_flag = 0;
    MC_flag = 0;
    
    n = zeros(1, k); % initial sample size
    
    % Generate random problem instance
    % mu_i ~ -logN(0, (1.5*delta)^2) and sigma_i^2 ~ ChiSquared(4)
    mu = -lognrnd(0, 1.5*delta, 1, k);
    sigma = sqrt(chi2rnd(4, 1, k));
    
    % True configuration is slippage configuration with difference Delta
%     Delta = 1; % difference in performance for the slippage configuration
%     mu = zeros(1,k);
%     mu(1) = Delta; % System 1 is the best
%     sigma = sqrt(chi2rnd(4, 1, k));

    while stopping_flags < 3
        
        if sum(n) == 0

            % Generate initial samples from all systems
            data = mvnrnd(mu, diag(sigma.^2), n0);
            sample_means = mean(data, 1);
            sample_vars = var(data, 0, 1);
            n = n0*ones(1,k); % cumulative sample sizes (to be updated)

        else
            
            % Allocation Rule: Choose a set of systems "allocate" to sample next
            if strcmp(known_var, 'known') == 1
                next_alt = Allocate(alloc_rule, k_brack, sample_means, (sigma.^2)./n, n-1, known_var, delta);
            elseif strcmp(known_var, 'unknown') == 1
                next_alt = Allocate(alloc_rule, k_brack, sample_means, sample_vars./n, n-1, known_var, delta);
            end

            for ell_index = 1:length(next_alt)
                ell = next_alt(ell_index);
                % Take a new sample from System ell
                new_observation = normrnd(mu(ell), sigma(ell));

                % Update the sample variance, sample means, and sample size of System ell
                sample_vars(ell) = ((n(ell)-1)/n(ell))*sample_vars(ell) + (new_observation - sample_means(ell))^2/(n(ell)+1); 
                sample_means(ell) = n(ell)/(n(ell)+1)*sample_means(ell) + (1/(n(ell)+1))*new_observation; % shortcut update formula
                n(ell) = n(ell) + 1;
            end

        end
        
        % Check stopping conditions
        
        % Find the index of the best-looking alternative
        tic;
        [~, k_brack] = max(sample_means);
        increment = toc;
        % Add time to all methods still running
        if bound_flag == 0
            check_times_bound(m) = check_times_bound(m) + increment;
        end
        if exact_flag == 0
            check_times_exact(m) = check_times_exact(m) + increment;
        end
        if MC_flag == 0
            check_times_MC(m) = check_times_MC(m) + increment; 
        end
        
        % Calculate Bonf bound on posterior EOC for best-looking
        % alternative
        if bound_flag == 0
            tic;
            if strcmp(known_var, 'known') == 1
                post_qty_bound = CalcPostQty('EOC-Bonf', k_brack, sample_means, (sigma.^2)./n, n-1, known_var, delta);
            elseif strcmp(known_var, 'unknown') == 1
                post_qty_bound = CalcPostQty('EOC-Bonf', k_brack, sample_means, sample_vars./n, n-1, known_var, delta);
            end
            check_times_bound(m) = check_times_bound(m) + toc;

            % Check if bound method stops
            if minmax*post_qty_bound > minmax*threshold
                stopping_flags = stopping_flags + 1;
                bound_flag = 1;
                total_samples_bound(m) = sum(n);
            end
        end
        
        % Calculate exact posterior EOC for best-looking alternative
        if exact_flag == 0
            tic;
            if strcmp(known_var, 'known') == 1
                best_post_qty = CalcPostQty('EOC-Exact1', k_brack, sample_means, (sigma.^2)./n, n-1, known_var, delta);
            elseif strcmp(known_var, 'unknown') == 1
                best_post_qty = CalcPostQty('EOC-Exact1', k_brack, sample_means, sample_vars./n, n-1, known_var, delta);
            end
            check_times_exact(m) = check_times_exact(m) + toc;

            % Check if exact method stops
            if minmax*best_post_qty > minmax*threshold
                stopping_flags = stopping_flags + 1;
                exact_flag = 1;
                total_samples_exact(m) = sum(n);
            end
        end
        
        % Calculate MC posterior EOC for best-looking alternative
        if MC_flag == 0
            tic;
            %W_draw_matrix = zeros(R, k);
            if strcmp(known_var, 'known') == 1 % known variances --> normal distribution
                W_draw_matrix = mvnrnd(sample_means, diag(sigma.^2./n), R);
%                 for r = 1:R
%                     W_draw_matrix(r,:) = normrnd(sample_means, sqrt((sigma.^2)./n));
%                 end
            elseif strcmp(known_var, 'unknown') == 1 % unknown variances --> t distribution
                %W_draw_matrix = sample_means + sqrt(sample_vars./n).*mvtrnd(eye(k), n-1, R);
                %for r = 1:R
                %    W_draw_matrix(r,:) = sample_means + sqrt(sample_vars./n).*trnd(n-1);
                %end
                for i = 1:k
                    W_draw_matrix(:, i) = sample_means(i) + sqrt(sample_vars(i)/n(i)).*trnd(n(i)-1, [R, 1]);
                end
            end
            pEOC_MC = mean(max(W_draw_matrix, [], 2) - W_draw_matrix(:,k_brack));
            check_times_MC(m) = check_times_MC(m) + toc;
            
            if minmax*pEOC_MC > minmax*threshold % --> Compute exact quantity
                
                tic;
                if strcmp(known_var, 'known') == 1
                    best_post_qty = CalcPostQty('EOC-Exact1', k_brack, sample_means, (sigma.^2)./n, n-1, known_var, delta);
                elseif strcmp(known_var, 'unknown') == 1
                    best_post_qty = CalcPostQty('EOC-Exact1', k_brack, sample_means, sample_vars./n, n-1, known_var, delta);
                end
                check_times_MC(m) = check_times_MC(m) + toc;

                if minmax*best_post_qty > minmax*threshold
                    stopping_flags = stopping_flags + 1;
                    MC_flag = 1;
                    total_samples_MC(m) = sum(n);
                end
            end
        end
    end
end