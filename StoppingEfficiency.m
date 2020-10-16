function [total_samples_bound, total_samples_exact] = StoppingEfficiency(post_obj, alloc_rule, known_var)
% Run different allocation procedures (EA, TS, OCBA) for different stopping
% rule (PGS and EOC). Compare runtimes of procedures that use bounds to
% check the stopping rule vs procedures that use the exact quantity.

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

delta = 1; % the good selection parameter (unused if EOC is selected)

if strcmp(post_obj, 'PGS') == 1
    threshold = 0.9; % 1 - alpha = confidence
    minmax = 1; % do NOT flip the direction of inequalities/maxima
elseif strcmp(post_obj, 'EOC') == 1    
    threshold = 0.5; % beta = EOC tolerance
    minmax = -1; % flip the direction of inequalities/maxima
end

% Run the procedure many times
M = 100; % number of macroreplications
S = 50; % number of splits per macroreplication (for variance reduction)
% S = 1 corresponds to standard Monte Carlo

% Track the overall number of samples taken in each macroreplication
total_samples_bound = zeros(1,M*S); % ... using bound
total_samples_exact = zeros(1,M*S); % ... using numerical integration

% counter = 0;
% iter = 0;
% eoc = 0;
% eoc_bound = 0;

for m = 1:M % loop over macroreplications
    
    fprintf('Running macroreplication %d ...\n',m);

    % Generate random problem instance
    % mu_i ~ -logN(0, (1.5*delta)^2) and sigma_i^2 ~ ChiSquared(4)
    mu = -lognrnd(0, 1.5*delta, 1, k);
    sigma = sqrt(chi2rnd(4, 1, k));
    
    % True configuration is slippage configuration with difference Delta
%     Delta = 1; % difference in performance for the slippage configuration
%     mu = zeros(1,k);
%     mu(1) = Delta; % System 1 is the best
%     sigma = sqrt(chi2rnd(4, 1, k));
    
    % Generate initial samples from all systems
    data = mvnrnd(mu, diag(sigma.^2), n0);
    sample_means = mean(data, 1);
    sample_vars = var(data, 0, 1);
    n = n0*ones(1,k); % cumulative sample sizes (to be updated)

    % Find the index of the best-looking alternative
    [~, k_brack] = max(sample_means);

    % Find the best posterior quantity
    if strcmp(post_obj, 'PGS') == 1

        % Calculate exact posterior PGS for candidate alternatives
        if strcmp(known_var, 'known') == 1
            CandidateIndex = FindMaxPGSCandidatesSort(delta, k_brack, sample_means, (sigma.^2)./n);
        elseif strcmp(known_var, 'unknown') == 1
            CandidateIndex = FindMaxPGSCandidatesSort(delta, k_brack, sample_means, sample_vars./n);
        end

        PGS_exact_vector = zeros(1,length(CandidateIndex));
        for i = 1:length(CandidateIndex)
            if strcmp(known_var, 'known') == 1
                PGS_exact_vector(i) = CalcPostQty('PGS-Exact', CandidateIndex(i), sample_means, (sigma.^2)./n, n-1, known_var, delta);
            elseif strcmp(known_var, 'unknown') == 1
                PGS_exact_vector(i) = CalcPostQty('PGS-Exact', CandidateIndex(i), sample_means, sample_vars./n, n-1, known_var, delta);
            end
        end
        best_post_qty = max(PGS_exact_vector); % Find the highest posterior PGS

    elseif strcmp(post_obj, 'EOC') == 1
        
        % Calculate exact posterior EOC for best-looking alternative
        if strcmp(known_var, 'known') == 1
            best_post_qty = CalcPostQty('EOC-Exact1', k_brack, sample_means, (sigma.^2)./n, n-1, known_var, delta);
        elseif strcmp(known_var, 'unknown') == 1
            best_post_qty = CalcPostQty('EOC-Exact1', k_brack, sample_means, sample_vars./n, n-1, known_var, delta);
        end
    end
    
    
    while minmax*best_post_qty < minmax*threshold % ... until the exact PGS stopping condition is met

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

        % Find the index of the best-looking alternative
        [~, k_brack] = max(sample_means);

         % Find the best posterior quantity
        if strcmp(post_obj, 'PGS') == 1

            % Calculate exact posterior PGS for candidate alternatives
            if strcmp(known_var, 'known') == 1
                CandidateIndex = FindMaxPGSCandidatesSort(delta, k_brack, sample_means, (sigma.^2)./n);
            elseif strcmp(known_var, 'unknown') == 1
                CandidateIndex = FindMaxPGSCandidatesSort(delta, k_brack, sample_means, sample_vars./n);
            end

            PGS_exact_vector = zeros(1,length(CandidateIndex));
            for i = 1:length(CandidateIndex)
                if strcmp(known_var, 'known') == 1
                    PGS_exact_vector(i) = CalcPostQty('PGS-Exact', CandidateIndex(i), sample_means, (sigma.^2)./n, n-1, known_var, delta);
                elseif strcmp(known_var, 'unknown') == 1
                    PGS_exact_vector(i) = CalcPostQty('PGS-Exact', CandidateIndex(i), sample_means, sample_vars./n, n-1, known_var, delta);
                end
            end
            best_post_qty = max(PGS_exact_vector); % Find the highest posterior PGS

        elseif strcmp(post_obj, 'EOC') == 1

            % Calculate exact posterior EOC for best-looking alternative
            if strcmp(known_var, 'known') == 1
                best_post_qty = CalcPostQty('EOC-Exact1', k_brack, sample_means, (sigma.^2)./n, n-1, known_var, delta);
            elseif strcmp(known_var, 'unknown') == 1
                best_post_qty = CalcPostQty('EOC-Exact1', k_brack, sample_means, sample_vars./n, n-1, known_var, delta);
            end
        end

%         counter = counter + 1;
%         iter = [iter, counter];
%         eoc = [eoc, best_post_qty];
%         eoc_bound = [eoc_bound, 0];
    end

    % Fix stats at the time of the exact PGS stopping condition
    sample_means_stop = sample_means;
    sample_vars_stop = sample_vars;
    n_stop = n;
    n_exact = n;

    % SPLITTING

    for s = 1:S % loop over the splits

        if mod(s, 10) == 0
            fprintf('  Split %d ...\n',s);
        end

        % Reset the stats from the stopping time
        sample_means = sample_means_stop;
        sample_vars = sample_vars_stop;
        n = n_stop;

        % Identify the best looking system and calculate bound on its
        % posterior quantity
        [~, k_brack] = max(sample_means);
        
        if strcmp(post_obj, 'PGS') == 1
            if strcmp(known_var, 'known') == 1
                post_qty_bound = CalcPostQty('PGS-Slep', k_brack, sample_means, (sigma.^2)./n, n-1, known_var, delta);
            elseif strcmp(known_var, 'unknown') == 1
                post_qty_bound = CalcPostQty('PGS-Slep', k_brack, sample_means, sample_vars./n, n-1, known_var, delta);
            end
        elseif strcmp(post_obj, 'EOC') == 1
            if strcmp(known_var, 'known') == 1
                post_qty_bound = CalcPostQty('EOC-Bonf', k_brack, sample_means, (sigma.^2)./n, n-1, known_var, delta);
            elseif strcmp(known_var, 'unknown') == 1
                post_qty_bound = CalcPostQty('EOC-Bonf', k_brack, sample_means, sample_vars./n, n-1, known_var, delta);
            end
        end
        
        while minmax*post_qty_bound < minmax*threshold % ... until the bound stopping condition is met

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

                % Update the sample variance, sample mean, and sample size of System ell
                sample_vars(ell) = ((n(ell)-1)/n(ell))*sample_vars(ell) + (new_observation - sample_means(ell))^2/(n(ell)+1); 
                sample_means(ell) = n(ell)/(n(ell)+1)*sample_means(ell) + (1/(n(ell)+1))*new_observation; % shortcut update formula
                n(ell) = n(ell) + 1;
            end

            % Identify the best looking system and calculate bound on its
            % posterior quantity
            [~, k_brack] = max(sample_means);

            if strcmp(post_obj, 'PGS') == 1
                if strcmp(known_var, 'known') == 1
                    post_qty_bound = CalcPostQty('PGS-Slep', k_brack, sample_means, (sigma.^2)./n, n-1, known_var, delta);
                elseif strcmp(known_var, 'unknown') == 1
                    post_qty_bound = CalcPostQty('PGS-Slep', k_brack, sample_means, sample_vars./n, n-1, known_var, delta);
                end
            elseif strcmp(post_obj, 'EOC') == 1
                if strcmp(known_var, 'known') == 1
                    post_qty_bound = CalcPostQty('EOC-Bonf', k_brack, sample_means, (sigma.^2)./n, n-1, known_var, delta);
                elseif strcmp(known_var, 'unknown') == 1
                    post_qty_bound = CalcPostQty('EOC-Bonf', k_brack, sample_means, sample_vars./n, n-1, known_var, delta);
                end
            end

%             counter = counter + 1;
%             iter = [iter, counter];
%             eoc = [eoc, best_post_qty];
%             eoc_bound = [eoc_bound, post_qty_bound];

            
        end

        rep_index = S*(m-1) + s;
        total_samples_exact(rep_index) = sum(n_exact);
        total_samples_bound(rep_index) = sum(n);

        %fprintf('\t (Exact): \t Total # of Samples = %d\n', total_samples_exact(rep_index));
        %fprintf('\t (Bound): \t Total # of Samples = %d\n', total_samples_bound(rep_index)); 
    end
    
end

% hold on
% plot(iter, eoc, 'b-')
% plot(iter, eoc_bound, 'r-')
% hold off