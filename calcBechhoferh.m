function h = calcBechhoferh(k, alpha)

% Use simulation to estimate Bechhofer h
% h is the 1-alpha quantile of the maximum of a (k-1) dimensional MVN with mean 0,
% variance 1, and positive correlation 1/2.

R = 1000000; % Number of replications
Mu = zeros(1,k-1); % Mean vector
Sigma = 0.5*ones(k-1,k-1) + 0.5*diag(ones(k-1,1)); % Covariance matrix
samples = mvnrnd(Mu, Sigma, R); % Generate random vectors
max_samples = max(samples, [], 2); % Identify max
sort_max = sort(max_samples); % Sort max
h = sort_max(floor(R*(1-alpha))); % Find (1-alpha)-quantile.