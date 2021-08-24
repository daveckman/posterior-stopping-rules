function [next_alt] = Allocate(alloc_rule, k_brack, post_means, post_vars, post_dfs, known_var, delta)
% Determine the next alternative(s) to simulate given an allocation rule
% Assumes outputs of simulation replications are normally distributed

% INPUTS
% -------------------------------------------------------------------------
% alloc_rule    : string specifying the allocation rule
%                   - 'EA'          : Equal allocation
%                   - 'OCBA-PGS'    : OCBA with PGS Slepian
%                   - 'OCBA-EOC'    : OCBA with EOC Bonferroni
%                   - 'TS'          : Thompson sampling
% k_brack       : index of alternative with the highest posterior mean
% post_means    : a (1 x k) vector of posterior means for W_1, ..., W_k
% post_vars     : a (1 x k) vector of posterior variances for W_1, ..., W_k
% post_dfs      : a (1 x k) vector of posterior degrees of freedom for W_1, ..., W_k
%               (post_dfs will not be used unless variances are unknown)
% known_var     : string specifying whether the true variances are known
%               ('known' if known and 'unknown' if unknown)

% OUTPUTS
% -------------------------------------------------------------------------
% next_alt     : index(es) of next alternative(s) to simulate

% MAIN
% -------------------------------------------------------------------------

k = length(post_means);
%[~, k_brack] = max(post_means);

if strcmp(alloc_rule, 'EA') == 1
 
    next_alt = 1:k;
    
elseif strcmp(alloc_rule, 'OCBA-PGS') == 1
    
    % OCBA allocation rule with PGS
    EAPGS_slep_vector = zeros(1,k);
    
    for i = 1:k
        
        % Pretend that an extra sample is taken from Alternative i
        % Update the posterior variance, but leave the posterior mean
        % unchanged % df = n_i - 1
        new_post_vars = post_vars;
        new_post_vars(i) = post_vars(i)*((post_dfs(i)+1)/(post_dfs(i)+2));
        
        EAPGS_slep_vector(i) = CalcPostQty('PGS-Slep', k_brack, post_means, new_post_vars, post_dfs, known_var, delta);

    end
    
    % Pick the alternative that leads to the largest expected approximate 
    % PGS of Alternative k_brack
    [~, next_alt] = max(EAPGS_slep_vector);

elseif strcmp(alloc_rule, 'OCBA-EOC') == 1
    
    % OCBA allocation rule with EOC
    EAEOC_bonf_vector = zeros(1,k);
    
    for i = 1:k
        
        % Pretend that an extra sample is taken from Alternative i
        % Update the posterior variance, but leave the posterior mean
        % unchanged % df = n_i - 1
        new_post_vars = post_vars;
        new_post_vars(i) = post_vars(i)*((post_dfs(i)+1)/(post_dfs(i)+2));
        
        % Calculate the expected approximate EOC of Alternative k_brack
        EAEOC_bonf_vector(i) = CalcPostQty('EOC-Bonf', k_brack, post_means, new_post_vars, post_dfs, known_var, 0);
    end

    % Pick the alternative that leads to the smallest expected approximate 
    % EOC of Alternative k_brack
    [~, next_alt] = min(EAEOC_bonf_vector);
    
elseif strcmp(alloc_rule, 'TS') == 1
    
    % Thompson sampling
    
    if strcmp(known_var, 'known') == 1 % known variances --> normal distribution

        thompson_draw = normrnd(post_means, post_vars);

    elseif strcmp(known_var, 'unknown') == 1 % unknown variances --> t distribution

        thompson_draw = post_means + sqrt(post_vars).*trnd(post_dfs);

    end
    
    % Pick the best-looking alternative from the sample
    [~, next_alt] = max(thompson_draw);

end % end cases

end % end function




