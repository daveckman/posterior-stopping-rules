function [post_qty] = CalcPostQty(qty_type, i, post_means, post_vars, post_dfs, known_var, delta)
% Calculate a posterior quantity of interest for Alternative i
% Assumes outputs of simulation replications are normally distributed

% INPUTS
% -------------------------------------------------------------------------
% qty_type      : string specifying the type of posterior quantity
%                   - 'PGS-Bonf'    : Bonferroni bound for posterior PGS
%                   - 'PGS-Slep'    : Slepian bound for posterior PGS
%                   - 'PGS-Exact'   : Exact posterior PGS
%                   - 'EOC-Bonf'    : Bonferroni bound for posterior EOC
%                   - 'EOC-Slep'    : Slepian bound for posterior EOC
%                   - 'EOC-Exact1'  : Exact posterior EOC ( 1 integral )
%                   - 'EOC-Exact2'  : Exact posterior EOC ( 2 integral )
% i             : index of alternative for calculating posterior quantity
% post_means    : a (1 x k) vector of posterior means for W_1, ..., W_k
% post_vars     : a (1 x k) vector of posterior variances for W_1, ..., W_k
% post_dfs      : a (1 x k) vector of posterior degrees of freedom for W_1, ..., W_k
%               (post_dfs will not be used unless variances are unknown)
% known_var     : string specifying whether the true variances are known
%               ('known' if known and 'unknown' if unknown)

% OUTPUTS
% -------------------------------------------------------------------------
% post_qty     : posterior quantity of interest of Alternative i

% MAIN
% -------------------------------------------------------------------------

k = length(post_means);

offset_sd = 5; % number of standard deviations to offset +/- in product integrals

if strcmp(qty_type, 'PGS-Bonf') == 1
    
    % Calculate means and variances of the differences W_i - W_j for all j != i
    diff_means = post_means(i) - post_means;
    diff_means = diff_means([1:i-1,i+1:k]);
    sum_var = post_vars(i) + post_vars;
    sum_var = sum_var([1:i-1,i+1:k]);
    
    if strcmp(known_var, 'known') == 1 % known variances --> normal distribution
    
        % Calculate Bonferroni bound for posterior PGS of Alternative i
        post_qty = 1 - sum(normcdf((-delta - diff_means)./sqrt(sum_var)));

    elseif strcmp(known_var, 'unknown') == 1 % unknown variances --> t distribution

        % Calculate df of the differences W_i - W_j for all j != i
        % Welch approximation
        df_denominator = post_vars(i)^2/post_dfs(i) + post_vars.^2./post_dfs;
        df_denominator = df_denominator([1:i-1,i+1:k]);
        diff_df = (sum_var.^2)./df_denominator;

        % Calculate Bonferroni bound for posterior PGS of Alternative i
        post_qty = 1 - sum(tcdf((-delta - diff_means)./sqrt(sum_var), diff_df));

    else

        disp('Argument for unknown/known variance should be "known" or "unknown".')
        return

    end
    
elseif strcmp(qty_type, 'PGS-Slep') == 1
    
    % Calculate means and variances of the differences W_i - W_j for all j != i
    diff_means = post_means(i) - post_means;
    diff_means = diff_means([1:i-1,i+1:k]);
    sum_var = post_vars(i) + post_vars;
    sum_var = sum_var([1:i-1,i+1:k]);
    
    if strcmp(known_var, 'known') == 1 % known variances --> normal distribution
    
        % Calculate Slepian bound for posterior PGS of Alternative i
        post_qty = prod(normcdf((delta + diff_means)./sqrt(sum_var)));

    elseif strcmp(known_var, 'unknown') == 1 % unknown variances --> t distribution

        % Calculate df of the differences W_i - W_j for all j != i
        % Welch approximation
        df_denominator = post_vars(i)^2/post_dfs(i) + post_vars.^2./post_dfs;
        df_denominator = df_denominator([1:i-1,i+1:k]);
        diff_df = (sum_var.^2)./df_denominator;

        % Calculate Slepian bound for posterior PGS of Alternative i
        post_qty = prod(tcdf((delta + diff_means)./sqrt(sum_var), diff_df));

    else

        disp('Argument for unknown/known variance should be "known" or "unknown".')
        return
    
    end

elseif strcmp(qty_type, 'PGS-Exact') == 1
    
    % Exclude Alternative i
    minus_means = post_means([1:i-1,i+1:k]);
    minus_var = post_vars([1:i-1,i+1:k]);
    
    if strcmp(known_var, 'known') == 1 % known variances --> normal distribution
        
        % Calculate exact posterior PGS of Alternative i
        prodfun = @(w) prod(normcdf((delta + w - minus_means)./sqrt(minus_var))).*normpdf(w, post_means(i), sqrt(post_vars(i)));
        post_qty = integral(@(w)prodfun(w), -Inf, Inf, 'ArrayValued', true, 'Waypoints', [post_means(i) - offset_sd*sqrt(post_vars(i)), post_means(i) + offset_sd*sqrt(post_vars(i))]);
        %post_qty = integral(@(w)prodfun(w), -Inf, Inf, 'ArrayValued', true); %, 'AbsTol', 1e-6);

    elseif strcmp(known_var, 'unknown') == 1 % unknown variances --> t distribution
    
        % Exclude Alternative i
        minus_df = post_dfs([1:i-1,i+1:k]);
        
        % Calculate exact posterior PGS of Alternative i
        % Integrate the product of k-1 3-parameter t cdfs and one 3-parameter t pdf
        prodfun = @(w) prod(tcdf((delta + w - minus_means)./sqrt(minus_var), minus_df)).*(1/sqrt(post_vars(i))).*tpdf((w - post_means(i))/sqrt(post_vars(i)), post_dfs(i));
        post_qty = integral(@(w)prodfun(w), -Inf, Inf, 'ArrayValued', true, 'Waypoints', [post_means(i) - offset_sd*sqrt(post_vars(i)*post_dfs(i)/(post_dfs(i) - 2)), post_means(i) + offset_sd*sqrt(post_vars(i)*post_dfs(i)/(post_dfs(i) - 2))]);
        %post_qty = integral(@(w)prodfun(w), -Inf, Inf, 'ArrayValued', true); %, 'AbsTol', 1e-6);

        
    else 
        
        disp('Argument for unknown/known variance should be "known" or "unknown".')
        return
    
    end

elseif strcmp(qty_type, 'EOC-Bonf') == 1
    
    % Calculate means and variances of the differences W_i - W_j for all j != i
    diff_means = post_means(i) - post_means;
    diff_means = diff_means([1:i-1,i+1:k]);
    sum_var = post_vars(i) + post_vars;
    sum_var = sum_var([1:i-1,i+1:k]);
    
    % Argument for Psi_nu[s] function in Branke et al. (2005)
    s = diff_means./sqrt(sum_var); % a (1 x k-1) vector

    if strcmp(known_var, 'known') == 1 % known variances --> normal distribution
    
        % Calculate Bonferroni bound for posterior EOC of Alternative i
        post_qty = sum(sqrt(sum_var).*(normpdf(s) - s.*normcdf(-s)));
    
    elseif strcmp(known_var, 'unknown') == 1 % unknown variances --> t distribution
    
        % Calculate df of the differences W_i - W_j for all j != i
        % Welch approximation
        df_denominator = post_vars(i)^2/post_dfs(i) + post_vars.^2./post_dfs;
        df_denominator = df_denominator([1:i-1,i+1:k]);
        diff_df = (sum_var.^2)./df_denominator;

        % Calculate Bonferroni bound for posterior EOC of Alternative i
        post_qty = sum(sqrt(sum_var).*(((diff_df + s.^2)./(diff_df - 1)).*tpdf(s, diff_df) - s.*tcdf(-s, diff_df)));

    else
    
        disp('Argument for unknown/known variance should be "known" or "unknown".')
        return
    
    end
    
elseif strcmp(qty_type, 'EOC-Slep') == 1
    
    % Calculate means and variances of the differences W_i - W_j for all j != i
    diff_means = post_means(i) - post_means;
    diff_means = diff_means([1:i-1,i+1:k]);
    sum_var = post_vars(i) + post_vars;
    sum_var = sum_var([1:i-1,i+1:k]);
    
    if strcmp(known_var, 'known') == 1 % known variances --> normal distribution
        
        % Calculate Slepian bound for posterior EOC of Alternative i
        SlepPBS = @(delta) 1 - prod(normcdf((delta + diff_means)./sqrt(sum_var)));
        post_qty = integral(@(delta)SlepPBS(delta), 0, Inf, 'ArrayValued', true); 

    elseif strcmp(known_var, 'unknown') == 1 % unknown variances --> t distribution
                
        % Calculate df of the differences W_i - W_j for all j != i
        % Welch approximation
        df_denominator = post_vars(i)^2/post_dfs(i) + post_vars.^2./post_dfs;
        df_denominator = df_denominator([1:i-1,i+1:k]);
        diff_df = (sum_var.^2)./df_denominator;
        
        % Calculate Slepian bound for posterior EOC of Alternative i
        SlepPBS = @(delta) 1 - prod(tcdf((delta + diff_means)./sqrt(sum_var), diff_df));
        post_qty = integral(@(delta)SlepPBS(delta), 0, Inf, 'ArrayValued', true); 
    
    else
    
        disp('Argument for unknown/known variance should be "known" or "unknown".')
        return
    
    end
        
elseif strcmp(qty_type, 'EOC-Exact1') == 1

    % Sum of k 1d integrals

    if strcmp(known_var, 'known') == 1 % known variances --> normal distribution

        E_Wk = 0;

        for j = 1:k

            % Exclude Alternative j
            minus_means = post_means([1:j-1,j+1:k]);
            minus_var = post_vars([1:j-1,j+1:k]);

            prodfun = @(w) w.*prod(normcdf((w - minus_means)./sqrt(minus_var))).*normpdf(w, post_means(j), sqrt(post_vars(j)));
            j_int = integral(@(w)prodfun(w), -Inf, Inf, 'ArrayValued', true, 'Waypoints', [post_means(j) - offset_sd*sqrt(post_vars(j)), post_means(j) + offset_sd*sqrt(post_vars(j))]);
            %j_int = integral(@(w)prodfun(w), -Inf, Inf, 'ArrayValued', true);
            E_Wk = E_Wk + j_int;

        end

        % Calculate exact posterior EOC of Alternative i
        post_qty = E_Wk - post_means(i);

    elseif strcmp(known_var, 'unknown') == 1 % unknown variances --> t distribution

        E_Wk = 0;

        for j = 1:k

            % Exclude Alternative j
            minus_means = post_means([1:j-1,j+1:k]);
            minus_var = post_vars([1:j-1,j+1:k]);
            minus_df = post_dfs([1:j-1,j+1:k]);

            % Integrate the product of k-1 3-parameter t cdfs and one 3-parameter t pdf
            prodfun = @(w) w.*prod(tcdf((w - minus_means)./sqrt(minus_var), minus_df)).*(1/sqrt(post_vars(j))).*tpdf((w - post_means(j))/sqrt(post_vars(j)), post_dfs(j));
            j_int = integral(@(w)prodfun(w), -Inf, Inf, 'ArrayValued', true, 'Waypoints', [post_means(j) - offset_sd*sqrt(post_vars(j)*post_dfs(j)/(post_dfs(j) - 2)), post_means(j) + offset_sd*sqrt(post_vars(j)*post_dfs(j)/(post_dfs(j) - 2))]);
            %j_int = integral(@(w)prodfun(w), -Inf, Inf, 'ArrayValued', true);

            E_Wk = E_Wk + j_int;

        end

        % Calculate exact posterior EOC of Alternative i
        post_qty = E_Wk - post_means(i);

    else
        disp('Argument for unknown/known variance should be "known" or "unknown".')
        return

    end
    
elseif strcmp(qty_type, 'EOC-Exact2') == 1
            
    % Double integral

    % Exclude Alternative i
    minus_means = post_means([1:i-1,i+1:k]);
    minus_var = post_vars([1:i-1,i+1:k]);

    if strcmp(known_var, 'known') == 1 % known variances --> normal distribution

        % Calculate exact posterior EOC of Alternative i
        prodfun = @(w, delta) (1 - prod(normcdf((delta + w - minus_means)./sqrt(minus_var)))).*normpdf(w, post_means(i), sqrt(post_vars(i)));
        inner_int = @(delta) integral(@(w)prodfun(w, delta), -Inf, Inf, 'ArrayValued', true, 'Waypoints', [post_means(i) - offset_sd*sqrt(post_vars(i)), post_means(i) + offset_sd*sqrt(post_vars(i))]);
        %inner_int = @(delta) integral(@(w)prodfun(w, delta), -Inf, Inf);
        post_qty = integral(@(delta)inner_int(delta), 0, Inf, 'ArrayValued', true);

    elseif strcmp(known_var, 'unknown') == 1 % unknown variances --> t distribution

        % Exclude Alternative i
        minus_df = post_dfs([1:i-1,i+1:k]);

        % Calculate exact posterior EOC of Alternative i
        prodfun = @(w, delta) (1 - prod(tcdf((delta + w - minus_means)./sqrt(minus_var), minus_df))).*(1/sqrt(post_vars(i))).*tpdf((w - post_means(i))/sqrt(post_vars(i)), post_dfs(i));
        inner_int = @(delta) integral(@(w)prodfun(w, delta), -Inf, Inf, 'ArrayValued', true, 'Waypoints', [post_means(i) - offset_sd*sqrt(post_vars(i)*post_dfs(i)/(post_dfs(i) - 2)), post_means(i) + offset_sd*sqrt(post_vars(i)*post_dfs(i)/(post_dfs(i) - 2))]);
        %inner_int = @(delta) integral(@(w)prodfun(w, delta), -Inf, Inf, 'ArrayValued', true);
        post_qty = integral(@(delta)inner_int(delta), 0, Inf, 'ArrayValued', true);

    end
    
end % end cases

end % end function

