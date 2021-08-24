function CandidateIndex = FindMaxPGSCandidatesSort(delta, k_brack, post_means, post_vars)
% Return a vector of indexes for systems that could have the max posterior
% PGS

% INPUTS
% -------------------------------------------------------------------------
% delta         : good selection parameter
% k_brack       : index of alternative with the highest posterior mean
% post_means    : a (1 x k) vector of posterior means for W_1, ..., W_k
% post_vars     : a (1 x k) vector of posterior variances for W_1, ..., W_k

% OUTPUTS
% -------------------------------------------------------------------------
% CandidateIndex    : indexes of candidate alternatives

% MAIN
% -------------------------------------------------------------------------

% First eliminate systems with posterior means at least delta below the
% best
CandidateIndex = find(post_means > post_means(k_brack) - delta);
NumCandidates = length(CandidateIndex);

% Then eliminate systems that are "softly" dominated in terms of posterior
% means and variances
CandMeans = post_means(CandidateIndex);
CandVar = post_vars(CandidateIndex);

% Sort the candidate systems by their posterior means
[SortedMeans, SortedIndex] = sort(CandMeans, 'descend');
CandidateIndex = CandidateIndex(SortedIndex);
SortedVariances = CandVar(SortedIndex);

KeepVector = ones(1, NumCandidates);

for i = 1:NumCandidates
    if KeepVector(i) == 1 % if System i is still in contention
        for j = i:NumCandidates
            if ((SortedMeans(i) - delta/2 > SortedMeans(j)) && (SortedVariances(i) >= SortedVariances(j))) || ((SortedMeans(i) - delta/2 >= SortedMeans(j)) && (SortedVariances(i) > SortedVariances(j)))
                KeepVector(j) = 0; % System j is dominated by some System i
            end        
        end        
    end
end

CandidateIndex = CandidateIndex(logical(KeepVector));

end

