% SS_TO_UCPWGC  Pairwise unconditional GC matrix from SS parameters.
%
% Computes unconditional pairwise GC F(ti,si) = F_{sources(si) -> targets(ti)}
% for specified target and source channels. This is the pairwise scalar
% specialisation of ss_to_ucgc, analogous to how ss_to_pwcgc relates to
% ss_to_mvgc. It is more efficient than calling ss_to_ucgc for each pair
% because this would recompute the target reduced error covariance for
% each source (Step 1 below).
%
% Each entry is computed as:
%
%   F(ti,si) = log Sigma_i^{(i)} - log Sigma_i^{(i,j)}
%
% where Sigma_i^{(r)} is the innovations variance of channel i in the
% marginal model restricted to observing variables r, derived from the
% full joint SS model via the DARE.
%
% IMPORTANT: You MUST run ss_info before using this function!
%
% Input:
%
%   A,C,K,V   - innovations-form state-space parameters (n variables)
%   targets   - (optional) vector of target channel indices. Default: 1:n.
%   sources   - (optional) vector of source channel indices. Default: 1:n.
%
% Output:
%
%   F          - n_targets x n_sources matrix of unconditional pairwise GCs.
%               F(ti,si) = GC from sources(si) to targets(ti), unconditional.
%               Entries where target == source are NaN.
%
% See also: ss_to_ucgc, ss_to_pwcgc, ss_to_mvgc

function F = ss_to_ucpwgc(A,C,K,V,targets,sources)

[n,~,L] = ss_parms(A,C,K,V);

if nargin < 5 || isempty(targets), targets = 1:n; end
if nargin < 6 || isempty(sources), sources = 1:n; end

targets = targets(:)';
sources = sources(:)';

assert(all(targets >= 1 & targets <= n), 'some target indices out of range');
assert(all(sources >= 1 & sources <= n), 'some source indices out of range');

nt = length(targets);
ns = length(sources);

F = nan(nt, ns);

KL = K*L;

% --- Step 1: Univariate marginal innovations variances for each target ---

V_solo = nan(nt, 1);

for ti = 1:nt
    i = targets(ti);
    [~,VRi,rep] = mdare(A,C(i,:),KL*KL',V(i,i),K*V(:,i));
    if sserror(rep,i), continue; end
    V_solo(ti) = VRi;  % scalar
end

% --- Step 2: Bivariate marginal innovations variances for each pair ------

for ti = 1:nt
    if isnan(V_solo(ti)), continue; end
    i = targets(ti);
    
    for si = 1:ns
        j = sources(si);
        if i == j, continue; end
        
        r = [i, j];
        [~,VRij,rep] = mdare(A,C(r,:),KL*KL',V(r,r),K*V(:,r));
        if sserror(rep), continue; end
        
        % VRij(1,1) = innovations variance of i when observing {i, j}
        F(ti,si) = log(V_solo(ti)) - log(VRij(1,1));
    end
end
