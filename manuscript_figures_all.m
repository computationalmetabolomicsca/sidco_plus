% ======================================================================================
% Matlab code used for creation of figures in the manuscript entitled: 
% Unbiased distance correlation with the consideration of sample size dependence as a data driven network 
% comparative analysis methodology 
% by: Miroslava Cuperlovic-Culf, Anuradha Surendra, Irina Alecu, Finn Archinuk, Hosna Jabbari
%
% MCuperlovic-Culf, Ottawa, 2026
% ======================================================================================
% Step 1 — generate the data  MC simulations
r1 = mc_functional_dependence('n', 1000, 'noise_std', 1, 'seed', 42, 'plot', false);

% Step 2 — run the focused dCor analysis (produces all 4 figures automatically)
results = dcor_convergence_analysis(r1, ...
    'sample_sizes', [5 10 20 50 100 200 500], ...
    'n_reps',       200, ...
    'B_bb',         100, ...
    'alpha',        0.003);



function results = dcor_convergence_analysis(r1, varargin)
% DCOR_CONVERGENCE_ANALYSIS
%   For each functional model in r1, draws random subsamples of increasing
%   size from N=200 simulated observations, computes the UNBIASED distance
%   correlation (U-statistic, Szekely & Rizzo 2014) and three confidence
%   intervals, repeating n_reps times to characterise the sampling
%   distribution and produces figures for the manuscript.
%
% =========================================================================
% COMPUTATION RESULTS
% =========================================================================
%   For each model m, sample size n, replicate r:
%     1. Draw ns random observations from the N=1000 full dataset.
%     2. Compute dCor(U) — unbiased U-statistic.
%     3. Compute p-value via Szekely-Rizzo (2013) t-test:
%           T_n = sqrt(n(n-3)/2) * R / sqrt(1-R^2)  ~  t(n(n-3)/2 - 1)
%     4. Compute three CIs on the replicate distribution:
%           Fisher z   — atanh transform, SE = 1/sqrt(n-3)
%           Bayes BB   — Dirichlet-weighted bootstrap (Rubin 1981)
%           Bernstein  — Maurer & Pontil (2009) empirical bound
%     5. Apply BH-FDR correction across all n_reps replicates at each
%        sample size independently.
%
% =========================================================================
%  FIGURES 
% =========================================================================
%   P-VALUE / POWER (one panel per model)
%     Distribution of raw p-values as violins vs sample size.
%     Empirical power (fraction p < alpha) annotated above each violin.
%     Fraction surviving BH-FDR annotated below.
%
%   VIOLIN + CI BANDS (one panel per model)
%     Violin of dCor(U) distribution vs sample size.
%     Three CI bands overlaid: Fisher z (black), Bayes BB (blue),
%     Bernstein (green).  Red dashed = full-N reference.
%     Individual replicate dots coloured by BH-FDR q-value:
%       green = q < alpha (significant after FDR)
%       red   = q >= alpha (not significant after FDR)
%
%
%   CI COMPARISON (one panel per model)
%     Width of each CI method vs sample size on log scale.
%     Shows how quickly each CI narrows with increasing n.
%
%   SAMPLE SIZE EFFECT SUMMARY
%     Bias (mean dCor − true dCor) and SD vs log(n), all models.
%
% =========================================================================
%  USAGE
% =========================================================================
%   results = dcor_convergence_analysis(r1)
%   results = dcor_convergence_analysis(r1, ...
%       'sample_sizes', [5 10 20 50 100 200 500], ...
%       'n_reps', 200, 'B_bb', 500, 'alpha', 0.05)
%
% INPUTS:
%   r1  — struct from mc_functional_dependence() with fields:
%          .linear / .quadratic / .cubic / .sin / .random
%          each containing .X, .Y  [N x 1]
%
% OPTIONAL:
%   'sample_sizes'  — subset sizes to test       (default: [5 10 20 50 100 200 500])
%   'n_reps'        — resampling replicates      (default: 200)
%   'B_bb'          — Bayesian bootstrap reps    (default: 500)
%   'alpha'         — significance/CI level      (default: 0.05)
%   'seed'          — random seed                (default: 0)
%
% OUTPUTS:
%   results  — struct with fields per model:
%     .dcor_u          [n_reps x n_sizes]  raw dCor U-stat values
%     .pval            [n_reps x n_sizes]  raw SR13 p-values
%     .qval            [n_reps x n_sizes]  BH-FDR corrected q-values
%     .fdr_sig         [n_reps x n_sizes]  logical: q < alpha
%     .ci_fz           [n_sizes x 2]       Fisher z CI [lo hi]
%     .ci_bb           [n_sizes x 2]       Bayes BB CI
%     .ci_bern         [n_sizes x 2]       Bernstein CI
%     .ref_dcor        scalar              full-N dCor U-stat
%     .ref_pval        scalar              full-N p-value
%   results.params     — input parameters
%
% REFERENCES:
%   Szekely & Rizzo (2014). Ann. Statist. 42(6) — U-statistic dCor
%   Szekely & Rizzo (2013). J. Multivariate Anal. 117 — t-test
%   Maurer & Pontil (2009). ICML — Bernstein CI
%   Rubin (1981). Ann. Statist. 9(1) — Bayesian bootstrap
%   Benjamini & Hochberg (1995). J. R. Stat. Soc. B — BH-FDR


ip = inputParser;
addRequired(ip,  'r1');
addParameter(ip, 'sample_sizes', [5 10 20 50 100 200 500]);
addParameter(ip, 'n_reps',       200);
addParameter(ip, 'B_bb',         500);
addParameter(ip, 'B_perm',       499);
addParameter(ip, 'perm_max_n',   100);
addParameter(ip, 'alpha',        0.05);
addParameter(ip, 'seed',         0);
parse(ip, r1, varargin{:});
pm = ip.Results;

rng(pm.seed);

sample_sizes = pm.sample_sizes(:)';
n_sizes      = numel(sample_sizes);
n_reps       = pm.n_reps;
alpha        = pm.alpha;
B_bb         = pm.B_bb;
B_perm       = pm.B_perm;
perm_max_n   = pm.perm_max_n;

model_names  = {'linear','quadratic','cubic','sin','random'};
model_labels = {'Linear','Quadratic','Cubic','Sin','Random'};
n_models     = numel(model_names);

fprintf('\n=================================================================\n');
fprintf('  dCor U-stat Convergence Analysis\n');
fprintf('  Sample sizes : %s\n', num2str(sample_sizes));
fprintf('  Replicates   : %d    B_bb=%d    alpha=%.3f\n', n_reps, B_bb, alpha);
fprintf('  CI methods   : Fisher z / Bayes BB / Bernstein\n');
fprintf('  P-value      : Szekely-Rizzo (2013) t-test\n');
fprintf('  FDR          : Benjamini-Hochberg per sample size\n');
fprintf('=================================================================\n\n');

results = struct();

for m = 1:n_models
    mname  = model_names{m};
    X_full = r1.(mname).X(:);
    Y_full = r1.(mname).Y(:);
    N_full = numel(X_full);

    fprintf('Model: %-12s  N=%d\n', model_labels{m}, N_full);

    % Full-N reference
    ref_dc = dcor_U_stat(X_full, Y_full, N_full);
    ref_pv = pval_sr13(ref_dc, N_full);
    ref_pval_chi2 = dcor_chi2_test(ref_dc,N_full);
    ref_pval_perm = pval_permutation(X_full, Y_full,ref_dc,N_full,B_perm);  
   

    dc_mat    = NaN(n_reps, n_sizes);
    pval_mat  = NaN(n_reps, n_sizes);   % SR13 t-test
    pval_perm = NaN(n_reps, n_sizes);   % permutation
    pval_chi2 = NaN(n_reps, n_sizes);   % chi2
    qval_mat  = NaN(n_reps, n_sizes);
    fdr_mat   = false(n_reps, n_sizes);

    ci_fz   = NaN(n_sizes, 2);
    ci_bb   = NaN(n_sizes, 2);
    ci_bern = NaN(n_sizes, 2);

    for s = 1:n_sizes
        ns = sample_sizes(s);
        if ns >= N_full
            warning('n=%d >= N_full=%d for %s — skipped.', ns, N_full, mname);
            continue;
        end
        if ns < 5
            % dCor U requires n >= 5
            continue;
        end

        use_perm = (ns <= perm_max_n);
        for r = 1:n_reps
            idx            = randperm(N_full, ns);
            Xs             = X_full(idx);
            Ys             = Y_full(idx);
            dc             = dcor_U_stat(Xs, Ys, ns);
            dc_mat(r, s)   = dc;
            pval_mat(r, s) = pval_sr13(dc, ns);
            pval_chi2(r, s) = dcor_chi2_test(dc,ns);
            if use_perm
                pval_perm(r, s) = pval_permutation(Xs, Ys, ns, dc, B_perm);
            else
                pval_perm(r, s) = pval_chi2(r, s);  % reuse chi-test for large n
            end
        end

        % BH-FDR across the n_reps replicates at this sample size
        [h_col, ~, ~, q_col]=fdr_bh(pval_chi2(:, s),alpha);

        qval_mat(:, s)  = q_col;
        fdr_mat(:, s)   = h_col;

        % CI on the replicate distribution
        d = dc_mat(:, s);
        d = d(~isnan(d));
        if numel(d) >= 2
            ci_fz(s,:)   = fisher_z_ci(median(d), ns, alpha);
            ci_bb(s,:)   = bayes_bb_ci(d, B_bb, alpha);
            ci_bern(s,:) = bernstein_ci(d, ns, alpha);
        end

        fprintf('  n=%4d | mean=%.3f  SD=%.3f  FDR-sig=%.0f%%  p_t=%.3f  p_perm=%.3f\n', ...
            ns, nanmean(dc_mat(:,s)), nanstd(dc_mat(:,s)), ...
            100*mean(fdr_mat(:,s)), ...
            nanmean(pval_mat(:,s)), nanmean(pval_perm(:,s)));
    end

    results.(mname).dcor_u    = dc_mat;
    results.(mname).pval      = pval_mat;    % SR13 t-test
    results.(mname).pval_perm = pval_perm;   % permutation
    results.(mname).pval_chi2 = pval_chi2;
    results.(mname).qval     = qval_mat;
    results.(mname).fdr_sig  = fdr_mat;
    results.(mname).ci_fz    = ci_fz;
    results.(mname).ci_bb    = ci_bb;
    results.(mname).ci_bern  = ci_bern;
    results.(mname).ref_dcor = ref_dc;
    results.(mname).ref_pval = ref_pv;
    fprintf('\n');
end

results.model_names  = model_names;
results.model_labels = model_labels;
results.sample_sizes = sample_sizes;
results.params       = pm;

% -------------------------------------------------------------------------
%    figures
% -------------------------------------------------------------------------
dcor_convergence_plot(results);

fprintf('Done.\n');
end


% =========================================================================
%  dCor U-STATISTIC  (Szekely & Rizzo 2014)
% =========================================================================
function dc = dcor_U_stat(X, Y, n)
X=normalize(X);
Y=normalize(Y);
if n < 5, dc = NaN; return; end
A  = squareform(pdist(X(:)));
B  = squareform(pdist(Y(:)));
Au = u_center(A, n);
Bu = u_center(B, n);
mask     = ~eye(n);
dF       = n * (n-3);
dXY      = sum(Au(:) .* Bu(:) .* mask(:)) / dF;
dXX      = sum(Au(:).^2        .* mask(:)) / dF;
dYY      = sum(Bu(:).^2        .* mask(:)) / dF;
denom    = sqrt(dXX * dYY);
if denom <= 0, dc = 0; return; end
val = dXY / denom;
dc  = abs(max(0, sign(val) * sqrt(abs(val))));
end

function Au = u_center(A, n)
rs = sum(A, 2);  gs = sum(rs);
Au = A - rs*ones(1,n)/(n-2) - ones(n,1)*rs'/(n-2) + gs/((n-1)*(n-2));
Au(1:n+1:end) = 0;
end


% =========================================================================
%  P-VALUE  — Szekely & Rizzo (2013) t-test
% =========================================================================
function p = pval_sr13(dc, n)
if n < 5 || isnan(dc), p = NaN; return; end
dc = max(0, min(0.9999, dc));
nu = n*(n-3)/2 - 1;
if nu <= 0, p = NaN; return; end
Tn = sqrt(nu) * dc / sqrt(max(1e-12, 1 - dc^2));
p  = 1- tcdf((Tn), nu);
end


% =========================================================================
%  CI METHODS
% =========================================================================
function ci = fisher_z_ci(r, n, alpha)
if n <= 3 || isnan(r), ci = [NaN NaN]; return; end
r  = max(-0.9999, min(0.9999, r));
z  = atanh(r);
se = 1 / sqrt(n - 3);
zc = norminv(1 - alpha/2);
ci = max(-1, min(1, tanh([z - zc*se,  z + zc*se])));
end

function ci = bayes_bb_ci(d, B_bb, alpha)
d = d(:);
K = numel(d);
if K < 2, ci = [NaN NaN]; return; end
G = -log(rand(B_bb, K));
W = G ./ sum(G, 2);
ci = max(-1, min(1, quantile(W * d, [alpha/2,  1-alpha/2])));
end

function ci = bernstein_ci(d, ns, alpha)
% Maurer & Pontil (2009)
d     = d(:);
d_bar = mean(d);
s2    = mean((d - d_bar).^2);
K     = ns;
if K < 2, ci = [NaN NaN]; return; end
lnt  = log(2 / alpha);
eps  = sqrt(2 * s2 * lnt / K) + (2/3) * lnt / K;
ci   = max(-1, min(1, [d_bar - eps,  d_bar + eps]));
end


% =========================================================================
%  PERMUTATION P-VALUE for dCor U-stat
% =========================================================================
function p = pval_permutation(X, Y, n, dc_obs, B_perm)
if n < 5 || isnan(dc_obs) || (dc_obs==0) || B_perm < 1, p = NaN; return; end
count = 0;
for b = 1:B_perm
    Yp = Y(randperm(n));
    if dcor_U_stat(X, Yp, n) >= dc_obs
        count = count + 1;
    end
end
p = (1 + count) / (1 + B_perm);
end
% =========================================================================
%  CHI-SQUARE P-VALUE for dCor U-stat
% =========================================================================
function [pval] = dcor_chi2_test(dcor,n)
% Distance correlation chi-square test 
% X: n×p   Y: n×q   (n samples)
% Returns:
%   dcor : distance correlation
%   T    : test statistic (n * dCor^2)
%   pval : chi-square mixture approximation p-value

    % --- Test statistic ---
    T = n * dcor^2;
    k = 1/2;      % shape
    theta = 2;    % scale

    pval = 1 - chi2cdf(T+1,1); %gamcdf(T, k, theta);
end


%=========================================================================
%  BH-FDR  (Benjamini & Hochberg 1995)
%  Applied across the n_reps replicates at a single sample size.
% =========================================================================
% fdr_bh() - Executes the Benjamini & Hochberg (1995) and the Benjamini &
%            Yekutieli (2001) procedure for controlling the false discovery 
%            rate (FDR) of a family of hypothesis tests. 
%
%            Also return the false coverage-statement rate 
%            (FCR)-adjusted selected confidence interval coverage (i.e.,
%            the coverage needed to construct multiple comparison corrected
%            confidence intervals that correspond to the FDR-adjusted p-values).
%
%
% Usage:
%  >> [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pvals,q,method,report);
%
% Required Input:
%   pvals - A vector or matrix (two dimensions or more) containing the
%           p-value of each individual test in a family of tests.
%
% Optional Inputs:
%   q       - The desired false discovery rate. {default: 0.05}
%   method  - ['pdep' or 'dep'] If 'pdep,' the original Bejnamini & Hochberg
%             FDR procedure is used, which is guaranteed to be accurate if
%             the individual tests are independent or positively dependent
%             (e.g., Gaussian variables that are positively correlated or
%             independent).  If 'dep,' the FDR procedure
%             described in Benjamini & Yekutieli (2001) that is guaranteed
%             to be accurate for any test dependency structure (e.g.,
%             Gaussian variables with any covariance matrix) is used. 'dep'
%             is always appropriate to use but is less powerful than 'pdep.'
%             {default: 'pdep'}
%   report  - ['yes' or 'no'] If 'yes', a brief summary of FDR results are
%             output to the MATLAB command line {default: 'no'}
%
%
% Outputs:
%   h       - A binary vector or matrix of the same size as the input "pvals."
%             If the ith element of h is 1, then the test that produced the 
%             ith p-value in pvals is significant (i.e., the null hypothesis
%             of the test is rejected).
%   crit_p  - All uncorrected p-values less than or equal to crit_p are 
%             significant (i.e., their null hypotheses are rejected).  If 
%             no p-values are significant, crit_p=0.
%   adj_ci_cvrg - The FCR-adjusted BH- or BY-selected 
%             confidence interval coverage. For any p-values that 
%             are significant after FDR adjustment, this gives you the
%             proportion of coverage (e.g., 0.99) you should use when generating
%             confidence intervals for those parameters. In other words,
%             this allows you to correct your confidence intervals for
%             multiple comparisons. You can NOT obtain confidence intervals 
%             for non-significant p-values. The adjusted confidence intervals
%             guarantee that the expected FCR is less than or equal to q
%             if using the appropriate FDR control algorithm for the  
%             dependency structure of your data (Benjamini & Yekutieli, 2005).
%             FCR (i.e., false coverage-statement rate) is the proportion 
%             of confidence intervals you construct
%             that miss the true value of the parameter. adj_ci=NaN if no
%             p-values are significant after adjustment.
%   adj_p   - All adjusted p-values less than or equal to q are significant
%             (i.e., their null hypotheses are rejected). Note, adjusted 
%             p-values can be greater than 1.
%
%
% References:
%   Benjamini, Y. & Hochberg, Y. (1995) Controlling the false discovery
%     rate: A practical and powerful approach to multiple testing. Journal
%     of the Royal Statistical Society, Series B (Methodological). 57(1),
%     289-300.
%
%   Benjamini, Y. & Yekutieli, D. (2001) The control of the false discovery
%     rate in multiple testing under dependency. The Annals of Statistics.
%     29(4), 1165-1188.
%
%   Benjamini, Y., & Yekutieli, D. (2005). False discovery rate?adjusted 
%     multiple confidence intervals for selected parameters. Journal of the 
%     American Statistical Association, 100(469), 71?81. doi:10.1198/016214504000001907
%
%
% Example:
%  nullVars=randn(12,15);
%  [~, p_null]=ttest(nullVars); %15 tests where the null hypothesis
%  %is true
%  effectVars=randn(12,5)+1;
%  [~, p_effect]=ttest(effectVars); %5 tests where the null
%  %hypothesis is false
%  [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh([p_null p_effect],.05,'pdep','yes');
%  data=[nullVars effectVars];
%  fcr_adj_cis=NaN*zeros(2,20); %initialize confidence interval bounds to NaN
%  if ~isnan(adj_ci_cvrg),
%     sigIds=find(h);
%     fcr_adj_cis(:,sigIds)=tCIs(data(:,sigIds),adj_ci_cvrg); % tCIs.m is available on the
%     %Mathworks File Exchagne
%  end
%
%
% Updated from:
% Author:
% David M. Groppe
% Kutaslab
% Dept. of Cognitive Science
% University of California, San Diego
% March 24, 2010
%%%%%%%%%%%%%%%% REVISION LOG %%%%%%%%%%%%%%%%%
%
% 5/7/2010-Added FDR adjusted p-values
% 5/14/2013- D.H.J. Poot, Erasmus MC, improved run-time complexity
% 10/2015- Now returns FCR adjusted confidence intervals

function [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pvals,q,method,report)
if nargin<1
    error('You need to provide a vector or matrix of p-values.');
else
    if ~isempty(find(pvals<0,1))
        error('Some p-values are less than 0.');
    elseif ~isempty(find(pvals>1,1))
        error('Some p-values are greater than 1.');
    end
end
if nargin<2
    q=.05;
end
if nargin<3
    method='pdep';
end
if nargin<4
    report='no';
end
s=size(pvals);
if (length(s)>2) || s(1)>1
    [p_sorted, sort_ids]=sort(reshape(pvals,1,prod(s)));
else
    %p-values are already a row vector
    [p_sorted, sort_ids]=sort(pvals);
end
[dummy, unsort_ids]=sort(sort_ids); %indexes to return p_sorted to pvals order
m=length(p_sorted); %number of tests
if strcmpi(method,'pdep')
    %BH procedure for independence or positive dependence
    thresh=(1:m)*q/m;
    wtd_p=m*p_sorted./(1:m);
    
elseif strcmpi(method,'dep')
    %BH procedure for any dependency structure
    denom=m*sum(1./(1:m));
    thresh=(1:m)*q/denom;
    wtd_p=denom*p_sorted./[1:m];
    %Note, it can produce adjusted p-values greater than 1!
    %compute adjusted p-values
else
    error('Argument ''method'' needs to be ''pdep'' or ''dep''.');
end
if nargout>3
    %compute adjusted p-values; This can be a bit computationally intensive
    adj_p=zeros(1,m)*NaN;
    [wtd_p_sorted, wtd_p_sindex] = sort( wtd_p );
    nextfill = 1;
    for k = 1 : m
        if wtd_p_sindex(k)>=nextfill
            adj_p(nextfill:wtd_p_sindex(k)) = wtd_p_sorted(k);
            nextfill = wtd_p_sindex(k)+1;
            if nextfill>m
                break;
            end
        end
    end
    adj_p=reshape(adj_p(unsort_ids),s);
end
rej=p_sorted<=thresh;
max_id=find(rej,1,'last'); %find greatest significant pvalue
if isempty(max_id)
    crit_p=0;
    h=pvals*0;
    adj_ci_cvrg=NaN;
else
    crit_p=p_sorted(max_id);
    h=pvals<=crit_p;
    adj_ci_cvrg=1-thresh(max_id);
end
if strcmpi(report,'yes')
    n_sig=sum(p_sorted<=crit_p);
    if n_sig==1
        fprintf('Out of %d tests, %d is significant using a false discovery rate of %f.\n',m,n_sig,q);
    else
        fprintf('Out of %d tests, %d are significant using a false discovery rate of %f.\n',m,n_sig,q);
    end
    if strcmpi(method,'pdep')
        fprintf('FDR/FCR procedure used is guaranteed valid for independent or positively dependent tests.\n');
    else
        fprintf('FDR/FCR procedure used is guaranteed valid for independent or dependent tests.\n');
    end
end
end

function dcor_convergence_plot(results, varargin)
% DCOR_CONVERGENCE_PLOT
%   Publication-ready figures for the dCor U-stat convergence analysis.
%   Called automatically by dcor_convergence_analysis(), or separately:
%
%     dcor_convergence_plot(results)
%     dcor_convergence_plot(results, 'export', true, 'fmt', 'pdf')
%
% OPTIONAL NAME-VALUE PAIRS:
%   'export'  — save figures (default: false)
%   'fmt'     — 'pdf'|'png'|'svg'  (default: 'pdf')
%   'outdir'  — output directory   (default: pwd)
%   'names'   — group labels for Set1/Set2 in figure titles (default: {'','...'})
%
% FIGURES:
%   Fig 1 — Violin distributions + three CI bands + FDR-coloured dots
%   Fig 2 — P-value / power distributions  (raw p and FDR-corrected)
%   Fig 3 — CI width comparison vs sample size
%   Fig 4 — Bias and variance summary

ip = inputParser;
addRequired(ip, 'results');
addParameter(ip, 'export', false);
addParameter(ip, 'fmt',    'pdf');
addParameter(ip, 'outdir', pwd);
parse(ip, results, varargin{:});
pm = ip.Results;

% ------------------------------------------------------------------
% Figure properties
% ------------------------------------------------------------------
S.font      = 'Helvetica';
S.fsz       = 8;          % base font
S.fsz_title = 9;
S.fsz_ax    = 7;
S.lw        = 0.85;
S.lw_ax     = 0.6;
S.ms        = 3.5;

% Wong (2011) colourblind-safe
S.blue   = [  0 114 178]/255;
S.orange = [230 159   0]/255;
S.green  = [  0 158 115]/255;
S.red    = [213  94   0]/255;
S.purple = [204 121 167]/255;
S.sky    = [ 86 180 233]/255;
S.black  = [  0   0   0];
S.grey   = [0.50 0.50 0.50];
S.lgrey  = [0.88 0.88 0.88];

% Per-model colours
S.mcols = [  0 114 178;
           213  94   0;
             0 158 115;
           204 121 167;
           120 120 120] / 255;

% CI colours
S.fz_col    = S.black;          S.fz_fill   = [0.80 0.80 0.80];
S.bb_col    = S.blue;           S.bb_fill   = [0.72 0.87 0.96];
S.bern_col  = S.green;          S.bern_fill = [0.70 0.93 0.82];

% FDR significance dot colours
S.col_sig    = S.green;    % q < alpha — significant
S.col_nonsig = S.red;      % q >= alpha — not significant

S.export = pm.export;
S.fmt    = pm.fmt;
S.outdir = pm.outdir;
S.params = results.params;   % used by plot_pval_comparison, plot_ci_coverage

model_names  = results.model_names;
model_labels = results.model_labels;
sample_sizes = results.sample_sizes;
n_models     = numel(model_names);
n_sizes      = numel(sample_sizes);
alpha        = results.params.alpha;
n_reps       = results.params.n_reps;
x_pos        = 1:n_sizes;
x_tick       = arrayfun(@num2str, sample_sizes, 'UniformOutput', false);
Y_LIM        = [-0.05, 1.05];

fig1 = plot_violins_ci(results, S, model_names, model_labels, ...
    sample_sizes, x_pos, x_tick, Y_LIM, alpha, n_reps, n_sizes);
fig2 = plot_pvalues(results, S, model_names, model_labels, ...
    sample_sizes, x_pos, x_tick, alpha, n_reps, n_sizes);
fig3 = plot_ci_widthsA(results, S, model_names, model_labels, ...
    sample_sizes, alpha, n_reps, n_sizes);

fig4 = plot_bias_sd(results, S, model_names, model_labels, ...
    sample_sizes, n_sizes);
fig5 = plot_pval_comparison(results, S, model_names, model_labels, ...
    sample_sizes, n_sizes);

fig5A = plot_pval_comparisonA(results, S, model_names, model_labels, ...
    sample_sizes, n_sizes,0);

fig5B = plot_pval_comparisonA(results, S, model_names, model_labels, ...
    sample_sizes, n_sizes,1);


fig6 = plot_ci_coverage(results, S, model_names, model_labels, ...
    sample_sizes, n_sizes);
fig7 = plot_dcor_diff_significance(results, S, model_names, model_labels, ...
    sample_sizes, n_sizes);

fig8 = plot_power_type1(results, S, model_names, model_labels, ...
               sample_sizes, n_sizes);

if S.export
    save_fig(fig1, 'Fig4_dCor_Violins_CI',      S);
    save_fig(fig2, 'Fig4_dCor_Pvalues',         S);
    save_fig(fig3, 'Fig5_dCor_CI_width',        S);
    save_fig(fig4, 'Fig6_dCor_Bias_SD',         S);
    save_fig(fig5, 'Fig4_dCor_Pval_Compare',    S);
    save_fig(fig5A, 'Fig4A_dCor_Pval_Compare',    S);
    save_fig(fig5B, 'Fig4A_dCor_Qval_Compare',    S);

    save_fig(fig6, 'Fig5_dCor_CI_Coverage',     S);
    save_fig(fig7, 'Fig5_dCor_Diff_Signif',     S);
    save_fig(fig8, 'Fig8_Power_Type1', S);


end
end


% =========================================================================
%  FIGURE 1 — VIOLIN + CI BANDS + FDR-COLOURED DOTS
%  Rows = models, single column.  Each panel: dCor(U) distribution vs n.
%  Violin shape from KSD. Three CI bands overlaid.
%  Individual replicate dots: green = FDR q<alpha, red = q>=alpha.
% =========================================================================
function fig = plot_violins_ci(res, S, model_names, model_labels, ...
        sample_sizes, x_pos, x_tick, Y_LIM, alpha, n_reps, n_sizes)

n_models = numel(model_names);
ci_pct   = 100 * (1 - alpha);

fig = figure('Name','Fig1_dCor_Violins_CI','NumberTitle','off', ...
    'Color','w','Units','centimeters','Position',[1 1 22 28]);

sgtitle(sprintf(['dCor (U-stat) Sampling Distribution vs Sample Size  |  ' ...
    '%d replicates  |  %.0f%% CI'], n_reps, ci_pct), ...
    'FontName',S.font,'FontSize',S.fsz_title,'FontWeight','bold','Color',S.black);

panel_letters = 'ABCDE';

for m = 1:n_models
    ax = subplot(n_models, 1, m);
    hold(ax, 'on');

    mname = model_names{m};
    mcol  = S.mcols(m,:);
    ref   = res.(mname).ref_dcor;

    % --- CI bands (drawn first as background) ---
    draw_ci_band(ax, x_pos, res.(mname).ci_fz,   S.fz_fill,   S.fz_col,   1.1, 0.22, '-');
    draw_ci_band(ax, x_pos, res.(mname).ci_bb,   S.bb_fill,   S.bb_col,   1.3, 0.26, '--');
    draw_ci_band(ax, x_pos, res.(mname).ci_bern, S.bern_fill, S.bern_col, 1.3, 0.30, '-.');

    % --- Violin + FDR dots ---
    for s = 1:n_sizes
        d    = res.(mname).dcor_u(:, s);
        q    = res.(mname).qval(:, s);
        sig  = res.(mname).fdr_sig(:, s);
        ok   = ~isnan(d);
        if sum(ok) < 3, continue; end
        d_ok = d(ok);  q_ok = q(ok);  sig_ok = sig(ok);

        % Violin shape
        d_c = max(Y_LIM(1), min(Y_LIM(2), d_ok));
        try
            [f, xi] = ksdensity(d_c, 'Support', [Y_LIM(1)-0.01, Y_LIM(2)+0.01], ...
                'BoundaryCorrection', 'reflection', 'NumPoints', 120);
        catch
            continue
        end
        in_r = xi >= Y_LIM(1) & xi <= Y_LIM(2);
        f = f(in_r);  xi = xi(in_r);
        if isempty(f) || max(f) == 0, continue; end
        f = f / max(f) * 0.28;
        patch(ax, [x_pos(s)+f, x_pos(s)-fliplr(f)], [xi, fliplr(xi)], mcol, ...
            'FaceAlpha', 0.18, 'EdgeColor', mcol, 'LineWidth', 0.5);

        % Median line and IQR bar
        med = median(d_ok);
        plot(ax, [x_pos(s)-0.15, x_pos(s)+0.15], [med, med], '-', ...
            'Color', mcol*0.50, 'LineWidth', 1.8);
        plot(ax, [x_pos(s), x_pos(s)], [prctile(d_ok,25), prctile(d_ok,75)], ...
            '-', 'Color', mcol*0.60, 'LineWidth', 2.5);

        % FDR-coloured individual dots: green = q<alpha, red = q>=alpha
        sig_mask = ~isnan(q_ok) & (q_ok < alpha);

        % Jitter x slightly so stacked points are visible
        jit = (rand(sum(ok), 1) - 0.5) * 0.18;

        % Draw non-significant (red) first so significant (green) sits on top
        if any(~sig_mask)
            scatter(ax, x_pos(s) + jit(~sig_mask), d_ok(~sig_mask), 6, ...
                S.col_nonsig, 'filled', ...
                'MarkerFaceAlpha', 0.35, 'MarkerEdgeColor', 'none', ...
                'HandleVisibility', 'off');
        end
        if any(sig_mask)
            scatter(ax, x_pos(s) + jit(sig_mask), d_ok(sig_mask), 6, ...
                S.col_sig, 'filled', ...
                'MarkerFaceAlpha', 0.55, 'MarkerEdgeColor', 'none', ...
                'HandleVisibility', 'off');
        end
    end

    % Full-N reference line
    yline(ax, ref, '--', 'Color', [0.80 0.05 0.05], 'LineWidth', 1.4, ...
        'Label', sprintf('N=1000  R=%.3f', ref), ...
        'LabelHorizontalAlignment', 'right', ...
        'FontName', S.font, 'FontSize', S.fsz_ax, ...
        'HandleVisibility', 'off');
    yline(ax, 0, ':', 'Color', S.grey, 'LineWidth', 0.6, 'HandleVisibility', 'off');

    ylim(ax, Y_LIM);
    xlim(ax, [0.4, n_sizes + 0.6]);
    set(ax, 'XTick', x_pos, 'XTickLabel', x_tick);
    xlabel(ax, 'Sample size  n', 'FontName', S.font, 'FontSize', S.fsz);
    ylabel(ax, 'dCor (U-stat)', 'FontName', S.font, 'FontSize', S.fsz);
    title(ax, model_labels{m}, 'FontName', S.font, 'FontSize', S.fsz_title, ...
        'FontWeight', 'bold', 'Color', mcol);

    % Random model: light background to distinguish null
    if strcmp(mname, 'random'), ax.Color = [0.96 0.95 0.92]; end

    panel_label(ax, panel_letters(m), S);
    pub_ax(ax, S);
    hold(ax, 'off');
end

% Shared legend annotation
annotation(fig, 'textbox', [0.02 0.001 0.96 0.025], ...
    'String', sprintf(['Black = Fisher z %.0f%% CI  |  Blue dashed = Bayes BB %.0f%% CI  |  ' ...
        'Green dash-dot = Bernstein %.0f%% CI  |  ' ...
        'Dots: green=FDR q<%.2f (significant), red=q\\geq%.2f'], ...
        ci_pct, ci_pct, ci_pct, alpha, alpha), ...
    'FontName', S.font, 'FontSize', 6, 'EdgeColor', 'none', ...
    'Color', S.grey, 'HorizontalAlignment', 'center', 'Interpreter', 'tex');

% q-value colourbar on right
add_qvalue_colourbar(fig, S, alpha);

set(fig, 'PaperUnits', 'centimeters', 'PaperSize', [22 28]);
end


% =========================================================================
%  FIGURE 2 — P-VALUE / POWER
%  Two columns per model: raw p-values (left) and FDR q-values (right).
% =========================================================================
function fig = plot_pvalues(res, S, model_names, model_labels, ...
        sample_sizes, x_pos, x_tick, alpha, n_reps, n_sizes)

n_models = numel(model_names);

fig = figure('Name','Fig2_dCor_Pvalues','NumberTitle','off', ...
    'Color','w','Units','centimeters','Position',[1 1 22 28]);

sgtitle(sprintf(['dCor (U-stat) P-value and FDR q-value Distributions  |  ' ...
    '%d replicates  |  \\alpha=%.2f'], n_reps, alpha), ...
    'FontName', S.font, 'FontSize', S.fsz_title, 'FontWeight', 'bold', ...
    'Color', S.black, 'Interpreter', 'tex');

panel_letters = 'ABCDEFGHIJ';
pl = 0;

for m = 1:n_models
    mname = model_names{m};
    mcol  = S.mcols(m,:);

    for col_idx = 1:2   % 1=raw p, 2=FDR q
        pl  = pl + 1;
        ax  = subplot(n_models, 2, (m-1)*2 + col_idx);
        hold(ax, 'on');

        if col_idx == 1
            data  = res.(mname).pval;
            ylab  = 'Raw p-value';
            titl  = sprintf('%s  —  raw p', model_labels{m});
            col   = mcol;
        else
            data  = res.(mname).qval;
            ylab  = 'BH-FDR  q-value';
            titl  = sprintf('%s  —  BH-FDR q', model_labels{m});
            col   = S.blue * 0.7 + mcol * 0.3;
        end

        for s = 1:n_sizes
            d  = data(:, s);
            ok = ~isnan(d);
            if sum(ok) < 3, continue; end
            d_ok = d(ok);

            % Skip violin if ALL replicates are below threshold —
            % the distribution is entirely significant; nothing informative to show.
            if all(d_ok < alpha)
                % Place a small star to indicate '100% pass'
                text(ax, x_pos(s), alpha * 0.6, '*', ...
                    'HorizontalAlignment', 'center', 'FontSize', 8, ...
                    'Color', S.green, 'Interpreter', 'tex');
                continue;
            end

            % Violin
            d_c = max(0, min(1, d_ok));
            try
                [f, xi] = ksdensity(d_c, 'Support', [-0.01, 1.01], ...
                    'BoundaryCorrection', 'reflection', 'NumPoints', 80);
            catch; continue; end
            in_r = xi >= 0 & xi <= 1;
            f = f(in_r); xi = xi(in_r);
            if isempty(f) || max(f) == 0, continue; end
            f = f / max(f) * 0.28;
            patch(ax, [x_pos(s)+f, x_pos(s)-fliplr(f)], [xi, fliplr(xi)], col, ...
                'FaceAlpha', 0.22, 'EdgeColor', col, 'LineWidth', 0.5);
            med = median(d_ok);
            plot(ax, [x_pos(s)-0.14,x_pos(s)+0.14],[med,med],'-','Color',col*0.50,'LineWidth',1.8);
            plot(ax,[x_pos(s),x_pos(s)],[prctile(d_ok,25),prctile(d_ok,75)], ...
                '-','Color',col*0.60,'LineWidth',2.5);
        end

        % Alpha / significance threshold
        yline(ax, alpha, '--', 'Color', [0.82 0.05 0.05], 'LineWidth', 1.2, ...
            'Label', sprintf('\\alpha=%.2f', alpha), ...
            'LabelHorizontalAlignment', 'right', ...
            'FontName', S.font, 'FontSize', S.fsz_ax, ...
            'HandleVisibility', 'off', 'Interpreter', 'tex');

        % Power / FDR-pass annotations
        for s = 1:n_sizes
            d  = data(:, s);
            ok = ~isnan(d);
            if ~any(ok), continue; end
            frac = mean(d(ok) < alpha);
            txt  = sprintf('%.0f%%', 100*frac);
            % Colour: green>=80%, amber 50-80%, red<50%
            if strcmp(mname,'random')
                tcol = S.grey;
            elseif frac >= 0.80, tcol = S.green;
            elseif frac >= 0.50, tcol = S.orange;
            else,                tcol = S.red;
            end
            text(ax, x_pos(s), 0.95, txt, 'Units', 'data', ...
                'HorizontalAlignment', 'center', ...
                'FontName', S.font, 'FontSize', S.fsz_ax, ...
                'FontWeight', 'bold', 'Color', tcol);
        end

        % Random model: uniform reference lines
        if strcmp(mname, 'random')
            yline(ax, 0.25, ':', 'Color', S.grey, 'LineWidth', 0.7, 'HandleVisibility', 'off');
            yline(ax, 0.75, ':', 'Color', S.grey, 'LineWidth', 0.7, 'HandleVisibility', 'off');
            ax.Color = [0.96 0.95 0.92];
        end

        ylim(ax, [0 1]);  yticks(ax, 0:0.25:1);
        xlim(ax, [0.4, n_sizes + 0.6]);
        set(ax, 'XTick', x_pos, 'XTickLabel', x_tick);
        ylabel(ax, ylab,  'FontName', S.font, 'FontSize', S.fsz);
        if m == n_models
            xlabel(ax, 'Sample size  n', 'FontName', S.font, 'FontSize', S.fsz);
        end
        if m == 1
            title(ax, titl, 'FontName', S.font, 'FontSize', S.fsz_title, ...
                'FontWeight', 'bold', 'Color', col);
        else
            title(ax, model_labels{m}, 'FontName', S.font, ...
                'FontSize', S.fsz_ax + 0.5, 'Color', mcol);
        end

        panel_label(ax, panel_letters(pl), S);
        pub_ax(ax, S);
        hold(ax, 'off');
    end
end

annotation(fig, 'textbox', [0.02 0.001 0.96 0.025], ...
    'String', ['Left: raw p-values (SR13 t-test)  |  Right: BH-FDR corrected q-values  |  ' ...
        'Numbers above violins: fraction < alpha  |  Green>=80%  Amber 50-80%  Red<50%'], ...
    'FontName', S.font, 'FontSize', 6, 'EdgeColor', 'none', ...
    'Color', S.grey, 'HorizontalAlignment', 'center');

set(fig, 'PaperUnits', 'centimeters', 'PaperSize', [22 28]);
end


% =========================================================================
%  FIGURE  — CI WIDTH COMPARISON
%  One panel per model.  Log-x axis.
%  Three lines: Fisher z, Bayes BB, Bernstein — width vs sample size.
%  Shows how quickly each CI narrows and which is tightest/widest.
% =========================================================================
function fig = plot_ci_widths(res, S, model_names, model_labels, ...
        sample_sizes, alpha, n_reps, n_sizes)

n_models = numel(model_names);
ci_pct   = 100 * (1 - alpha);

fig = figure('Name','Fig3_dCor_CI_width','NumberTitle','off', ...
    'Color','w','Units','centimeters','Position',[1 1 24 7]);

sgtitle(sprintf('CI Width vs Sample Size  |  dCor U-stat  |  %.0f%% CI  |  %d replicates', ...
    ci_pct, n_reps), ...
    'FontName', S.font, 'FontSize', S.fsz_title, 'FontWeight', 'bold', 'Color', S.black);

panel_letters = 'ABCDE';

for m = 1:n_models
    ax  = subplot(n_models,1,  m);
    hold(ax, 'on');

    mname = model_names{m};
    ref   = res.(mname).ref_dcor;

    w_fz   = res.(mname).ci_fz(:,2)   - res.(mname).ci_fz(:,1);
    w_bb   = res.(mname).ci_bb(:,2)   - res.(mname).ci_bb(:,1);
    w_bern = res.(mname).ci_bern(:,2) - res.(mname).ci_bern(:,1);

    valid = ~isnan(w_fz) & ~isnan(w_bb) & ~isnan(w_bern);
    if any(valid)
        plot(ax, sample_sizes(valid), w_fz(valid),   '-o',  ...
            'Color', S.fz_col,   'LineWidth', S.lw,   'MarkerSize', S.ms, ...
            'MarkerFaceColor', S.fz_col,   'DisplayName', 'Fisher z');
        plot(ax, sample_sizes(valid), w_bb(valid),   '--s', ...
            'Color', S.bb_col,   'LineWidth', S.lw,   'MarkerSize', S.ms, ...
            'MarkerFaceColor', S.bb_col,   'DisplayName', 'Bayes BB');
        plot(ax, sample_sizes(valid), w_bern(valid), '-.^', ...
            'Color', S.bern_col, 'LineWidth', S.lw,   'MarkerSize', S.ms, ...
            'MarkerFaceColor', S.bern_col, 'DisplayName', 'Bernstein');
    end

    % 1/sqrt(n) reference
    ns_ref = linspace(sample_sizes(1), sample_sizes(end), 100);
    w0     = nanmean([w_fz(find(valid,1)), w_bb(find(valid,1)), w_bern(find(valid,1))]);
    if ~isnan(w0) && w0 > 0
        plot(ax, ns_ref, w0 * sqrt(sample_sizes(find(valid,1))) ./ sqrt(ns_ref), ...
            ':', 'Color', S.grey, 'LineWidth', 0.7, 'DisplayName', '1/\surd n ref');
    end

    % True dCor as horizontal context
    yline(ax, ref, '--', 'Color', [0.80 0.05 0.05], 'LineWidth', 0.8, ...
        'Label', sprintf('R=%.2f', ref), ...
        'LabelHorizontalAlignment', 'right', ...
        'FontName', S.font, 'FontSize', S.fsz_ax - 1, 'HandleVisibility', 'off');

    set(ax, 'XScale', 'log', 'XTick', sample_sizes);
    xlim(ax, [sample_sizes(1)*0.8, sample_sizes(end)*1.3]);
    ylim(ax, [0, Inf]);
    xlabel(ax, 'Sample size  n', 'FontName', S.font, 'FontSize', S.fsz);
    ylabel(ax, 'CI width', 'FontName', S.font, 'FontSize', S.fsz);
    title(ax, model_labels{m}, 'FontName', S.font, 'FontSize', S.fsz_title, ...
        'FontWeight', 'bold', 'Color', S.mcols(m,:));

    if m == 1
        leg = legend(ax, 'Location', 'northeast', 'FontName', S.font, ...
            'FontSize', S.fsz_ax, 'Box', 'off');
    end

    panel_label(ax, panel_letters(m), S);
    pub_ax(ax, S);
    hold(ax, 'off');
end

annotation(fig, 'textbox', [0.02 0.001 0.96 0.030], ...
    'String', ['CI width = upper CI − lower CI  |  log x-axis  |  ' ...
        'Dashed red = true dCor (full N=1000)  |  ' ...
        'Narrower is better — compare methods at same n'], ...
    'FontName', S.font, 'FontSize', 6, 'EdgeColor', 'none', ...
    'Color', S.grey, 'HorizontalAlignment', 'center');

set(fig, 'PaperUnits', 'centimeters', 'PaperSize', [24 7]);
end

% ==========================================================================
% Figure  Alternative way to show effect of CI selection
% =========================================================================
function fig = plot_ci_widthsA(res, S, model_names, model_labels, ...
        sample_sizes, alpha, n_reps, n_sizes)
% Bar chart showing the count of dCor values that fall WITHIN the CI
% interval [lo, hi] for each CI method and sample size.
%
% Count per panel cell = sum(lo <= dCor <= hi) across all pairs/replicates
%
% Rows    = models
% Columns = CI methods (Fisher z, Bayes BB, Bernstein)
% Within each panel: grouped bars, one group per sample size
%
% Data layout:
%   res.(mname).dcor_u     [n_pairs x n_sizes]  observed dCor
%   res.(mname).ci_fz      [n_pairs x n_sizes x 2]  (:,:,1)=lo, (:,:,2)=hi
%   res.(mname).ci_bb      same
%   res.(mname).ci_bern    same

n_models   = numel(model_names);
ci_pct     = 100 * (1 - alpha);
n_methods  = 3;
method_names = {'Fisher z', 'Bayes BB', 'Bernstein'};
method_cols  = [S.fz_col; S.bb_col; S.bern_col];

panel_letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
pl = 0;

bw     = 0.25;   % single bar per group (one method per subplot)
x_base = 1:n_sizes;

fig = figure('Name','Fig3_dCor_within_CI','NumberTitle','off', ...
    'Color','w','Units','centimeters', ...
    'Position',[1 1 6*n_methods 5*n_models]);

sgtitle(sprintf('N dCor values within CI bounds  |  %.0f%% CI  |  %d replicates', ...
    ci_pct, n_reps), ...
    'FontName',S.font,'FontSize',S.fsz_title,'FontWeight','bold', ...
    'Color',S.black);

for m = 1:n_models
    mname = model_names{m};

    % [n_pairs x n_sizes]
    dc      = res.(mname).dcor_u;

    % [n_pairs x n_sizes x 2]
    ci_fz   = res.(mname).ci_fz;
    ci_bb   = res.(mname).ci_bb;
    ci_bern = res.(mname).ci_bern;

    lo_all = {ci_fz(:,1),   ci_bb(:,1),   ci_bern(:,1)};
    hi_all = {ci_fz(:,2),   ci_bb(:,2),   ci_bern(:,2)};

    % Total valid pairs per sample size (denominator for annotation)
    n_valid = zeros(1, n_sizes);
    for s = 1:n_sizes
        n_valid(s) = sum(~isnan(dc(:,s)));
    end

    for tp = 1:n_methods
        pl  = pl + 1;
        ax  = subplot(n_models, n_methods, (m-1)*n_methods + tp);
        hold(ax,'on');

        lo_mat = lo_all{tp};   % [n_pairs x n_sizes]
        hi_mat = hi_all{tp};   % [n_pairs x n_sizes]
        col    = method_cols(tp,:);

        cnt_within  = zeros(1, n_sizes);
        cnt_below   = zeros(1, n_sizes);
        cnt_above   = zeros(1, n_sizes);

        for s = 1:n_sizes
            lo_s  = lo_mat(s);
            hi_s  = hi_mat(s);
            dc_s  = dc(:,s);
            valid = ~isnan(dc_s); %~isnan(lo_s); % & ~isnan(hi_s); % & 

            cnt_within(s) = sum( dc_s(valid) >= lo_s & ...
                                 dc_s(valid) <= hi_s );
            cnt_below(s)  = sum( dc_s(valid) <  lo_s );
            cnt_above(s)  = sum( dc_s(valid) >  hi_s );
        end

        % Stacked bars: within (solid), below (light), above (dark)
        col_within = col;
        col_below  = min(col + 0.35, 1);   % lighter = below lo
        col_above  = col * 0.45;            % darker  = above hi

        b1 = bar(ax, x_base-0.25, cnt_within, bw, ...
            'FaceColor', col_within, 'EdgeColor','none', 'FaceAlpha',0.80, ...
            'DisplayName','Within CI');
        b2 = bar(ax, x_base, cnt_below, bw, ...
            'FaceColor', col_below, 'EdgeColor','none', 'FaceAlpha',0.60, ...
            'DisplayName','Below lo CI');
        b3 = bar(ax, x_base+0.25, cnt_above, bw, ...
            'FaceColor', col_above, 'EdgeColor','none', 'FaceAlpha',0.50, ...
            'DisplayName','Above hi CI');

        % Shift below/above bars slightly so they don't overlap within bars
        b2.XOffset = -0.18;
        b3.XOffset =  0.18;

        % Value labels on within-CI bars
        for s = 1:n_sizes
            text(ax, x_base(s), cnt_within(s) + 0.3, ...
                num2str(cnt_within(s)), ...
                'HorizontalAlignment','center', ...
                'FontName',S.font,'FontSize',max(S.fsz_ax-2,5), ...
                'Color',col*0.65);
        end

        % Total valid reference line
        yline(ax, mean(n_valid), '--', 'Color',[0.5 0.5 0.5], ...
            'LineWidth',0.8, ...
            'Label',sprintf('N valid=%.0f', mean(n_valid)), ...
            'LabelHorizontalAlignment','right', ...
            'FontName',S.font,'FontSize',S.fsz_ax-1, ...
            'HandleVisibility','off');

        % Axes
        ax.XTick               = x_base;
        ax.XTickLabel          = arrayfun(@(n) sprintf('n=%d',n), ...
                                     sample_sizes,'UniformOutput',false);
        ax.XTickLabelRotation  = 45;
        ax.XTickMode           = 'manual';
        ax.XTickLabelMode      = 'manual';
        ax.TickLabelInterpreter = 'none';
        ax.FontSize            = S.fsz_ax;
        ax.XLim                = [0.5, n_sizes + 0.5];
        ax.YLim                = [0, max(n_valid)*1.15 + 1];
        ax.GridAlpha           = 0.13;
        grid(ax,'on'); box(ax,'on');

        ylabel(ax, 'N dCor values', 'FontName',S.font,'FontSize',S.fsz);
        xlabel(ax, 'Sample size',   'FontName',S.font,'FontSize',S.fsz);

        % Titles
        if m == 1
            title(ax, method_names{tp}, ...
                'FontName',S.font,'FontSize',S.fsz_title, ...
                'FontWeight','bold','Color',col*0.75);
        else
            title(ax, sprintf('%s  —  %s', model_labels{m}, method_names{tp}), ...
                'FontName',S.font,'FontSize',S.fsz_title, ...
                'FontWeight','bold','Color',S.mcols(m,:));
        end

        % Legend once only
        if m == 1 && tp == 1
            legend(ax, [b1 b2 b3], ...
                {'Within CI  [lo \leq dCor \leq hi]', ...
                 'Below lo CI  (dCor < lo)', ...
                 'Above hi CI  (dCor > hi)'}, ...
                'Location','northwest','FontName',S.font, ...
                'FontSize',S.fsz_ax,'Box','off','Interpreter','tex');
        end

        if strcmp(mname,'random'), ax.Color = [0.96 0.95 0.92]; end
        panel_label(ax, panel_letters(pl), S);
        pub_ax(ax, S);
        hold(ax,'off');
    end
end

annotation(fig,'textbox',[0.02 0.001 0.96 0.025], ...
    'String',['Solid bar = N dCor values within CI [lo, hi]  |  ' ...
              'Light bar (left) = N values below lower CI bound  |  ' ...
              'Dark bar (right) = N values above upper CI bound  |  ' ...
              'Dashed line = total valid pairs  |  Raw counts, not normalised  |  ' ...
              'Columns = CI method  |  Rows = model'], ...
    'FontName',S.font,'FontSize',6,'EdgeColor','none', ...
    'Color',S.grey,'HorizontalAlignment','center');

set(fig,'PaperUnits','centimeters','PaperSize',[6*n_methods, 5*n_models]);
end
% =========================================================================
%  FIGURE  — BIAS AND SD SUMMARY
%  Left panels: bias (mean dCor − true dCor) vs log(n) for all models.
%  Right panels: SD of dCor vs log(n).
%  All models overlaid in a single panel each.
% =========================================================================
function fig = plot_bias_sd(res, S, model_names, model_labels, ...
        sample_sizes, n_sizes)

n_models = numel(model_names);

fig = figure('Name','Fig4_dCor_Bias_SD','NumberTitle','off', ...
    'Color','w','Units','centimeters','Position',[1 1 16 7]);

sgtitle('dCor (U-stat): Bias and SD vs Sample Size  |  All models', ...
    'FontName', S.font, 'FontSize', S.fsz_title, 'FontWeight', 'bold', ...
    'Color', S.black);

% Left: bias
ax1 = subplot(1, 2, 1);
hold(ax1, 'on');
yline(ax1, 0, '-', 'Color', S.grey, 'LineWidth', 0.7, 'HandleVisibility', 'off');

for m = 1:n_models
    mname = model_names{m};
    ref   = res.(mname).ref_dcor;
    bias  = nanmean(res.(mname).dcor_u, 1) - ref;
    valid = ~isnan(bias);
    if ~any(valid), continue; end
    col = S.mcols(m,:);
    plot(ax1, sample_sizes(valid), bias(valid), '-o', ...
        'Color', col, 'LineWidth', S.lw, 'MarkerSize', S.ms, ...
        'MarkerFaceColor', col, 'DisplayName', model_labels{m});
end

set(ax1, 'XScale', 'log', 'XTick', sample_sizes);
xlim(ax1, [sample_sizes(1)*0.8, sample_sizes(end)*1.3]);
xlabel(ax1, 'Sample size  n', 'FontName', S.font, 'FontSize', S.fsz);
ylabel(ax1, 'Bias  (R_{sub} - R_{true})', 'FontName', S.font, ...
    'FontSize', S.fsz, 'Interpreter', 'tex');
title(ax1, 'Bias', 'FontName', S.font, 'FontSize', S.fsz_title, 'FontWeight', 'bold');
leg = legend(ax1, 'Location', 'northeast', 'FontName', S.font, ...
    'FontSize', S.fsz_ax, 'Box', 'off');
panel_label(ax1, 'A', S);
pub_ax(ax1, S);
hold(ax1, 'off');

% Right: SD
ax2 = subplot(1, 2, 2);
hold(ax2, 'on');

for m = 1:n_models
    mname = model_names{m};
    sd    = nanstd(res.(mname).dcor_u, 0, 1);
    valid = ~isnan(sd);
    if ~any(valid), continue; end
    col = S.mcols(m,:);
    plot(ax2, sample_sizes(valid), sd(valid), '-o', ...
        'Color', col, 'LineWidth', S.lw, 'MarkerSize', S.ms, ...
        'MarkerFaceColor', col, 'DisplayName', model_labels{m});
end

% 1/sqrt(n) reference scaled to first point of first model
m1   = res.(model_names{1});
sd1  = nanstd(m1.dcor_u(:, find(~isnan(nanstd(m1.dcor_u,0,1)), 1)), 0, 1);
ns_ref = linspace(sample_sizes(1), sample_sizes(end), 100);
if ~isnan(sd1) && sd1 > 0
    plot(ax2, ns_ref, sd1 * sqrt(sample_sizes(1)) ./ sqrt(ns_ref), ...
        ':', 'Color', S.grey, 'LineWidth', 0.8, 'DisplayName', '1/\surd n ref');
end

set(ax2, 'XScale', 'log', 'XTick', sample_sizes);
xlim(ax2, [sample_sizes(1)*0.8, sample_sizes(end)*1.3]);
ylim(ax2, [0, Inf]);
xlabel(ax2, 'Sample size  n', 'FontName', S.font, 'FontSize', S.fsz);
ylabel(ax2, 'SD(R_{sub})', 'FontName', S.font, 'FontSize', S.fsz, 'Interpreter', 'tex');
title(ax2, 'Standard deviation', 'FontName', S.font, 'FontSize', S.fsz_title, 'FontWeight', 'bold');
panel_label(ax2, 'B', S);
pub_ax(ax2, S);
hold(ax2, 'off');

annotation(fig, 'textbox', [0.02 0.001 0.96 0.030], ...
    'String', ['Bias = mean(subsample dCor) - full-N dCor  |  ' ...
        'Positive bias = upward overestimation at small n  |  ' ...
        'SD expected to decay as 1/\surd n (dotted)'], ...
    'FontName', S.font, 'FontSize', 6, 'EdgeColor', 'none', ...
    'Color', S.grey, 'HorizontalAlignment', 'center', 'Interpreter', 'tex');

set(fig, 'PaperUnits', 'centimeters', 'PaperSize', [16 7]);
end




% =========================================================================
%  FIGURE  — P-VALUE COMPARISON: t-test vs Permutation
%
%  Layout: n_models rows x 1 column.
%  Each panel: scatter of t-test p-value vs permutation p-value for every
%  replicate at every sample size.  Points coloured by sample size.
%  Diagonal y=x line = perfect agreement.
%  Also shows marginal histograms as insets, and annotates Spearman
%  correlation between the two p-value methods per sample size.
% =========================================================================
function fig = plot_pval_comparison(res, S, model_names, model_labels, ...
        sample_sizes, n_sizes)
% Two-column layout per model:
%   Left  — raw p-values: t-test (x) vs permutation (y)
%   Right — FDR q-values: t-test q (x) vs permutation q (y)
% Both panels share the same quadrant shading and diagonal reference.
% Spearman rho annotated per sample size.

n_models  = numel(model_names);
alpha     = S.params.alpha;
size_cmap = parula(n_sizes);
panel_letters = 'ABCDEFGHIJ';
pl = 0;

fig = figure('Name','Fig5_dCor_Pval_Compare','NumberTitle','off', ...
    'Color','w','Units','centimeters','Position',[1 1 28 32]);

sgtitle(sprintf(['p-value Comparison: SR13 t-test, Chi2 and Permutation  |  \alpha=%.2f'], alpha), ...
    'FontName',S.font,'FontSize',S.fsz_title,'FontWeight','bold','Color',S.black, ...
    'Interpreter','tex');

for m = 1:n_models
    mname = model_names{m};

    % Compute FDR q-values for permutation p-values (per sample size, same
    % as done for t-test q-values in the analysis step).
    qperm_mat = NaN(size(res.(mname).pval_perm));
    for s = 1:n_sizes
        pv = res.(mname).pval_perm(:,s);
        ok = ~isnan(pv);
        if sum(ok) < 2, continue; end
        [~, ~, ~, q_col]=fdr_bh(pv,alpha);
    %     [~, q_col]=mafdr(pv); %MCC &&&

        qperm_mat(:,s) = q_col;
    end

    x_data   = (res.(mname).pval_perm);       % t-test raw p pval

    for col_idx = 1:2   % 1=Chi vs perm raw p, 2=Ttest vs perm
        pl  = pl + 1;
        ax  = subplot(2,5, m+(col_idx-1)*n_models);
        hold(ax,'on');

        if col_idx == 1
            y_data   = (res.(mname).pval_chi2);  % permutation raw p  
            
            xlab     = 'Permutation p-value'; %'Permutation p-value;
            ylab     = 'Chi2  p-value';
            thresh   = (alpha);
            col_title = 'Raw p-value';
        else
            y_data   = (res.(mname).pval);  % permutation raw p  
            
            xlab     = 'Permutation p-value'; %'Permutation p-value;
            ylab     = 't-test  p-value (SR13)';
            thresh   = (alpha);
            col_title = 'Raw p-value';
        end

        % Diagonal y=x reference
        plot(ax,[0 1],[0 1],'-','Color',S.lgrey,'LineWidth',0.8,'HandleVisibility','off');

        % Significance threshold lines
        xline(ax, thresh,'--','Color',[0.82 0.05 0.05],'LineWidth',0.9,'HandleVisibility','off');
        yline(ax, thresh,'--','Color',[0.82 0.05 0.05],'LineWidth',0.9,'HandleVisibility','off');

        % Quadrant shading
        % Bottom-left: both significant
        patch(ax,[0 thresh thresh 0],[0 0 thresh thresh], ...
            S.green,'FaceAlpha',0.08,'EdgeColor','none','HandleVisibility','off');
        % Bottom-right: t-test sig only
        patch(ax,[thresh 1 1 thresh],[0 0 thresh thresh], ...
            S.orange,'FaceAlpha',0.08,'EdgeColor','none','HandleVisibility','off');
        % Top-left: permutation sig only
        patch(ax,[0 thresh thresh 0],[thresh thresh 1 1], ...
            S.sky,'FaceAlpha',0.08,'EdgeColor','none','HandleVisibility','off');

        h_sz = gobjects(n_sizes,1);
            y_dataOK1   = (res.(mname).pval_chi2);
            y_dataOK2   = (res.(mname).pval); 
            
            for s = 1:n_sizes
            xv = x_data(:,s);
            yv1 = y_dataOK1(:,s);
            yv2 = y_dataOK2(:,s);
            yv  = y_data(:,s);
            ok = ~isnan(xv) & ~isnan(yv1)  & ~isnan(yv2);
            if ~any(ok), continue; end

            col   = size_cmap(s,:);
            h_sz(s) = scatter(ax, xv(ok), yv(ok), 16, col, 'filled', ...
                'MarkerFaceAlpha',0.45,'MarkerEdgeColor','k', ...
                'DisplayName',sprintf('n=%d',sample_sizes(s)));

            % Spearman rho annotation
            rho=nan;
            if sum(ok) >= 3 
                rho=1;
                if sum(xv(ok).*yv(ok))>0.01  %&&& try
                rho =corr(xv(ok), yv(ok), 'type','Pearson');
                end
          %      if (nnz(xv(ok).*yv(ok)==0)/nnz(xv(ok).*yv(ok)>=0))>0.9 , rho=1; end  %&&&
                
             %   rho=(0.01+nnz(xv(ok)<=alpha))/(0.01+nnz(yv(ok)<=alpha));

                text(ax, 0.97, 0.03 + (s-1)*0.09, ...
                    sprintf('n=%d: r=%.2f',sample_sizes(s),rho), ...
                    'Units','normalized','HorizontalAlignment','right', ...
                    'FontName',S.font,'FontSize',S.fsz_ax-0.5, ...
                    'Color',col*0.75,'Interpreter','tex');
            end
        end

        xlim(ax,[0 1]); ylim(ax,[0 1]);
        axis(ax,'square');
        xlabel(ax,xlab,'FontName',S.font,'FontSize',S.fsz);
        ylabel(ax,ylab,'FontName',S.font,'FontSize',S.fsz);

        if col_idx == 1
            title(ax, sprintf('%s', model_labels{m}), ...
                'FontName',S.font,'FontSize',S.fsz_title, ...
                'FontWeight','bold','Color',S.mcols(m,:));
        else
            title(ax, sprintf('%s ', model_labels{m}), ...
                'FontName',S.font,'FontSize',S.fsz_title, ...
                'FontWeight','bold','Color',S.mcols(m,:)*0.8);
        end

        if m == 1
            title(ax, sprintf('%s  — %s', model_labels{m}, col_title), ...
                'FontName',S.font,'FontSize',S.fsz_title, ...
                'FontWeight','bold','Color',S.mcols(m,:));
        end

        if strcmp(mname,'random'), ax.Color=[0.96 0.95 0.92]; end

        % Legend on first model 
  %      if m == 1 && col_idx == 2
  %          valid_h = arrayfun(@(h) isvalid(h) && isgraphics(h,'scatter'), h_sz);
  %          if any(valid_h)
  %              leg = legend(ax, h_sz(valid_h), ...
  %                  arrayfun(@(s) sprintf('n=%d',sample_sizes(s)), find(valid_h), 'UniformOutput',false), ...
  %                  'Location','eastoutside','FontName',S.font,'FontSize',S.fsz_ax,'Box','off');
  %              leg.Title.String = 'Sample size';
  %          end
  %      end

        panel_label(ax, panel_letters(pl), S);
        pub_ax(ax, S);
        hold(ax,'off');
    end
end

annotation(fig,'textbox',[0.02 0.001 0.96 0.025], ...
    'String',['Left: raw p-values  |  Right: BH-FDR q-values  |  ' ...
              'Each dot = one replicate  |  Colour = sample size  |  ' ...
              'Diagonal = perfect agreement between methods  |  ' ...
              'Green zone: both sig  |  Orange: t-test only  |  Blue: permutation only'], ...
    'FontName',S.font,'FontSize',6,'EdgeColor','none', ...
    'Color',S.grey,'HorizontalAlignment','center');

set(fig,'PaperUnits','centimeters','PaperSize',[28 32]);
end
% =========================================================================
% FIGURE Plot proportion agreeing on significant at different alpha values
% =========================================================================
function fig = plot_pval_comparisonA(res, S, model_names, model_labels, ...
        sample_sizes, n_sizes,pnotq)
% Bar/line plot showing number of significant edges per p-value type
% (Chi2, SR13 t-test, Permutation) across sample sizes, for each model.
% Rows    = models
% Columns = alpha levels
% Within each panel: grouped bars per sample size, one group per p-type.

alpha_levels = [0.001,0.003, 0.01, 0.05];
n_alpha      = numel(alpha_levels);
n_models     = numel(model_names);

% Colours for the three p-value types
col_chi2  = [0.20 0.50 0.80];
col_ttest = [0.15 0.65 0.35];
col_perm  = [0.75 0.22 0.18];
type_cols  = [col_chi2; col_ttest; col_perm];
type_names = {'Chi2','SR13 t-test','Permutation'};
n_types    = numel(type_names);

size_cmap  = parula(n_sizes);

fig = figure('Name','Fig5_SigCount_Compare','NumberTitle','off', ...
    'Color','w','Units','centimeters', ...
    'Position',[1 1 8*n_alpha 6*n_models]);

sgtitle('Number of significant edges by p-value type, sample size and \alpha', ...
    'FontName',S.font,'FontSize',S.fsz_title,'FontWeight','bold', ...
    'Color',S.black,'Interpreter','tex');

panel_letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
pl = 0;



for m = 1:n_models
    mname = model_names{m};

    for a = 1:n_alpha
        al  = alpha_levels(a);
        pl  = pl + 1;
        ax  = subplot(n_models, n_alpha, (m-1)*n_alpha + a);
        hold(ax, 'on');

        % For each sample size, count sig edges per p-value type
        % Matrix: n_sizes x n_types
        sig_counts = zeros(n_sizes, n_types);

        for s = 1:n_sizes
            % --- Chi2 ---
            pv = res.(mname).pval_chi2(:, s);
            ok = ~isnan(pv);
            if any(ok)
                [~,~,~,qv] = fdr_bh(pv(ok), al);
                sig_counts(s, 1) = sum(qv < al);
                if pnotq==1, sig_counts(s, 1) = sum(pv(ok) < al); end
            end

            % --- SR13 t-test ---
            pv = res.(mname).pval(:, s);
            ok = ~isnan(pv);
            if any(ok)
                [~,~,~,qv] = fdr_bh(pv(ok), al);
                sig_counts(s, 2) = sum(qv < al);
                if pnotq==1, sig_counts(s, 2) = sum(pv(ok) < al); end

            end

            % --- Permutation ---
            pv = res.(mname).pval_perm(:, s);
            ok = ~isnan(pv);
            if any(ok)
                [~,~,~,qv] = fdr_bh(pv(ok), al);
                sig_counts(s, 3) = sum(qv < al);
                if pnotq==1, sig_counts(s, 3) = sum(pv(ok) < al); end
        
            end
        end

        % --- Grouped bar chart ---
        % Each group = one sample size; bars within group = p-value type
        bw     = 0.22;
        offset = linspace(-(n_types-1)/2, (n_types-1)/2, n_types) * bw;
        x_base = 1:n_sizes;

        h_type = gobjects(n_types, 1);
        for tp = 1:n_types
            xpos = x_base + offset(tp);
            h_type(tp) = bar(ax, xpos, sig_counts(:, tp), bw, ...
                'FaceColor', type_cols(tp,:), ...
                'EdgeColor', 'none', ...
                'FaceAlpha', 0.82, ...
                'DisplayName', type_names{tp});
        end

        % Overlay lines connecting same-type counts across sample sizes
        for tp = 1:n_types
            xpos = x_base + offset(tp);
            plot(ax, xpos, sig_counts(:, tp), '-o', ...
                'Color', type_cols(tp,:)*0.7, ...
                'LineWidth', 1.0, 'MarkerSize', 4, ...
                'MarkerFaceColor', type_cols(tp,:), ...
                'HandleVisibility', 'off');
        end

        % Axes formatting
        ax.XTick      = x_base;
        ax.XTickLabel = arrayfun(@(n) sprintf('n=%d', n), sample_sizes, ...
                            'UniformOutput', false);
        ax.XTickLabelRotation = 45;
        ax.XTickMode          = 'manual';
        ax.XTickLabelMode     = 'manual';
        ax.FontSize           = S.fsz_ax;
        ax.TickLabelInterpreter = 'none';
        ax.XLim = [0.5, n_sizes + 0.5];
        ax.YLim = [0, max(sig_counts(:)) * 1.18 + 1];
        ax.GridAlpha = 0.13;
        grid(ax, 'on'); box(ax, 'on');

        ylabel(ax, 'Sig. edges (FDR)', 'FontName', S.font, 'FontSize', S.fsz);
        if pnotq==1,ylabel(ax, 'Sig. edges (P)', 'FontName', S.font, 'FontSize', S.fsz); end

        xlabel(ax, 'Sample size',       'FontName', S.font, 'FontSize', S.fsz);

        % Title: model on first alpha column only
        if a == 1
            title(ax, sprintf('%s  |  \\alpha=%.3f', model_labels{m}, al), ...
                'FontName', S.font, 'FontSize', S.fsz_title, ...
                'FontWeight', 'bold', 'Color', S.mcols(m,:), ...
                'Interpreter', 'tex');
        else
            title(ax, sprintf('\\alpha=%.3f', al), ...
                'FontName', S.font, 'FontSize', S.fsz_title, ...
                'FontWeight', 'bold', 'Color', S.black, ...
                'Interpreter', 'tex');
        end

        % Legend on first panel only
        if m == 1 && a == 1
            leg = legend(ax, h_type, type_names, ...
                'Location', 'northwest', 'FontName', S.font, ...
                'FontSize', S.fsz_ax, 'Box', 'off');
            leg.Title.String = 'p-value type';
        end

        if strcmp(mname, 'random'), ax.Color = [0.96 0.95 0.92]; end

        panel_label(ax, panel_letters(pl), S);
        pub_ax(ax, S);
        hold(ax, 'off');
    end
end

annotation(fig, 'textbox', [0.02 0.001 0.96 0.022], ...
    'String', ['Bar height = number of FDR-corrected significant edges  |  ' ...
               'Columns = alpha levels  |  Rows = models  |  ' ...
               'Blue = Chi2  |  Green = SR13 t-test  |  Red = Permutation'], ...
    'FontName', S.font, 'FontSize', 6, 'EdgeColor', 'none', ...
    'Color', S.grey, 'HorizontalAlignment', 'center');
 if pnotq==1
     annotation(fig, 'textbox', [0.02 0.001 0.96 0.022], ...
    'String', ['Bar height = number of significant edges  |  ' ...
               'Columns = alpha levels  |  Rows = models  |  ' ...
               'Blue = Chi2  |  Green = SR13 t-test  |  Red = Permutation'], ...
    'FontName', S.font, 'FontSize', 6, 'EdgeColor', 'none', ...
    'Color', S.grey, 'HorizontalAlignment', 'center');
 
 end
set(fig, 'PaperUnits', 'centimeters', 'PaperSize', [8*n_alpha, 6*n_models]);
end

% =========================================================================
%  FIGURE  — CI COVERAGE OF BIAS


% =========================================================================
%  FIGURE 6 — CI COVERAGE OF BIAS
%
%  For each model and sample size, shows:
%    - The distribution of subsample dCor values (box)
%    - The three CI bands (Fisher z, Bayes BB, Bernstein) computed from
%      those replicates
%    - The TRUE dCor (full N=1000) as a reference line
%    - The BIAS as a filled segment between mean(subsample) and true dCor
%    - A coverage indicator: green if the CI contains the true value,
%      red if it does not
%
%  Bottom row per model: empirical coverage rate (fraction of CIs that
%  contain the true dCor) vs sample size for each CI method.
% =========================================================================
function fig = plot_ci_coverage(res, S, model_names, model_labels, ...
        sample_sizes, n_sizes)

n_models    = numel(model_names);
alpha       = S.params.alpha;
ci_pct      = 100*(1-alpha);
x_pos       = 1:n_sizes;
x_tick      = arrayfun(@num2str, sample_sizes, 'UniformOutput',false);
panel_letters = 'ABCDEFGHIJ';
pl = 0;

fig = figure('Name','Fig6_dCor_CI_Coverage','NumberTitle','off', ...
    'Color','w','Units','centimeters','Position',[1 1 24 32]);

sgtitle(sprintf(['CI Coverage of True dCor  |  %.0f%% CI  |  ' ...
    'Does each CI method bracket the true value?'], ci_pct), ...
    'FontName',S.font,'FontSize',S.fsz_title,'FontWeight','bold','Color',S.black);

for m = 1:n_models
    mname = model_names{m};
    ref   = res.(mname).ref_dcor;

    ci_fz   = res.(mname).ci_fz;     % [n_sizes x 2]
    ci_bb   = res.(mname).ci_bb;
    ci_bern = res.(mname).ci_bern;

    % Per-CI coverage: fraction of individual replicate CIs that contain ref.
    % For each replicate r and sample size s, compute a per-replicate CI
    % (using the single observed dCor value) and check containment.
    % We use the stored distribution-level CIs for the band plot, but
    % compute per-replicate containment for the coverage rate.
    cov_fz   = NaN(1, n_sizes);
    cov_bb   = NaN(1, n_sizes);
    cov_bern = NaN(1, n_sizes);

    for s = 1:n_sizes
        d  = res.(mname).dcor_u(:,s);
        ok = ~isnan(d);
        if sum(ok) < 2, continue; end
        d_ok = d(ok);
        ns   = sample_sizes(s);
        n_ok = sum(ok);

        % Per-replicate CIs — Fisher z on each individual value
        lo_fz = NaN(n_ok,1); hi_fz = NaN(n_ok,1);
        lo_bb_r = NaN(n_ok,1); hi_bb_r = NaN(n_ok,1);
        lo_bn = NaN(n_ok,1); hi_bn = NaN(n_ok,1);

        for r = 1:n_ok
            r_val = d_ok(r);
            % Fisher z on single value
            if ns > 3 && ~isnan(r_val)
                rv = max(-0.9999,min(0.9999,r_val));
                z  = atanh(rv); se = 1/sqrt(max(ns-3,1));
                zc = norminv(1-alpha/2);
                lo_fz(r) = tanh(z-zc*se); hi_fz(r) = tanh(z+zc*se);
            end
            % Bernstein on single value (s2=0 → pure radius)
            if ns >= 2 && ~isnan(r_val)
                lnt = log(2/alpha);
                eps = (2/3)*lnt/ns;
                lo_bn(r) = max(-1,r_val-eps); hi_bn(r) = min(1,r_val+eps);
            end
        end
        % Bayes BB: on the whole replicate distribution
        ci_bb_s = res.(mname).ci_bb(s,:);
        lo_bb_r(:) = ci_bb_s(1); hi_bb_r(:) = ci_bb_s(2);

        cov_fz(s)   = mean(lo_fz   <= ref & hi_fz   >= ref, 'omitnan');
        cov_bb(s)   = mean(lo_bb_r <= ref & hi_bb_r >= ref, 'omitnan');
        cov_bern(s) = mean(lo_bn   <= ref & hi_bn   >= ref, 'omitnan');
    end

    % ---- Top panel: distribution + CI bands + bias ----
    pl  = pl + 1;
    ax1 = subplot(n_models, 2, (m-1)*2 + 1);
    hold(ax1,'on');

    % CI bands
    draw_ci_band(ax1,x_pos,ci_fz,  S.fz_fill,  S.fz_col,  1.1,0.22,'-');
    draw_ci_band(ax1,x_pos,ci_bb,  S.bb_fill,  S.bb_col,  1.3,0.26,'--');
    draw_ci_band(ax1,x_pos,ci_bern,S.bern_fill,S.bern_col,1.3,0.30,'-.');

    % Bias segment: mean dCor → true dCor per sample size
    for s = 1:n_sizes
        d   = res.(mname).dcor_u(:,s);
        mn  = nanmean(d);
        if isnan(mn), continue; end
        % Bias bar: thick line from mean to ref
        bias_col = S.orange;
        if mn > ref
            bias_col = S.red;    % overestimation
        elseif mn < ref
            bias_col = S.blue;   % underestimation
        end
        plot(ax1,[x_pos(s),x_pos(s)],[mn,ref],'-', ...
            'Color',bias_col,'LineWidth',3.0,'HandleVisibility','off');
        % Mark the mean
        plot(ax1,x_pos(s),mn,'o','MarkerSize',4, ...
            'MarkerFaceColor',bias_col,'MarkerEdgeColor','w','LineWidth',0.5, ...
            'HandleVisibility','off');
    end

    % IQR box for the replicate distribution
    for s = 1:n_sizes
        d = res.(mname).dcor_u(:,s);
        d = d(~isnan(d));
        if numel(d)<4, continue; end
        q1 = prctile(d,25); q3 = prctile(d,75);
        patch(ax1,[x_pos(s)-0.15,x_pos(s)+0.15,x_pos(s)+0.15,x_pos(s)-0.15], ...
            [q1,q1,q3,q3],S.mcols(m,:), ...
            'FaceAlpha',0.15,'EdgeColor',S.mcols(m,:),'LineWidth',0.8, ...
            'HandleVisibility','off');
    end

    % True dCor reference
    yline(ax1,ref,'--','Color',[0.80 0.05 0.05],'LineWidth',1.4, ...
        'Label',sprintf('True  R=%.3f',ref), ...
        'LabelHorizontalAlignment','right', ...
        'FontName',S.font,'FontSize',S.fsz_ax,'HandleVisibility','off');
    yline(ax1,0,':','Color',S.grey,'LineWidth',0.6,'HandleVisibility','off');

    ylim(ax1,[-0.05 1.05]);
    xlim(ax1,[0.4, n_sizes+0.6]);
    set(ax1,'XTick',x_pos,'XTickLabel',x_tick);
    xlabel(ax1,'Sample size  n','FontName',S.font,'FontSize',S.fsz);
    ylabel(ax1,'dCor (U-stat)','FontName',S.font,'FontSize',S.fsz);
    title(ax1, sprintf('%s  — CI bands & bias', model_labels{m}), ...
        'FontName',S.font,'FontSize',S.fsz_title,'FontWeight','bold', ...
        'Color',S.mcols(m,:));

    if strcmp(mname,'random'), ax1.Color=[0.96 0.95 0.92]; end
    panel_label(ax1, panel_letters(pl), S);
    pub_ax(ax1, S);
    hold(ax1,'off');

    % ---- Right panel: coverage rate vs sample size ----
    pl  = pl + 1;
    ax2 = subplot(n_models, 2, (m-1)*2 + 2);
    hold(ax2,'on');

    % Nominal level
    yline(ax2,1-alpha,'--','Color',[0.82 0.05 0.05],'LineWidth',1.2, ...
        'Label',sprintf('nominal %.0f%%',100*(1-alpha)), ...
        'LabelHorizontalAlignment','right', ...
        'FontName',S.font,'FontSize',S.fsz_ax,'HandleVisibility','off');

    valid = ~isnan(cov_fz);
    if any(valid)
        plot(ax2,sample_sizes(valid),cov_fz(valid)*100,  '-o', ...
            'Color',S.fz_col,  'LineWidth',S.lw,'MarkerSize',S.ms, ...
            'MarkerFaceColor',S.fz_col,  'DisplayName','Fisher z');
        plot(ax2,sample_sizes(valid),cov_bb(valid)*100,  '--s', ...
            'Color',S.bb_col,  'LineWidth',S.lw,'MarkerSize',S.ms, ...
            'MarkerFaceColor',S.bb_col,  'DisplayName','Bayes BB');
        plot(ax2,sample_sizes(valid),cov_bern(valid)*100,'-.^', ...
            'Color',S.bern_col,'LineWidth',S.lw,'MarkerSize',S.ms, ...
            'MarkerFaceColor',S.bern_col,'DisplayName','Bernstein');
    end

    % Shade over-coverage zone
    patch(ax2,[sample_sizes(1)*0.8 sample_sizes(end)*1.3 ...
               sample_sizes(end)*1.3 sample_sizes(1)*0.8], ...
        [100*(1-alpha) 100*(1-alpha) 101 101], ...
        S.green,'FaceAlpha',0.06,'EdgeColor','none','HandleVisibility','off');

    set(ax2,'XScale','log','XTick',sample_sizes);
    xlim(ax2,[sample_sizes(1)*0.8, sample_sizes(end)*1.3]);
    ylim(ax2,[0 101]);
    xlabel(ax2,'Sample size  n','FontName',S.font,'FontSize',S.fsz);
    ylabel(ax2,'Coverage (%)','FontName',S.font,'FontSize',S.fsz);
    title(ax2, sprintf('%s  — coverage rate', model_labels{m}), ...
        'FontName',S.font,'FontSize',S.fsz_title,'FontWeight','bold', ...
        'Color',S.mcols(m,:));

    if m==1
        leg = legend(ax2,'Location','southeast','FontName',S.font, ...
            'FontSize',S.fsz_ax,'Box','off');
    end

    if strcmp(mname,'random'), ax2.Color=[0.96 0.95 0.92]; end
    panel_label(ax2, panel_letters(pl), S);
    pub_ax(ax2, S);
    hold(ax2,'off');
end

annotation(fig,'textbox',[0.02 0.001 0.96 0.025], ...
    'String',['Left: CI bands (3 methods) + bias bar (mean subsample dCor to true dCor)  |  ' ...
              'Red bias = overestimate, Blue = underestimate  |  Box = IQR of replicates  |  ' ...
              'Right: fraction of per-replicate CIs that contain the true dCor  |  ' ...
              'Green zone = above nominal coverage (conservative)'], ...
    'FontName',S.font,'FontSize',6,'EdgeColor','none', ...
    'Color',S.grey,'HorizontalAlignment','center');

set(fig,'PaperUnits','centimeters','PaperSize',[24 32]);
end



% =========================================================================
%  FIGURE  — STATISTICAL SIGNIFICANCE OF dCor DIFFERENCE
%             (subsample vs full N)
%
%  For each model and sample size:
%    The question: is the subsample dCor significantly DIFFERENT from the
%    full-N dCor?  This is NOT about whether the correlation itself is
%    significant, but whether the estimation error due to small n is
%    statistically detectable.
%
%  Method: bootstrap hypothesis test.
%    H0: mean(subsample dCor) = true dCor (full N=1000)
%    Under H0, the distribution of (subsample mean - true) centred at 0.
%    We compute a two-sided p-value as:
%      p = 2 * min(fraction > 0, fraction < 0)
%    where the fraction is over the n_reps replicates.
%    This is equivalent to a sign test on the direction of bias.
%
%  Additionally shows:
%    - Effect size: Cohen's d = bias / SD(subsample dCor)
%    - BH-FDR correction across all (model x sample_size) cells
%
%  Layout:
%    Fig 7a: heatmap — -log10(q) for each (model x sample_size) cell
%    Fig 7b: effect size (bias/SD) heatmap
%    Fig 7c: line plots — bias and SD vs log(n) per model (same as bias
%            plot but now annotated with significance stars)
% =========================================================================
function fig = plot_dcor_diff_significance(res, S, model_names, model_labels, ...
        sample_sizes, n_sizes)

n_models = numel(model_names);
alpha    = S.params.alpha;
n_reps   = S.params.n_reps;
x_tick   = arrayfun(@num2str, sample_sizes, 'UniformOutput',false);
x_pos    = 1:n_sizes;

% ------------------------------------------------------------------
% Compute per-cell statistics
%   bias(m,s)    = mean(subsample) - true
%   sd_mat(m,s)  = std(subsample)
%   pval(m,s)    = two-sided sign test p-value
%   effsize(m,s) = bias / sd  (Cohen's d analog) change to Cohen's dav, 
% ------------------------------------------------------------------
bias_mat = NaN(n_models, n_sizes);
sd_mat   = NaN(n_models, n_sizes);
pval_mat = NaN(n_models, n_sizes);

for m = 1:n_models
    mname = model_names{m};
    ref   = res.(mname).ref_dcor;
    for s = 1:n_sizes
        d  = res.(mname).dcor_u(:,s);
        ok = ~isnan(d);
        if sum(ok) < 5, continue; end
        d_ok  = d(ok);
        b     = mean(d_ok) - ref;        % bias
        sd    = std(d_ok);
        % Two-sided sign test: fraction above/below ref
        frac_above = mean(d_ok > ref);
        frac_below = mean(d_ok < ref);
        p = 2 * min(frac_above, frac_below);
        p = max(p, 1/n_reps);            % floor at 1/n_reps
        bias_mat(m,s) = b;
        sd_mat(m,s)   = sd;
        pval_mat(m,s) = p;
    end
end

% Effect size: bias / SD (Cohen's d analog)
effsize_mat = abs(bias_mat) ./ (sd_mat+sd_mat(:,end))/2;

% BH-FDR across all cells jointly
p_flat      = pval_mat(:);
        [~, ~, ~, q_flat]=fdr_bh(p_flat,alpha);


% [~, q_flat]=mafdr(p_flat); %MCC &&&
qval_mat    = reshape(q_flat, n_models, n_sizes);


%qval_mat    = reshape(p_flat, n_models, n_sizes); %MCC try this &&&

logq_mat    = -log10(max(qval_mat, 1e-30));

sig_mat     = qval_mat < alpha;   % significant difference from full-N

% ------------------------------------------------------------------
% Figure
% ------------------------------------------------------------------
fig = figure('Name','Fig7_dCor_Diff_Signif','NumberTitle','off', ...
    'Color','w','Units','centimeters','Position',[1 1 24 22]);

sgtitle(sprintf(['Statistical Significance of dCor Difference:  ' ...
    'Subsample vs Full N=1000  |  BH-FDR \alpha=%.2f'], alpha), ...
    'FontName',S.font,'FontSize',S.fsz_title,'FontWeight','bold', ...
    'Color',S.black,'Interpreter','tex');

% ---- Panel A: -log10(q) heatmap ----
ax1 = subplot(2,3,[1 2]);
imagesc(ax1, logq_mat);
colormap(ax1, hot(256));
q_max = max(logq_mat(:));
clim(ax1,[0, max(q_max,0.1)]);
cb1 = colorbar(ax1);
cb1.Label.String = '-log_{10}(q)  [BH-FDR]';
cb1.Label.FontSize = S.fsz_ax;
cb1.FontSize = S.fsz_ax;
cb1.LineWidth = S.lw_ax;
% Mark significance threshold line on colourbar
q_thr_line = -log10(alpha);
if q_thr_line <= q_max
    yline_cb = q_thr_line / max(q_max, 0.1);
    annotation(fig,'line',...
        ax1.Position(1)+ax1.Position(3)+0.005 + [0.02 0.04], ...
        ax1.Position(2) + ax1.Position(4)*yline_cb * [1 1], ...
        'Color',[0.3 0.3 0.3],'LineWidth',0.8);
end
set(ax1,'YTick',1:n_models,'YTickLabel',model_labels, ...
    'XTick',x_pos,'XTickLabel',x_tick,'TickLength',[0 0],'FontSize',S.fsz_ax);
xlabel(ax1,'Sample size  n','FontName',S.font,'FontSize',S.fsz);
title(ax1,'Significance: -log_{10}(q)  [larger = more significant difference]', ...
    'FontName',S.font,'FontSize',S.fsz_title,'FontWeight','bold','Interpreter','tex');
% Overlay significance stars and value text
for m = 1:n_models
    for s = 1:n_sizes
        if isnan(qval_mat(m,s)), continue; end
        if sig_mat(m,s)
            txt = sprintf('q=%.3f*', qval_mat(m,s));
            tcol = [1 1 1];
        else
            txt = sprintf('q=%.3f', qval_mat(m,s));
            tcol = [0.6 0.6 0.6];
        end
        text(ax1, s, m, txt, 'HorizontalAlignment','center', ...
            'VerticalAlignment','middle', ...
            'FontName',S.font,'FontSize',max(4,S.fsz_ax-1),'Color',tcol);
    end
end
pub_ax(ax1,S); set(ax1,'Box','on','TickDir','in');
panel_label(ax1,'A',S);

% ---- Panel B: Effect size heatmap ----
ax2 = subplot(2,3,3);
% Diverging colourmap: blue=underestimate, white=no bias, red=overestimate
n_c   = 256;
cmap_div = [linspace(S.blue(1),1,n_c/2)', linspace(S.blue(2),1,n_c/2)', linspace(S.blue(3),1,n_c/2)'; ...
            linspace(1,S.red(1),n_c/2)',  linspace(1,S.red(2),n_c/2)',  linspace(1,S.red(3),n_c/2)'];
imagesc(ax2, effsize_mat);
colormap(ax2, cmap_div);
e_max = max(abs(effsize_mat(:)),[],'omitnan');
if isnan(e_max) || e_max == 0, e_max = 1; end
clim(ax2,[-e_max e_max]);
cb2 = colorbar(ax2);
cb2.Label.String = 'Effect size  (bias / SD)';
cb2.Label.FontSize = S.fsz_ax;
cb2.FontSize = S.fsz_ax;  cb2.LineWidth = S.lw_ax;
set(ax2,'YTick',1:n_models,'YTickLabel',model_labels, ...
    'XTick',x_pos,'XTickLabel',x_tick,'TickLength',[0 0],'FontSize',S.fsz_ax);
xlabel(ax2,'Sample size  n','FontName',S.font,'FontSize',S.fsz);
title(ax2,'Effect size  (bias / SD)','FontName',S.font,'FontSize',S.fsz_title,'FontWeight','bold');
for m = 1:n_models
    for s = 1:n_sizes
        if isnan(effsize_mat(m,s)), continue; end
        txt = sprintf('%.2f', effsize_mat(m,s));
        tcol = [0.15 0.15 0.15];
        if abs(effsize_mat(m,s)) > 0.6*e_max, tcol=[1 1 1]; end
        text(ax2,s,m,txt,'HorizontalAlignment','center','VerticalAlignment','middle', ...
            'FontName',S.font,'FontSize',max(4,S.fsz_ax-1),'Color',tcol);
    end
end
pub_ax(ax2,S); set(ax2,'Box','on','TickDir','in');
panel_label(ax2,'B',S);

% ---- Panel C: Bias lines with significance stars ----
ax3 = subplot(2,3,[4 5]);
hold(ax3,'on');
yline(ax3,0,'-','Color',S.grey,'LineWidth',0.7,'HandleVisibility','off');

for m = 1:n_models
    col   = S.mcols(m,:);
    valid = ~isnan(bias_mat(m,:));
    if ~any(valid), continue; end
    plot(ax3, sample_sizes(valid), bias_mat(m,valid), '-o', ...
        'Color',col,'LineWidth',S.lw,'MarkerSize',S.ms, ...
        'MarkerFaceColor',col,'DisplayName',model_labels{m});
    % Significance stars above/below each point
    for s = find(valid)
        if sig_mat(m,s)
            offset = sign(bias_mat(m,s)) * 0.012;
            text(ax3, sample_sizes(s), bias_mat(m,s)+offset, '*', ...
                'HorizontalAlignment','center','FontSize',9, ...
                'Color',col,'FontWeight','bold');
        end
    end
end

set(ax3,'XScale','log','XTick',sample_sizes,'XTickLabel',x_tick);
xlim(ax3,[sample_sizes(1)*0.8, sample_sizes(end)*1.3]);
xlabel(ax3,'Sample size  n','FontName',S.font,'FontSize',S.fsz);
ylabel(ax3,'Bias  (R_{sub} - R_{true})','FontName',S.font,'FontSize',S.fsz,'Interpreter','tex');
title(ax3,'Bias vs sample size  (* = FDR-significant difference from full N)', ...
    'FontName',S.font,'FontSize',S.fsz_title,'FontWeight','bold');
legend(ax3,'Location','northeast','FontName',S.font,'FontSize',S.fsz_ax,'Box','off');
pub_ax(ax3,S);
panel_label(ax3,'C',S);
hold(ax3,'off');

% ---- Panel D: SD with significance ----
ax4 = subplot(2,3,6);
hold(ax4,'on');

for m = 1:n_models
    col   = S.mcols(m,:);
    valid = ~isnan(sd_mat(m,:));
    if ~any(valid), continue; end
    plot(ax4, sample_sizes(valid), sd_mat(m,valid), '-o', ...
        'Color',col,'LineWidth',S.lw,'MarkerSize',S.ms, ...
        'MarkerFaceColor',col,'DisplayName',model_labels{m});
    % Mark which sample sizes show significant difference
    for s = find(valid)
        if sig_mat(m,s)
            plot(ax4,sample_sizes(s),sd_mat(m,s),'o', ...
                'MarkerSize',S.ms+3,'MarkerFaceColor','none', ...
                'MarkerEdgeColor',col,'LineWidth',1.2,'HandleVisibility','off');
        end
    end
end

% 1/sqrt(n) reference
all_sd = sd_mat(~isnan(sd_mat));
if ~isempty(all_sd)
    sc     = median(all_sd) * sqrt(sample_sizes(ceil(n_sizes/2)));
    n_ref  = linspace(sample_sizes(1), sample_sizes(end), 100);
    plot(ax4, n_ref, sc./sqrt(n_ref), ':', 'Color',S.grey, 'LineWidth',0.8, ...
        'DisplayName','1/\surd n ref');
end

set(ax4,'XScale','log','XTick',sample_sizes,'XTickLabel',x_tick);
xlim(ax4,[sample_sizes(1)*0.8, sample_sizes(end)*1.3]);
ylim(ax4,[0 Inf]);
xlabel(ax4,'Sample size  n','FontName',S.font,'FontSize',S.fsz);
ylabel(ax4,'SD(R_{sub})','FontName',S.font,'FontSize',S.fsz,'Interpreter','tex');
title(ax4,'SD  (open circle = sig. diff. from full N)', ...
    'FontName',S.font,'FontSize',S.fsz_title,'FontWeight','bold');
pub_ax(ax4,S);
panel_label(ax4,'D',S);
hold(ax4,'off');

annotation(fig,'textbox',[0.02 0.001 0.96 0.025], ...
    'String',['H0: subsample dCor = true dCor (full N=1000)  |  ' ...
              'Two-sided sign test across n_reps replicates + BH-FDR across all (model x n) cells  |  ' ...
              'Panel A: -log10(q): bright=significant difference from true  |  ' ...
              'Panel B: effect size (bias/SD): blue=underestimate, red=overestimate  |  ' ...
              'Panel C/D: * and open circles mark FDR-significant cells'], ...
    'FontName',S.font,'FontSize',6,'EdgeColor','none', ...
    'Color',S.grey,'HorizontalAlignment','center');

set(fig,'PaperUnits','centimeters','PaperSize',[24 22]);
end


% =========================================================================
%  HELPERS
% =========================================================================

function pub_ax(ax, S)
set(ax, 'FontName', S.font, 'FontSize', S.fsz_ax, ...
    'Box', 'off', 'TickDir', 'out', 'TickLength', [0.012 0.012], ...
    'LineWidth', S.lw_ax, 'XColor', S.black, 'YColor', S.black, ...
    'GridAlpha', 0.10, 'GridColor', [0.6 0.6 0.6], 'GridLineStyle', ':');
grid(ax, 'on');
end

function panel_label(ax, letter, S)
text(ax, -0.14, 1.08, letter, 'Units', 'normalized', ...
    'FontName', S.font, 'FontSize', S.fsz_title + 2, ...
    'FontWeight', 'bold', 'Color', S.black, ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
end

function draw_ci_band(ax, x_pos, bands, fill_col, line_col, lw, tick_w, ls)
lo = bands(:,1)';  hi = bands(:,2)';
valid = ~isnan(lo) & ~isnan(hi);
n_s   = numel(x_pos);
if any(valid)
    xv = [x_pos(valid), fliplr(x_pos(valid))];
    yv = [lo(valid),    fliplr(hi(valid))];
    patch(ax, xv, yv, fill_col, 'FaceAlpha', 0.30, 'EdgeColor', 'none');
end
if sum(valid) >= 2
    xs = x_pos(valid);
    xi = linspace(xs(1), xs(end), 120);
    lo_f = interp1(xs, lo(valid), xi, 'pchip');
    hi_f = interp1(xs, hi(valid), xi, 'pchip');
    plot(ax, xi, lo_f, ls, 'Color', line_col, 'LineWidth', lw);
    plot(ax, xi, hi_f, ls, 'Color', line_col, 'LineWidth', lw);
end
for s = 1:n_s
    if isnan(bands(s,1)), continue; end
    plot(ax, [x_pos(s)-tick_w, x_pos(s)+tick_w], [bands(s,1), bands(s,1)], ...
        ls, 'Color', line_col, 'LineWidth', lw);
    plot(ax, [x_pos(s)-tick_w, x_pos(s)+tick_w], [bands(s,2), bands(s,2)], ...
        ls, 'Color', line_col, 'LineWidth', lw);
end
end

function add_qvalue_colourbar(fig, S, alpha)
ax_cb = axes(fig, 'Position', [0.925 0.08 0.018 0.55]);
n_g   = 128;
grad  = linspace(0, 1, n_g)';
cdata = grad * S.col_sig + (1 - grad) * S.col_nonsig;
image(ax_cb, 1, linspace(0,1,n_g), permute(cdata,[1 3 2]));
q_thr = -log10(alpha);
q_ticks  = [alpha, alpha/5, alpha/50];
yt_norm  = min(-log10(q_ticks) / (q_thr*2), 1);
set(ax_cb, 'XTick', [], ...
    'YTick',      yt_norm, ...
    'YTickLabel', arrayfun(@(q) sprintf('q=%.3f',q), q_ticks, 'UniformOutput',false), ...
    'YDir', 'normal', 'FontName', S.font, 'FontSize', 6, ...
    'YAxisLocation', 'right', 'LineWidth', S.lw_ax);
yline(ax_cb, yt_norm(1), '-', 'Color',[0.3 0.3 0.3], 'LineWidth', 0.8);
ylabel(ax_cb, '-log_{10}(q)  [BH-FDR]', 'FontName', S.font, 'FontSize', 6.5, ...
    'Rotation', 270, 'VerticalAlignment', 'bottom', 'Interpreter', 'tex');
title(ax_cb, sprintf('q\x2264%.2f', alpha), 'FontName', S.font, 'FontSize', 6.5);
end

function save_fig(fig, name, S)
if ~S.export, return; end
% Create output directory if it does not exist
if ~isfolder(S.outdir)
    mkdir(S.outdir);
    fprintf('  Created output directory: %s\n', S.outdir);
end
fpath = fullfile(S.outdir, name);
switch lower(S.fmt)
    case 'pdf'
        set(fig, 'PaperPositionMode', 'auto');
        print(fig, fpath, '-dpdf', '-bestfit');
    case 'png'
        print(fig, fpath, '-dpng', '-r300');
    case 'svg'
        print(fig, fpath, '-dsvg');
end
fprintf('  Saved: %s.%s\n', fpath, S.fmt);
end

function results = mc_functional_dependence(varargin)
% MCuperlovicCulf, Ottawa 2026
% MC_FUNCTIONAL_DEPENDENCE
%   Monte Carlo simulation of 1000 sample points for two features X and Y
%   under five types of functional dependence with additive user-defined noise.
%
%   Dependence types:
%     Linear    : Y = a*X + b                  + ε
%     Quadratic : Y = a*X^2 + b*X + c          + ε
%     Cubic     : Y = a*X^3 + b*X^2 + c*X + d  + ε
%     sin  : Y = sin(A + B*X)         + ε
%     Random    : X ~ U(x_range), Y ~ U(x_range) independently (no dependence)
%
% USAGE:
%   results = mc_functional_dependence()
%   results = mc_functional_dependence('n', 1000, 'noise_std', 0.2, 'plot', true)
%   results = mc_functional_dependence('params', params_struct)
%
% OPTIONAL NAME-VALUE PAIRS:
%   'n'           - number of Monte Carlo samples (default: 1000)
%   'x_range'     - [xmin, xmax] for X domain (default: [-3, 3])
%   'x_dist'      - 'uniform' or 'normal' for X distribution (default: 'uniform')
%   'noise_std'   - scalar or 5-element vector of noise std devs per model
%                   Order: [linear, quad, cubic, gauss, random]
%                   (default: 0.2; random entry is ignored — Y is pure noise)
%   'noise_dist'  - 'normal', 'uniform', 'laplace', or 'student' (default: 'normal')
%   'noise_df'    - degrees of freedom if noise_dist='student' (default: 5)
%   'params'      - struct with fields: linear, quadratic, cubic, sin
%   'seed'        - random seed for reproducibility (default: 42)
%   'plot'        - true/false (default: true)
%   'compute_dcor'- true/false (default: true)
%
% OUTPUTS:
%   results  - struct with one field per model: linear/quadratic/cubic/sin/random
%       .X            (n x 1) input samples
%       .Y            (n x 1) noisy output
%       .Y_clean      (n x 1) noise-free output (= Y for random)
%       .noise        (n x 1) noise realization  (zeros for random)
%       .params       struct of parameters used
%       .noise_std    scalar noise std used
%       .dcor         distance correlation (U-statistic)
%       .pearson_r    Pearson correlation
%       .spearman_r   Spearman correlation
%
% REFERENCES:
%   [1] Székely, Rizzo & Bakirov (2007). Ann. Statist., 35(6), 2769-2794.
%   [2] Székely & Rizzo (2014). Ann. Statist., 42(6), 2382-2412.

% -------------------------------------------------------------------------
% Parse inputs
% -------------------------------------------------------------------------
p = inputParser;
addParameter(p, 'n',            1000);
addParameter(p, 'x_range',      [-3, 3]);
addParameter(p, 'x_dist',       'uniform', @(x) ismember(x,{'uniform','normal'}));
addParameter(p, 'noise_std',    0.2);
addParameter(p, 'noise_dist',   'normal',  @(x) ismember(x,{'normal','uniform','laplace','student'}));
addParameter(p, 'noise_df',     5);
addParameter(p, 'params',       []);
addParameter(p, 'seed',         42);
addParameter(p, 'plot',         true);
addParameter(p, 'compute_dcor', true);
parse(p, varargin{:});

n           = p.Results.n;
x_range     = p.Results.x_range;
x_dist      = p.Results.x_dist;
noise_std   = p.Results.noise_std;
noise_dist  = p.Results.noise_dist;
noise_df    = p.Results.noise_df;
user_params = p.Results.params;
seed        = p.Results.seed;
do_plot     = p.Results.plot;
do_dcor     = p.Results.compute_dcor;

rng(seed);

% Expand noise_std to 5-vector
if isscalar(noise_std)
    noise_std = repmat(noise_std, 1, 5);
elseif numel(noise_std) == 4
    noise_std = [noise_std, 0];
end
assert(numel(noise_std) == 5, 'noise_std must be scalar, 4- or 5-element vector.');

% -------------------------------------------------------------------------
% Default model parameters
% Chosen so that signal_std ~ 0.5-1.4, giving correlation spread from
% r~0.9 at low noise (sigma=0.1) down to r~0.05 at high noise (sigma=16)
% across all four functional models.  See noise_correlation_analysis.m.
% -------------------------------------------------------------------------
def.linear    = struct('a', 0.8,     'b', 0);
def.quadratic = struct('a', 0.20,    'b', 0,   'c', 0);
def.cubic     = struct('a', 0.05,    'b', 0,   'c', 0.3, 'd', 0);
def.sin  = struct('A', 0.2, 'B', 2);
def.random    = struct();

if ~isempty(user_params)
    fnames = fieldnames(user_params);
    for i = 1:numel(fnames)
        def.(fnames{i}) = user_params.(fnames{i});
    end
end

% -------------------------------------------------------------------------
% Generate shared X (used by all functional models)
% -------------------------------------------------------------------------
switch x_dist
    case 'uniform'
        X = x_range(1) + (x_range(2) - x_range(1)) * rand(n, 1);
    case 'normal'
        mu_x  = mean(x_range);
        sig_x = (x_range(2) - x_range(1)) / 6;
        X     = mu_x + sig_x * randn(n, 1);
        X     = max(x_range(1), min(x_range(2), X));
end

% -------------------------------------------------------------------------
% Model table  [label | fn | params | noise_sigma]
% -------------------------------------------------------------------------
models = {
    'Linear',    @(x,pr) pr.a*x + pr.b,                             def.linear,    noise_std(1);
    'Quadratic', @(x,pr) pr.a*x.^2 + pr.b*x + pr.c,                def.quadratic, noise_std(2);
    'Cubic',     @(x,pr) pr.a*x.^3 + pr.b*x.^2 + pr.c*x + pr.d,   def.cubic,     noise_std(3);
    'sin',  @(x,pr) sin((pr.A + pr.B*x)),                  def.sin,  noise_std(4);
    'Random',    [],                                                  def.random,    0;
};

model_labels = models(:,1);
model_fns    = models(:,2);
model_params = models(:,3);
model_noise  = cell2mat(models(:,4));
field_names  = {'linear','quadratic','cubic','sin','random'};
n_models     = 5;

% -------------------------------------------------------------------------
% Simulate
% -------------------------------------------------------------------------
results = struct();

fprintf('\n========================================================\n');
fprintf('  Monte Carlo Functional Dependence Simulation\n');
fprintf('  n = %d samples | noise: %s(0, sigma)\n', n, noise_dist);
fprintf('========================================================\n\n');
fprintf('%-12s  %-8s  %-10s  %-10s  %-10s\n', ...
    'Model', 'Noise s', 'Pearson', 'Spearman', 'dCor(U)');
fprintf('%s\n', repmat('-', 1, 58));

for m = 1:n_models
    fname = field_names{m};
    pr    = model_params{m};
    sig   = model_noise(m);

    if strcmp(fname, 'random')
        % Independent pair: both drawn fresh from uniform
        X_m     = x_range(1) + (x_range(2)-x_range(1)) * rand(n,1);
        Y_clean = x_range(1) + (x_range(2)-x_range(1)) * rand(n,1);
        eps     = zeros(n,1);
        Y       = Y_clean;
    else
        X_m     = X;
        fn      = model_fns{m};
        Y_clean = fn(X_m, pr);
        eps     = generate_noise(n, sig, noise_dist, noise_df);
        Y       = Y_clean + eps;
    end

    pearson_r  = corr(X_m, Y, 'type','Pearson');
    spearman_r = corr(X_m, Y, 'type','Spearman');
    dc = NaN;
    if do_dcor
        dc = dcor_U(squareform(pdist(X_m)), squareform(pdist(Y)), n);
    end

    results.(fname).X          = X_m;
    results.(fname).Y          = Y;
    results.(fname).Y_clean    = Y_clean;
    results.(fname).noise      = eps;
    results.(fname).params     = pr;
    results.(fname).noise_std  = sig;
    results.(fname).dcor       = dc;
    results.(fname).pearson_r  = pearson_r;
    results.(fname).spearman_r = spearman_r;

    fprintf('%-12s  %-8.3f  %-10.4f  %-10.4f  %-10.4f\n', ...
        model_labels{m}, sig, pearson_r, spearman_r, dc);
end
fprintf('\n');

if do_plot
    plot_results(results, field_names, model_labels, n, noise_dist);
end

end % main


% =========================================================================
%  NOISE GENERATOR
% =========================================================================
function eps = generate_noise(n, sigma, dist, df)
switch dist
    case 'normal'
        eps = sigma * randn(n,1);
    case 'uniform'
        eps = sigma * sqrt(3) * (2*rand(n,1) - 1);
    case 'laplace'
        b   = sigma / sqrt(2);
        u   = rand(n,1) - 0.5;
        eps = -b * sign(u) .* log(1 - 2*abs(u));
    case 'student'
        t_raw = trnd(df, n, 1);
        eps   = sigma * t_raw / std(t_raw);
end
end


% =========================================================================
%  DISTANCE CORRELATION — U-STATISTIC  (Székely & Rizzo 2014)
% =========================================================================
function dc = dcor_U(A, B, n)
Au = u_center(A, n);
Bu = u_center(B, n);
mask = ~eye(n);
dF   = n * (n-3);
dcov2_XY = sum(sum(Au .* Bu .* mask)) / dF;
dcov2_XX  = sum(sum(Au.^2  .* mask)) / dF;
dcov2_YY  = sum(sum(Bu.^2  .* mask)) / dF;
denom = sqrt(dcov2_XX * dcov2_YY);
if denom <= 0
    dc = 0;
else
    val = dcov2_XY / denom;
    dc  = max(0, sign(val)*sqrt(abs(val)));
end

dc=abs(dc); %just in case
end



% =========================================================================
%  PLOTTING
% =========================================================================
function plot_results(results, field_names, model_labels, n, noise_dist)

colors = [
    0.17  0.45  0.76;   % blue   - linear
    0.84  0.33  0.19;   % red    - quadratic
    0.22  0.64  0.40;   % green  - cubic
    0.58  0.27  0.69;   % purple - sin
    0.40  0.40  0.40;   % grey   - random
];

n_m = numel(field_names);

fig = figure('Name','MC Functional Dependence Simulation', ...
    'NumberTitle','off','Position',[60 60 1750 960],'Color',[0.97 0.97 0.97]);

sgtitle(sprintf(['Monte Carlo Simulation  |  n=%d  |  Noise: %s\n' ...
    'Linear  |  Quadratic  |  Cubic  |  sin  |  Random (X \x22A5 Y)'], ...
    n, noise_dist), 'FontSize',13,'FontWeight','bold','Color',[0.15 0.15 0.15]);

% Row 1: scatter
for m = 1:n_m
    fname = field_names{m};
    X = results.(fname).X;  Y = results.(fname).Y;
    Yc = results.(fname).Y_clean;
    dc = results.(fname).dcor;  pear = results.(fname).pearson_r;
    sig = results.(fname).noise_std;

    ax = subplot(3, n_m, m); hold on;
    scatter(X, Y, 9, colors(m,:),'filled','MarkerFaceAlpha',0.30);
    if ~strcmp(fname,'random')
        [Xs,idx] = sort(X);
        plot(Xs, Yc(idx), '-', 'Color', colors(m,:)*0.6, 'LineWidth', 2.0);
    end
    title(model_labels{m},'FontSize',11,'FontWeight','bold','Color',colors(m,:)*0.75);
    xlabel('X','FontSize',9); ylabel('Y','FontSize',9);
    ann = format_params(results.(fname).params, fname);
    if ~isempty(ann)
        text(0.04,0.97,ann,'Units','normalized','VerticalAlignment','top', ...
            'FontSize',7.5,'Color',[0.3 0.3 0.3],'BackgroundColor',[1 1 1 0.6]);
    end
    text(0.97,0.97,sprintf('dCor=%.3f\nr=%.3f\ns=%.2f',dc,pear,sig), ...
        'Units','normalized','VerticalAlignment','top','HorizontalAlignment','right', ...
        'FontSize',8,'FontWeight','bold','Color',colors(m,:)*0.75,'BackgroundColor',[1 1 1 0.6]);
    grid on; box on; ax.GridAlpha = 0.2; hold off;
end

% Row 2: noise / Y distributions
for m = 1:n_m
    fname = field_names{m};
    ax = subplot(3, n_m, n_m+m); hold on;
    if strcmp(fname,'random')
        histogram(results.(fname).Y, 40,'FaceColor',colors(m,:), ...
            'FaceAlpha',0.65,'EdgeColor','none','Normalization','pdf');
        xlabel('Y (indep.draw)','FontSize',8);
        title('Random: Y dist','FontSize',10,'FontWeight','bold','Color',colors(m,:)*0.75);
        text(0.97,0.97,'X \perp Y','Units','normalized','VerticalAlignment','top', ...
            'HorizontalAlignment','right','FontSize',9,'Color',[0.4 0.4 0.4], ...
            'BackgroundColor',[1 1 1 0.6]);
    else
        eps = results.(fname).noise;  sig = results.(fname).noise_std;
        histogram(eps,40,'FaceColor',colors(m,:),'FaceAlpha',0.65,'EdgeColor','none','Normalization','pdf');
        xr = linspace(min(eps),max(eps),200);
        plot(xr,normpdf(xr,0,sig),'k--','LineWidth',1.5);
        xlabel('epsilon','FontSize',9);
        title(sprintf('Noise: %s',model_labels{m}),'FontSize',10,'FontWeight','bold','Color',colors(m,:)*0.75);
        text(0.97,0.97,sprintf('s=%.3f\nskew=%.3f',std(eps),skewness(eps)), ...
            'Units','normalized','VerticalAlignment','top','HorizontalAlignment','right', ...
            'FontSize',8,'Color',[0.3 0.3 0.3],'BackgroundColor',[1 1 1 0.6]);
    end
    ylabel('Density','FontSize',9); grid on; box on; ax.GridAlpha=0.2; hold off;
end

% Row 3: correlation bar + SNR
subplot(3, n_m, [n_m*2+1, n_m*2+2, n_m*2+3]);
hold on;
pearson  = arrayfun(@(i) results.(field_names{i}).pearson_r,  1:n_m);
spearman = arrayfun(@(i) results.(field_names{i}).spearman_r, 1:n_m);
dcors    = arrayfun(@(i) results.(field_names{i}).dcor,       1:n_m);
x_pos = 1:n_m; bw = 0.25;
bar(x_pos-bw, abs(pearson),  bw,'FaceColor',[0.4 0.6 0.85],'EdgeColor','none');
bar(x_pos,    abs(spearman), bw,'FaceColor',[0.9 0.5 0.3], 'EdgeColor','none');
bar(x_pos+bw, dcors,         bw,'FaceColor',[0.3 0.7 0.4], 'EdgeColor','none');
set(gca,'XTick',x_pos,'XTickLabel',model_labels,'FontSize',10);
ylabel('|Correlation|','FontSize',10);
title('|Pearson|, |Spearman|, dCor (U-stat) — all five models','FontSize',11,'FontWeight','bold');
legend({'|Pearson r|','|Spearman rho|','dCor'},'Location','northeast','FontSize',9);
ylim([0 1.1]); grid on; box on; hold off;

subplot(3, n_m, [n_m*2+4, n_m*2+5]);
hold on;
snr_vals=[]; snr_labs={}; snr_cols=[];
for m = 1:n_m
    fname = field_names{m};
    if strcmp(fname,'random'), continue; end
    ss = std(results.(fname).Y_clean);
    ns = results.(fname).noise_std;
    if ns > 0
        snr_vals(end+1) = 20*log10(ss/ns); %#ok
    else
        snr_vals(end+1) = NaN; %#ok
    end
    snr_labs{end+1} = model_labels{m}; %#ok
    snr_cols = [snr_cols; colors(m,:)]; %#ok
end
b = bar(1:numel(snr_vals), snr_vals, 0.6,'EdgeColor','none');
b.FaceColor='flat';
for i=1:numel(snr_vals), b.CData(i,:)=snr_cols(i,:); end
set(gca,'XTick',1:numel(snr_vals),'XTickLabel',snr_labs,'FontSize',10);
ylabel('SNR (dB)','FontSize',10);
title('SNR per functional model','FontSize',11,'FontWeight','bold');
yline(0,'k--','LineWidth',1); grid on; box on; hold off;

set(findall(fig,'Type','axes'),'FontName','Helvetica');
end


% =========================================================================
%  PARAMETER ANNOTATION
% =========================================================================
function s = format_params(pr, fname)
switch fname
    case 'linear',    s = sprintf('a=%.2f, b=%.2f', pr.a, pr.b);
    case 'quadratic', s = sprintf('a=%.2f, b=%.2f, c=%.2f', pr.a, pr.b, pr.c);
    case 'cubic',     s = sprintf('a=%.2f, b=%.2f\nc=%.2f, d=%.2f', pr.a, pr.b, pr.c, pr.d);
    case 'sin',  s = sprintf('A=%.2f, B=%.2f', pr.A, pr.B);
    case 'random',    s = 'X indep. Y';
    otherwise,        s = '';
end
end


% =========================================================================
%  FIGURE  — POWER AND TYPE I ERROR
%
%  Four panels in a 2x2 layout:
%    A (top-left)  : Power curves vs sample size for each p-value method,
%                    one line per functional model (random excluded).
%                    Power = fraction of n_reps replicates with raw p < alpha.
%    B (top-right) : Type I error rate vs sample size — random model only.
%                    Horizontal dashed red line = nominal alpha.
%                    All methods shown; values should hover near alpha.
%    C (bot-left)  : Power heatmap [model x sample_size] per method,
%                    colour = empirical power, overlaid with numeric values.
%                    Three sub-panels side by side: t-test / Chi2 / Permutation.
%    D (bot-right) : Bernstein CI-based rejection rate:
%                    Rejection = CI excludes 0 (lower bound > 0).
%                    Shown for all models; random = Type I error proxy.
%
%  INPUTS:
%    res          — results struct from dcor_convergence_analysis
%    S            — style struct
%    model_names  — {'linear','quadratic','cubic','sin','random'}
%    model_labels — display labels
%    sample_sizes — vector of sample sizes tested
%    n_sizes      — numel(sample_sizes)
%
%  OUTPUTS:
%    fig  — figure handle
% =========================================================================
function fig = plot_power_type1(res, S, model_names, model_labels, ...
        sample_sizes, n_sizes)

alpha      = S.params.alpha;
n_reps     = S.params.n_reps;
n_models   = numel(model_names);
x_pos      = 1:n_sizes;
x_tick     = arrayfun(@num2str, sample_sizes, 'UniformOutput', false);
rand_idx   = find(strcmp(model_names, 'random'));
func_idx   = setdiff(1:n_models, rand_idx);

% ------------------------------------------------------------------
% Method definitions
%   name   : display label
%   field  : field in res.(mname) containing raw p-values
%   col    : line colour
%   ls     : line style
% ------------------------------------------------------------------
methods(1).name  = 'SR13 t-test';
methods(1).field = 'pval';
methods(1).col   = [0.15 0.65 0.35];
methods(1).ls    = '-o';

methods(2).name  = 'Chi2';
methods(2).field = 'pval_chi2';
methods(2).col   = [0.20 0.50 0.80];
methods(2).ls    = '-s';

methods(3).name  = 'Permutation';
methods(3).field = 'pval_perm';
methods(3).col   = [0.75 0.22 0.18];
methods(3).ls    = '-^';

n_methods = numel(methods);

% ------------------------------------------------------------------
% Pre-compute power and Type I error for each method, model, sample size
%   power_mat{k}(m, s) = fraction of reps with p < alpha
%   type1_mat{k}(s)    = same for random model
% ------------------------------------------------------------------
power_mat = cell(n_methods, 1);
type1_mat = cell(n_methods, 1);

for k = 1:n_methods
    pw = NaN(n_models, n_sizes);
    t1 = NaN(1, n_sizes);
    for m = 1:n_models
        mname = model_names{m};
        pv    = res.(mname).(methods(k).field);   % [n_reps x n_sizes]
        for s = 1:n_sizes
            col_pv = pv(:, s);
            ok     = ~isnan(col_pv);
            if sum(ok) < 2, continue; end
            pw(m, s) = mean(col_pv(ok) < alpha);
        end
    end
    t1 = pw(rand_idx, :);
    power_mat{k} = pw;
    type1_mat{k} = t1;
end

% ------------------------------------------------------------------
% Bernstein CI rejection rate: reject H0 if lower CI bound > 0
%   Uses res.(mname).ci_bern  [n_sizes x 2] (distribution-level CI)
%   For a per-replicate version we use the individual dCor values and
%   recompute per-replicate Bernstein CIs on the fly.
% ------------------------------------------------------------------
bern_reject = NaN(n_models, n_sizes);   % rejection rate

for m = 1:n_models
    mname = model_names{m};
    for s = 1:n_sizes
        ns    = sample_sizes(s);
        d     = res.(mname).dcor_u(:, s);
        ok    = ~isnan(d);
        if sum(ok) < 2, continue; end
        d_ok  = d(ok);
        n_ok  = numel(d_ok);

        % Per-replicate Bernstein CI lower bound
        % Maurer & Pontil (2009): eps = sqrt(2*s2*ln(2/alpha)/K) + (2/3)*ln(2/alpha)/K
        % where K = ns (sample size), s2 = variance of the *statistic* estimate.
        % Single-observation case: s2 = 0, so eps = (2/3)*ln(2/alpha)/ns
        lnt  = log(2 / alpha);
        eps_vec = (2/3) * lnt ./ ns;                    % scalar, same for all reps
        lo_bern = d_ok - eps_vec;                        % [n_ok x 1]

        % Reject if lo > 0  (CI entirely above zero)
        bern_reject(m, s) = mean(lo_bern > 0.1); %&&& change to 0.1 from 0
    end
end

% ------------------------------------------------------------------
% Figure layout
% ------------------------------------------------------------------
fig = figure('Name', 'Fig8_Power_Type1', 'NumberTitle', 'off', ...
    'Color', 'w', 'Units', 'centimeters', 'Position', [1 1 26 26]);

sgtitle(sprintf(['Power and Type I Error  |  dCor (U-stat)  |  \\alpha=%.3f  ' ...
    '|  %d replicates'], alpha, n_reps), ...
    'FontName', S.font, 'FontSize', S.fsz_title + 1, ...
    'FontWeight', 'bold', 'Color', S.black, 'Interpreter', 'tex');

% ---- Panel A: Power curves (functional models only) ---- 
ax_A = subplot(2, 2, 1);
hold(ax_A, 'on');

% Nominal alpha reference
yline(ax_A, alpha, '--', 'Color', [0.82 0.05 0.05], 'LineWidth', 1.0, ...
    'Label', sprintf('\\alpha=%.3f', alpha), ...
    'LabelHorizontalAlignment', 'right', ...
    'FontName', S.font, 'FontSize', S.fsz_ax, ...
    'HandleVisibility', 'off', 'Interpreter', 'tex');
yline(ax_A, 0.80, ':', 'Color', S.grey, 'LineWidth', 0.7, ...
    'HandleVisibility', 'off');   % 80% power guideline

% One linestyle per method, one hue per model — use marker shape for model
model_markers = {'o', 's', '^', 'd'};   % one per functional model
h_legend = [];
leg_labels = {};

for k = 1:n_methods
    for mi = 1:numel(func_idx)
        m     = func_idx(mi);
        pw    = power_mat{k}(m, :);
        valid = ~isnan(pw);
        if ~any(valid), continue; end

        % Blend method colour with model colour for unique appearance
        blend_col = methods(k).col; % try method color instead 0.55 * methods(k).col + 0.45 * S.mcols(m, :);

        mk  = [methods(k).ls(1), model_markers{mi}];   % e.g. '-o'
        h   = plot(ax_A, sample_sizes(valid), pw(valid), mk, ...
            'Color', blend_col, 'LineWidth', S.lw, ...
            'MarkerSize', S.ms + 0.5, ...
            'MarkerFaceColor', blend_col, ...
            'DisplayName', sprintf('%s / %s', methods(k).name, model_labels{m}));
        h_legend(end+1)  = h;         %#ok
        leg_labels{end+1} = sprintf('%s – %s', methods(k).name, model_labels{m}); %#ok
    end
end

% Power = 1 reference line
yline(ax_A, 1.0, '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 0.6, ...
    'HandleVisibility', 'off');

set(ax_A, 'XScale', 'log', 'XTick', sample_sizes, 'XTickLabel', x_tick);
xlim(ax_A, [sample_sizes(1)*0.8, sample_sizes(end)*1.3]);
ylim(ax_A, [-0.03, 1.05]);
xlabel(ax_A, 'Sample size  n', 'FontName', S.font, 'FontSize', S.fsz);
ylabel(ax_A, 'Empirical power  (P[p < \alpha])', ...
    'FontName', S.font, 'FontSize', S.fsz, 'Interpreter', 'tex');
title(ax_A, 'A   Power: functional models (H_1 true)', ...
    'FontName', S.font, 'FontSize', S.fsz_title, ...
    'FontWeight', 'bold', 'Interpreter', 'tex');
legend(ax_A, h_legend, leg_labels, 'Location', 'southeast', ...
    'FontName', S.font, 'FontSize', 5.5, 'Box', 'off', ...
    'NumColumns', 2);

annotation_text = sprintf('Dashed red = \\alpha  |  Dotted = 80%% power guideline');
text(ax_A, 0.02, 0.03, annotation_text, 'Units', 'normalized', ...
    'FontName', S.font, 'FontSize', S.fsz_ax - 1, ...
    'Color', S.grey, 'Interpreter', 'tex');
pub_ax(ax_A, S);
hold(ax_A, 'off');

% ---- Panel B: Type I error (random model) ----
ax_B = subplot(2, 2, 2);
hold(ax_B, 'on');

% Nominal alpha band: +/- 2 * sqrt(alpha*(1-alpha)/n_reps) — 95% MC interval
mc_se    = 2 * sqrt(alpha * (1-alpha) / n_reps);
patch(ax_B, [sample_sizes(1)*0.8, sample_sizes(end)*1.3, ...
             sample_sizes(end)*1.3, sample_sizes(1)*0.8], ...
    [alpha - mc_se, alpha - mc_se, alpha + mc_se, alpha + mc_se], ...
    [0.85 0.85 0.85], 'FaceAlpha', 0.45, 'EdgeColor', 'none', ...
    'HandleVisibility', 'off');

yline(ax_B, alpha, '--', 'Color', [0.82 0.05 0.05], 'LineWidth', 1.2, ...
    'Label', sprintf('\\alpha=%.3f', alpha), ...
    'LabelHorizontalAlignment', 'right', ...
    'FontName', S.font, 'FontSize', S.fsz_ax, ...
    'HandleVisibility', 'off', 'Interpreter', 'tex');

for k = 1:n_methods
    t1    = type1_mat{k};
    valid = ~isnan(t1);
    if ~any(valid), continue; end
    plot(ax_B, sample_sizes(valid), t1(valid), methods(k).ls, ...
        'Color', methods(k).col, 'LineWidth', S.lw + 0.2, ...
        'MarkerSize', S.ms + 1, 'MarkerFaceColor', methods(k).col, ...
        'DisplayName', methods(k).name);
end

% Bernstein Type I error for random model
bern_t1 = bern_reject(rand_idx, :);
valid   = ~isnan(bern_t1);
if any(valid)
    plot(ax_B, sample_sizes(valid), bern_t1(valid), '-p', ...
        'Color', S.bern_col, 'LineWidth', S.lw + 0.2, ...
        'MarkerSize', S.ms + 1, 'MarkerFaceColor', S.bern_col, ...
        'DisplayName', 'Bernstein CI');
end

set(ax_B, 'XScale', 'log', 'XTick', sample_sizes, 'XTickLabel', x_tick);
xlim(ax_B, [sample_sizes(1)*0.8, sample_sizes(end)*1.3]);
ylim(ax_B, [0,0.5]); % min(1, alpha * 6)]);   % zoom around alpha &&&
xlabel(ax_B, 'Sample size  n', 'FontName', S.font, 'FontSize', S.fsz);
ylabel(ax_B, 'Type I error rate  (H_0 true)', ...
    'FontName', S.font, 'FontSize', S.fsz, 'Interpreter', 'tex');
title(ax_B, 'B   Type I error: random model (H_0 true)', ...
    'FontName', S.font, 'FontSize', S.fsz_title, ...
    'FontWeight', 'bold', 'Interpreter', 'tex');
legend(ax_B, 'Location', 'northeast', 'FontName', S.font, ...
    'FontSize', S.fsz_ax, 'Box', 'off');

% Annotate with actual values
for k = 1:n_methods
    t1    = type1_mat{k};
    valid = ~isnan(t1);
    if ~any(valid), continue; end
    last_s = find(valid, 1, 'last');
    text(ax_B, sample_sizes(last_s) * 1.05, t1(last_s), ...
        sprintf('%.3f', t1(last_s)), ...
        'FontName', S.font, 'FontSize', S.fsz_ax - 1, ...
        'Color', methods(k).col * 0.8, 'VerticalAlignment', 'middle');
end

text(ax_B, 0.02, 0.05, 'Grey band = ±2 MC SE of nominal \alpha', ...
    'Units', 'normalized', 'FontName', S.font, 'FontSize', S.fsz_ax - 1, ...
    'Color', S.grey, 'Interpreter', 'tex');
pub_ax(ax_B, S);
ax_B.Color = [0.96 0.95 0.92];   % random model background
hold(ax_B, 'off');

% ---- Panel C: Power heatmaps (one per method, functional models only) ----
ax_C = subplot(2, 2, 3);
hold(ax_C, 'on');

%  tile 3 sub-heatmaps left-to-right inside this axes using imagesc
% with manual x offsets; simpler to use a dedicated tiled approach.
% Build a combined [n_func_models x n_sizes x n_methods] array and
% display as a single wide heatmap with model-group separators.

n_func   = numel(func_idx);
% Concatenate vertically: rows = [ttest models; chi2 models; perm models]
heat_mat = NaN(n_methods * n_func, n_sizes);

for k = 1:n_methods
    rows = (k-1)*n_func + (1:n_func);
    for mi = 1:n_func
        m = func_idx(mi);
        heat_mat(rows(mi), :) = power_mat{k}(m, :);
    end
end

imagesc(ax_C, heat_mat);
colormap(ax_C, flipud(gray(256)));   % dark = high power
clim(ax_C, [0, 1]);
cb_C = colorbar(ax_C, 'eastoutside');
cb_C.Label.String = 'Power (fraction p < \alpha)';
cb_C.Label.FontSize = S.fsz_ax;
cb_C.FontSize = S.fsz_ax;

% Y-axis labels: method + model
y_labs = cell(n_methods * n_func, 1);
for k = 1:n_methods
    for mi = 1:n_func
        m = func_idx(mi);
        y_labs((k-1)*n_func + mi) = ...
            {sprintf('%s|%s', methods(k).name, model_labels{m})};
    end
end
set(ax_C, 'YTick', 1:(n_methods*n_func), 'YTickLabel', y_labs, ...
    'XTick', x_pos, 'XTickLabel', x_tick, ...
    'TickLength', [0 0], 'FontSize', S.fsz_ax - 1);

% Method group separator lines
for k = 1:n_methods - 1
    yline(ax_C, k * n_func + 0.5, '-', 'Color', [0.5 0.5 0.5], ...
        'LineWidth', 1.2, 'HandleVisibility', 'off');
end

% Overlay numeric values
for r = 1:(n_methods * n_func)
    for s = 1:n_sizes
        v = heat_mat(r, s);
        if isnan(v), continue; end
        % White text on dark cells, black on light
        tcol = [1 1 1] * double(v < 0.55);
        text(ax_C, s, r, sprintf('%.2f', v), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
            'FontName', S.font, 'FontSize', max(4, S.fsz_ax - 2), ...
            'Color', tcol);
    end
end

xlabel(ax_C, 'Sample size  n', 'FontName', S.font, 'FontSize', S.fsz);
title(ax_C, 'C   Power heatmap: all methods × functional models', ...
    'FontName', S.font, 'FontSize', S.fsz_title, 'FontWeight', 'bold');
pub_ax(ax_C, S);
set(ax_C, 'Box', 'on', 'TickDir', 'in');
hold(ax_C, 'off');

% ---- Panel D: Bernstein CI rejection rate for all models ----
ax_D = subplot(2, 2, 4);
hold(ax_D, 'on');

yline(ax_D, alpha, '--', 'Color', [0.82 0.05 0.05], 'LineWidth', 1.2, ...
    'Label', sprintf('\\alpha=%.3f  (Type I ref)', alpha), ...
    'LabelHorizontalAlignment', 'right', ...
    'FontName', S.font, 'FontSize', S.fsz_ax, ...
    'HandleVisibility', 'off', 'Interpreter', 'tex');
yline(ax_D, 0.80, ':', 'Color', S.grey, 'LineWidth', 0.7, ...
    'HandleVisibility', 'off');

for m = 1:n_models
    br    = bern_reject(m, :);
    valid = ~isnan(br);
    if ~any(valid), continue; end
    col   = S.mcols(m, :);
    lw_m  = S.lw + 0.3 * double(m == rand_idx);   % thicker for random
    ls_m  = '-o';
    if m == rand_idx, ls_m = '-p'; end
    plot(ax_D, sample_sizes(valid), br(valid), ls_m, ...
        'Color', col, 'LineWidth', lw_m, ...
        'MarkerSize', S.ms + double(m == rand_idx), ...
        'MarkerFaceColor', col, ...
        'DisplayName', model_labels{m});
end

set(ax_D, 'XScale', 'log', 'XTick', sample_sizes, 'XTickLabel', x_tick);
xlim(ax_D, [sample_sizes(1)*0.8, sample_sizes(end)*1.3]);
ylim(ax_D, [-0.03, 1.05]);
xlabel(ax_D, 'Sample size  n', 'FontName', S.font, 'FontSize', S.fsz);
ylabel(ax_D, 'Rejection rate  (CI_{lo} > 0.1)', ...
    'FontName', S.font, 'FontSize', S.fsz);
title(ax_D, 'D   Bernstein CI rejection rate (all models)', ...
    'FontName', S.font, 'FontSize', S.fsz_title, ...
    'FontWeight', 'bold');
legend(ax_D, 'Location', 'southeast', 'FontName', S.font, ...
    'FontSize', S.fsz_ax, 'Box', 'off');

text(ax_D, 0.02, 0.05, ...
    'Reject H_0 when Bernstein CI lower bound > 0', ...
    'Units', 'normalized', 'FontName', S.font, ...
    'FontSize', S.fsz_ax - 1, 'Color', S.grey, 'Interpreter', 'tex');

ax_D.Color = [0.96 0.95 0.92];   % highlight random in background
pub_ax(ax_D, S);
hold(ax_D, 'off');

% ------------------------------------------------------------------
% Footer annotation
% ------------------------------------------------------------------
annotation(fig, 'textbox', [0.02 0.001 0.96 0.028], ...
    'String', [ ...
    'A: Power curves — each line = one method × functional model. ', ...
    'B: Type I error — random model (H0 true); all values should lie near \alpha (red dashed); grey band = ±2 MC SE. ', ...
    'C: Power heatmap — dark = high power; rows grouped by method; separator lines between method blocks. ', ...
    'D: Bernstein CI rejection rate — H0 rejected when lower CI bound > 0; random (pentagon) = Type I proxy.'], ...
    'FontName', S.font, 'FontSize', 6, 'EdgeColor', 'none', ...
    'Color', S.grey, 'HorizontalAlignment', 'center', 'Interpreter', 'tex');

set(fig, 'PaperUnits', 'centimeters', 'PaperSize', [26 26]);
end
