% combined analysis for two samples in the input for all things and then
% adding difference plots for degree and betweenness 


%two_samples_second_analysis("C:\Project\CorrSize_Inflammation\NC_data.csv", ...
%    "C:\Project\CorrSize_Inflammation\AD_data.csv",...
%    "C:\Project\CorrSize_Inflammation\",...
%    0.01, 0.1, 1, 1,22,50,100)

function two_samples_second_analysis(Filename1,Filename2, outputdirectory, alpha, ci_min, Normalize, Imputer,select_label,NSubset,NBoot, do_plot, layout, seed, B_perm)




% -------------------------------------------------------------------------
% Defaults
% -------------------------------------------------------------------------
if nargin < 3  || isempty(outputdirectory)
    out_dir = fullfile(pwd, 'dcor_figures');
    if ~exist(out_dir, 'dir'), mkdir(out_dir); end
else
    out_dir = outputdirectory;
end



if nargin < 4  || isempty(alpha),     alpha     = 0.05;  end
if nargin < 5  || isempty(ci_min),    ci_min    = 0.10;  end
if nargin < 6  || isempty(Normalize), Normalize = 1;     end
if nargin < 7  || isempty(Imputer),   Imputer   = 1;     end
if nargin < 8  || isempty(select_label),   select_label   = 1;  end

if nargin < 9  || isempty(NSubset),   NSubset   = 20;  end
if nargin < 10  || isempty(NBoot),    NBoot    = 20; end

if nargin < 11  || isempty(do_plot),   do_plot   = true;  end
if nargin < 12  || isempty(layout),    layout    = 'force'; end
if nargin < 13  || isempty(seed),      seed      = 1;     end
if nargin < 14 || isempty(B_perm),    B_perm    = 100;   end

rng(seed);

% Pack into struct
pm.alpha     = str2double(alpha); %
pm.ci_min    = str2double(ci_min);%
pm.Normalize = Normalize;
pm.Imputer   = Imputer;
pm.plot      = do_plot;
pm.layout    = layout;
pm.seed      = seed;
pm.B_perm    = B_perm;
pm.Filename  =Filename1;
pm.select_label=str2double(select_label);%
pm.NSubset=str2double(NSubset);
pm.NBoot=str2double(NBoot);


net1 = dcor_network_single_sample_second_analysis(Filename1, outputdirectory,pm,"Set1")

pm.Filename  =Filename2;

net2 = dcor_network_single_sample_second_analysis(Filename2, outputdirectory, pm,"Set2")


fig_colors = struct(...
    'fdr',    [0.20 0.50 0.80], ...
    'ci',     [0.15 0.65 0.35], ...
    'fdr_ci', [0.72 0.15 0.15]);
threshold_types  = {'fdr','ci','fdr_ci'};
threshold_labels = {sprintf('FDR q<%.2f', pm.alpha), ...
                    sprintf('CI lo>%.2f', pm.ci_min), ...
                    sprintf('FDR q<%.2f AND CI lo>%.2f', pm.alpha, pm.ci_min)};

comparison = compare_networks(net1, net2);




end


function net1 = dcor_network_single_sample_second_analysis(Filename, outputdirectory, pm,set)
% select_label is the number of metabolite from the drop down list
warning('OFF')
opts = detectImportOptions(Filename, 'ReadVariableNames', true, 'PreserveVariableNames', true);
opts.DataLines = [2, Inf];
opts.Delimiter = ",";
inputfile = readtable(Filename, opts);
clear opts

% -------------------------------------------------------------------------
% Defaults
% -------------------------------------------------------------------------
if nargin < 2  || isempty(outputdirectory)
    out_dir = fullfile(pwd, 'dcor_figures');
    if ~exist(out_dir, 'dir'), mkdir(out_dir); end
else
    out_dir = outputdirectory;
end
seed=pm.seed;

rng(seed);

% Pack into struct, rest already in structure
pm.set       =set;

% -------------------------------------------------------------------------
% Convert table to matrix
% -------------------------------------------------------------------------
col_names = {};
if istable(inputfile)
    col_names = inputfile.Properties.VariableNames(2:end)';
    X1 = table2array(inputfile(:, 2:end));
end
if ~isa(X1, 'double'), X1 = double(X1); end

if pm.Imputer == 1
    X1 = knnimpute(X1, 5);
end
if pm.Normalize == 1
    X1 = normalize(X1);
end

% Save preprocessed data
filenameout14 = 'Inputdata_with_preprocessing.xlsx';
output_values = inputfile;
writetable(output_values, fullfile(out_dir, filenameout14), 'Sheet', 'Original Input');
output_values{:, 2:end} = X1;
writetable(output_values, fullfile(out_dir, filenameout14), 'Sheet', 'With preprocessing');

[n1, p1] = size(X1);
p = p1;
fprintf('  Input dimensions: X1=[%d obs x %d feat]\n', n1, p1);

% Feature labels
if ~isempty(col_names)
    labels = col_names;
    assert(numel(labels) == p, ...
        'Column names have %d entries but table has %d columns.', numel(labels), p);
else
    labels = arrayfun(@(k) sprintf('F%d', k), 1:p, 'UniformOutput', false);
end
pm.labels = labels;

fprintf('  Features (%d): %s ... %s\n', p, labels{1}, labels{end});
fprintf('  alpha=%.3f  ci_min=%.3f\n', pm.alpha, pm.ci_min);

% -------------------------------------------------------------------------
% Compute network
% -------------------------------------------------------------------------
fprintf('Network 1 (n=%d observations):\n', n1);
net1 = compute_network(X1, n1, p, labels, pm, out_dir);
net1 = add_graph_metrics(net1);

if pm.plot
    plot_all(net1, out_dir, pm);
end

fprintf('\nDone.\n');
end


% =========================================================================
%  COMPUTE_NETWORK
% =========================================================================
function net = compute_network(X, n, p, labels, pm, out_dir)
if istable(X), X = table2array(X); end
if ~isa(X, 'double'), X = double(X); end
[Xr, Xc] = size(X);
assert(Xr == n, 'compute_network: X has %d rows but expected n=%d.', Xr, n);
assert(Xc == p, 'compute_network: X has %d columns but expected p=%d.', Xc, p);
alpha=pm.alpha;
ci_min=pm.ci_min;
% --- bootstrap parameters ---
n_boot    = pm.NBoot;
n_sub=pm.NSubset;
assert(n >= n_sub, 'Need at least %d samples, have %d.', n_sub, n);

n_pairs   = p*(p-1)/2;
fprintf('  Computing %d pairwise dCor values (%d boots x n=%d subset)...\n', ...
    n_pairs, n_boot, n_sub);

% accumulate dcor across boots
dcor_boot = zeros(p, p, n_boot);

for b = 1:n_boot
    idx  = randperm(n, n_sub);
    Xsub = X(idx, :);

    % u-centered distance matrices for subset
    D_sub = cell(p, 1);
    U_sub = cell(p, 1);
    for k = 1:p
        D_sub{k} = squareform(pdist(Xsub(:,k)));
        U_sub{k} = u_center(D_sub{k}, n_sub);
    end

    for i = 1:p
        for j = i+1:p
            dc = dcor_from_ucentered(U_sub{i}, U_sub{j}, n_sub);
            dcor_boot(i,j,b) = dc;
            dcor_boot(j,i,b) = dc;
        end
    end
    fprintf('    boot %2d/%d done\n', b, n_boot);
end

% --- mean dCor across boots ---
dcor_mat = mean(dcor_boot, 3);

% --- p-value and CI from mean dCor (treating n_sub as effective n) ---
pval_mat  = ones(p, p);
ci_lo_mat = zeros(p, p);
ci_hi_mat = zeros(p, p);
pval_vec  = NaN(n_pairs, 1);
pair_idx  = zeros(n_pairs, 2);
pair_count = 0;
% Use full n for p-value (recommended)
%The mean dCor is estimated from all data implicitly — use n for the test, n_sub only controls subsample noise:
for i = 1:p
    for j = i+1:p
        dc = dcor_mat(i,j);
        pv = dcor_chi2_test(dc, n);          % use n_sub as effective n
        ci = bernstein_ci_dcor(dc, n, alpha);
        pval_mat(i,j) = pv;  pval_mat(j,i) = pv;
        ci_lo_mat(i,j) = ci(1);  ci_lo_mat(j,i) = ci(1);
        ci_hi_mat(i,j) = ci(2);  ci_hi_mat(j,i) = ci(2);
        pair_count = pair_count + 1;
        pval_vec(pair_count) = pv;
        pair_idx(pair_count,:) = [i, j];
    end
end

% --- FDR correction ---
[~, ~, ~, qval_vec] = fdr_bh(pval_vec, alpha);
qval_mat = ones(p, p);
for k = 1:pair_count
    i = pair_idx(k,1);  j = pair_idx(k,2);
    qval_mat(i,j) = qval_vec(k);
    qval_mat(j,i) = qval_vec(k);
end

% --- adjacency matrices ---
adj_full   = double(dcor_mat > 0);
adj_fdr    = double(qval_mat  < alpha);
adj_ci     = double(ci_lo_mat > ci_min);
adj_fdr_ci = double(adj_fdr & adj_ci);
for k = 1:p
    adj_full(k,k)=0; adj_fdr(k,k)=0; adj_ci(k,k)=0; adj_fdr_ci(k,k)=0;
end

n_fdr   = sum(adj_fdr(:))/2;
n_ci    = sum(adj_ci(:))/2;
n_fdrci = sum(adj_fdr_ci(:))/2;
fprintf('  Edges — full: %d | FDR: %d | CI: %d | FDR+CI: %d\n', ...
    n_pairs, n_fdr, n_ci, n_fdrci);

% --- pack output ---
net.dcor       = dcor_mat;
net.dcor_boot  = dcor_boot;          % p x p x n_boot, if you need it later
net.dcor_std   = std(dcor_boot, 0, 3);  % boot standard deviation
net.pval       = pval_mat;
net.qval       = qval_mat;
net.ci_lo      = ci_lo_mat;
net.ci_hi      = ci_hi_mat;
net.adj_full   = adj_full;
net.adj_fdr    = adj_fdr;
net.adj_ci     = adj_ci;
net.adj_fdr_ci = adj_fdr_ci;
net.labels     = labels;
net.n_features = p;
net.alpha      = alpha;
net.ci_min     = ci_min;
net.n          = n;
net.n_boot     = n_boot;
net.n_sub      = n_sub;

% --- Excel output ---
filenameout24 = 'Single_Sample_UnbiasedDistanceCorrelation.xlsx';
write_net_matrix(net, 'dcor',     'Unbiased dCorr (mean)', filenameout24, out_dir);
write_net_matrix(net, 'dcor_std', 'dCorr boot SD',         filenameout24, out_dir);
write_net_matrix(net, 'pval',     'Chi2 p-value',          filenameout24, out_dir);
write_net_matrix(net, 'qval',     'FDR q-value',           filenameout24, out_dir);
write_net_matrix(net, 'ci_lo',    'Low limit CI',          filenameout24, out_dir);
write_net_matrix(net, 'ci_hi',    'High limit CI',         filenameout24, out_dir);
end



% =========================================================================
%  WRITE_NET_MATRIX  — helper to write one matrix sheet to Excel
% =========================================================================
function write_net_matrix(net, field, sheet_name, filename, out_dir)
T = array2table(net.(field));
T.Properties.VariableNames = string(net.labels');
name_col = table(string(net.labels), 'VariableNames', {'Name'});
T = [name_col, T];
writetable(T, fullfile(out_dir, filename), 'Sheet', sheet_name);
end


% =========================================================================
%  ADD_GRAPH_METRICS
% =========================================================================
function net = add_graph_metrics(net)
types = {'full','fdr','ci','fdr_ci'};
for t = 1:numel(types)
    adj = net.(['adj_' types{t}]);
    net.(['degree_' types{t}]) = sum(adj, 2);
    net.(['bc_'     types{t}]) = betweenness_centrality(adj);
end
end


% =========================================================================
%  COMPARE_NETWORKS
% =========================================================================
function comp = compare_networks(net1, net2)
comp = struct();
types = {'full','fdr','ci','fdr_ci'};

assert(net1.n_features == net2.n_features, ...
    'compare_networks: net1 has %d features, net2 has %d. Must be equal.', ...
    net1.n_features, net2.n_features);
assert(isequal(net1.labels, net2.labels), ...
    'compare_networks: feature labels differ between net1 and net2.');

comp.shared_labels = net1.labels;

for t = 1:numel(types)
    d1 = net1.(['degree_' types{t}]);
    d2 = net2.(['degree_' types{t}]);
    b1 = net1.(['bc_'     types{t}]);
    b2 = net2.(['bc_'     types{t}]);
    comp.(['delta_degree_' types{t}]) = d1 - d2;
    comp.(['delta_bc_'     types{t}]) = b1 - b2;
end
end


% =========================================================================
%  PLOT_ALL
% =========================================================================
function plot_all(net1, out_dir, pm)

fig_colors = struct(...
    'fdr',    [0.20 0.50 0.80], ...
    'ci',     [0.15 0.65 0.35], ...
    'fdr_ci', [0.72 0.15 0.15]);
threshold_types  = {'fdr','ci','fdr_ci'};
threshold_labels = {sprintf('FDR q<%.2f', pm.alpha), ...
                    sprintf('CI lo>%.2f', pm.ci_min), ...
                    sprintf('FDR q<%.2f AND CI lo>%.2f', pm.alpha, pm.ci_min)};

fprintf('\n  Saving figures to: %s\n', out_dir);

%fig1 = plot_heatmaps(net1);
%save_figure(fig1, out_dir, 'fig1_dcor_heatmaps', net1);

%fig3 = plot_degree_comparison(net1, threshold_types, threshold_labels, fig_colors);
%save_figure(fig3, out_dir, 'fig3_degree_comparison', net1);

%fig4 = plot_bc_comparison(net1, threshold_types, threshold_labels, fig_colors);
%save_figure(fig4, out_dir, 'fig4_betweenness_centrality', net1);

%fig5 = plot_edge_summary(net1, pm);
%save_figure(fig5, out_dir, 'fig5_edge_summary', net1);


% ---- Figure N: Network graphs (one column per threshold type) ----------
figN = plot_networks(net1,  threshold_types, threshold_labels, ...
    fig_colors, pm);

filenamefig = sprintf('figN_networks_%s' , pm.set);

save_figure(figN, out_dir, filenamefig, net1);


fprintf('  All figures saved.\n');
end


% =========================================================================
%  FIGURE 1 — HEATMAPS
% =========================================================================
function fig = plot_heatmaps(net1)
% Create figure without displaying it with ,'Visible', 'off');

fig = figure('Name','dCor Heatmaps','NumberTitle','off', ...
    'Position',[30 30 1600 780],'Color',[0.97 0.97 0.97],'Visible', 'off');
sgtitle('Pairwise dCor (U-stat)  |  p-value  |  Bernstein CI width', ...
    'FontSize',12,'FontWeight','bold');

net = net1;
p   = net.n_features;
lbl = net.labels;

ax1 = subplot(1, 4, 1);
imagesc(ax1, net.dcor, 'AlphaData', ~isnan(net.dcor));
set(ax1, 'Color', 'white');
colorbar(ax1); clim(ax1,[0 1]);
colormap(ax1, parula(256));
set(ax1,'XTick',1:p,'XTickLabel',lbl,'YTick',1:p,'YTickLabel',lbl, ...
    'XTickLabelRotation',45,'FontSize',6,'TickLength',[0 0]);
title(ax1,'Set 1  —  dCor','FontSize',9,'FontWeight','bold');
axis(ax1,'square');

tempdcor = net.dcor; tempdcor(tempdcor==0 | net.adj_fdr==0) = nan;
ax2 = subplot(1, 4, 2);
imagesc(ax2, tempdcor, 'AlphaData', ~isnan(tempdcor));
set(ax2, 'Color', 'white');
colorbar(ax2); clim(ax2,[0 1]);
colormap(ax2, parula(256));
set(ax2,'XTick',1:p,'XTickLabel',lbl,'YTick',1:p,'YTickLabel',lbl, ...
    'XTickLabelRotation',45,'FontSize',6,'TickLength',[0 0]);
title(ax2,'Set 1  —  dCor (FDR sig.)','FontSize',9,'FontWeight','bold');
axis(ax2,'square');

tempdcor = net.dcor; tempdcor(tempdcor==0 | net.adj_ci==0) = nan;
ax3 = subplot(1, 4, 3);
imagesc(ax3, tempdcor, 'AlphaData', ~isnan(tempdcor));
set(ax3, 'Color', 'white');
colorbar(ax3); clim(ax3,[0 1]);
colormap(ax3, parula(256));
set(ax3,'XTick',1:p,'XTickLabel',lbl,'YTick',1:p,'YTickLabel',lbl, ...
    'XTickLabelRotation',45,'FontSize',6,'TickLength',[0 0]);
title(ax3,'Set 1  —  dCor Bernstein CI pass','FontSize',9,'FontWeight','bold');
axis(ax3,'square');

tempdcor = net.dcor; tempdcor(tempdcor==0 | net.adj_fdr_ci==0) = nan;
ax4 = subplot(1, 4, 4);
imagesc(ax4, tempdcor, 'AlphaData', ~isnan(tempdcor));
set(ax4, 'Color', 'white');
colorbar(ax4); clim(ax4,[0 1]);
colormap(ax4, parula(256));
set(ax4,'XTick',1:p,'XTickLabel',lbl,'YTick',1:p,'YTickLabel',lbl, ...
    'XTickLabelRotation',45,'FontSize',6,'TickLength',[0 0]);
title(ax4,'Set 1  —  dCor FDR+ Bernstein CI pass','FontSize',9,'FontWeight','bold');
axis(ax4,'square');

annotation(fig,'textbox',[0.01 0.002 0.98 0.025], ...
    'String',['Diagonal zeroed  |  q-values: Benjamini-Hochberg FDR  |  ' ...
              'CI: Bernstein (Maurer & Pontil 2009)'], ...
    'FontSize',7.5,'EdgeColor','none','HorizontalAlignment','center');
set(findall(fig,'Type','axes'),'FontName','Helvetica');
end


% =========================================================================
%  FIGURE 3 — DEGREE COMPARISON
% =========================================================================
function fig = plot_degree_comparison(net1, thresh_types, thresh_labels, fig_colors)

n_thresh = numel(thresh_types);
n_shared = numel(net1.labels);
x_pos    = 1:n_shared;

fig = figure('Name','Degree Information','NumberTitle','off', ...
    'Position',[50 50 1600 900],'Color',[0.97 0.97 0.97],'Visible', 'off');
sgtitle('Node Degree: Set 1', 'FontSize',11,'FontWeight','bold');

for t = 1:n_thresh
    col = fig_colors.(thresh_types{t});
    d1  = net1.(['degree_' thresh_types{t}]);

    ax1 = subplot(1, n_thresh, t);
    hold(ax1,'on');
    bar(ax1, x_pos, d1, 0.35, 'FaceColor', col, 'EdgeColor','none', ...
        'FaceAlpha',0.85, 'DisplayName','Set 1');
    hold(ax1,'off');

    ax1.XLim               = [min(x_pos)-0.5, max(x_pos)+0.5];
    ax1.XTickMode          = 'manual';
    ax1.XTickLabelMode     = 'manual';
    ax1.XTick              = x_pos;
    ax1.XTickLabel         = net1.labels;
    ax1.XTickLabelRotation = 45;
    ax1.FontSize           = 7.5;
    ax1.TickLabelInterpreter = 'none';

    ylabel(ax1, 'Degree', 'FontSize', 9);
    title(ax1, thresh_labels{t}, 'FontSize',9, 'FontWeight','bold', 'Color',col*0.75);
    legend(ax1, 'Degree', 'Location','northwest', 'FontSize',7, 'Box','off');
    grid(ax1,'on'); box(ax1,'on'); ax1.GridAlpha = 0.13;
end

annotation(fig,'textbox',[0.01 0.002 0.98 0.025], ...
    'String','Bars = Node Degree', ...
    'FontSize',7.5,'EdgeColor','none','HorizontalAlignment','center');
set(findall(fig,'Type','axes'),'FontName','Helvetica');
end


% =========================================================================
%  FIGURE 4 — BETWEENNESS CENTRALITY
% =========================================================================
function fig = plot_bc_comparison(net1, thresh_types, thresh_labels, fig_colors)

n_thresh = numel(thresh_types);
n_shared = numel(net1.labels);
x_pos    = 1:n_shared;

fig = figure('Name','Betweenness Centrality','NumberTitle','off', ...
    'Position',[60 60 1600 900],'Color',[0.97 0.97 0.97],'Visible', 'off');
sgtitle('Betweenness Centrality: Set 1', 'FontSize',11,'FontWeight','bold');

for t = 1:n_thresh
    col = fig_colors.(thresh_types{t});
    b1  = net1.(['bc_' thresh_types{t}]);

    ax1 = subplot(1, n_thresh, t);
    hold(ax1,'on');
    bar(ax1, x_pos, b1, 0.35, 'FaceColor', col, 'EdgeColor','none', ...
        'FaceAlpha',0.85, 'DisplayName','Set 1');
    hold(ax1,'off');

    ax1.XLim               = [min(x_pos)-0.5, max(x_pos)+0.5];
    ax1.XTickMode          = 'manual';
    ax1.XTickLabelMode     = 'manual';
    ax1.XTick              = x_pos;
    ax1.XTickLabel         = net1.labels;
    ax1.XTickLabelRotation = 45;
    ax1.FontSize           = 7.5;
    ax1.TickLabelInterpreter = 'none';

    ylabel(ax1, 'Betweenness', 'FontSize', 9);
    title(ax1, thresh_labels{t}, 'FontSize',9, 'FontWeight','bold', 'Color',col*0.75);
    legend(ax1, 'Betweenness', 'Location','northwest', 'FontSize',7, 'Box','off');
    grid(ax1,'on'); box(ax1,'on'); ax1.GridAlpha = 0.13;
end

annotation(fig,'textbox',[0.01 0.002 0.98 0.025], ...
    'String',['Betweenness centrality = normalised count of shortest paths through each node  |  ' ...
              'Higher = more central / bridge node'], ...
    'FontSize',7.5,'EdgeColor','none','HorizontalAlignment','center');
set(findall(fig,'Type','axes'),'FontName','Helvetica');
end


% =========================================================================
%  FIGURE 5 — EDGE SUMMARY
% =========================================================================
function fig = plot_edge_summary(net1, pm)

fig = figure('Name','Edge Summary','NumberTitle','off', ...
    'Position',[70 70 1400 600],'Color',[0.97 0.97 0.97],'Visible', 'off');
sgtitle(sprintf('Edge-level summary: dCor vs significance  |  alpha=%.2f  |  ci_min=%.2f', ...
    pm.alpha, pm.ci_min), 'FontSize',11,'FontWeight','bold');

net = net1;
p   = net.n_features;

dc_vec  = [];  q_vec   = [];
ci_lo_v = [];  ci_hi_v = [];

for i = 1:p
    for j = i+1:p
        dc_vec(end+1)  = net.dcor(i,j);  %#ok
        q_vec(end+1)   = net.qval(i,j);  %#ok
        ci_lo_v(end+1) = net.ci_lo(i,j); %#ok
        ci_hi_v(end+1) = net.ci_hi(i,j); %#ok
    end
end

logq    = -log10(max(q_vec, 1e-30));
sig_fdr = q_vec < pm.alpha;
sig_ci  = ci_lo_v > pm.ci_min;

col_map = [0.82 0.20 0.20;
           0.05 0.40 0.75;
           0.85 0.60 0.10;
           0.10 0.62 0.22];
state = 1 + double(sig_fdr) + 2*double(sig_ci);

ax = subplot(1, 1, 1);
hold(ax,'on');

for st = 1:4
    mask = state == st;
    if ~any(mask), continue; end
    for k = find(mask)
        plot(ax, [dc_vec(k) dc_vec(k)], [ci_lo_v(k) ci_hi_v(k)], ...
            '-', 'Color',[0.65 0.65 0.65 0.40], 'LineWidth',0.6, ...
            'HandleVisibility','off');
    end
    scatter(ax, dc_vec(mask), logq(mask), 20, col_map(st,:), 'filled', ...
        'MarkerFaceAlpha',0.65, 'MarkerEdgeColor','none');
end

xline(ax, pm.ci_min, '--', 'Color',[0.40 0.40 0.40], 'LineWidth',1.2, ...
    'Label', sprintf('ci_min=%.2f', pm.ci_min), ...
    'LabelHorizontalAlignment','left','FontSize',7.5);
yline(ax, -log10(pm.alpha), '--', 'Color',[0.70 0.12 0.12], 'LineWidth',1.2, ...
    'Label', sprintf('q=%.2f', pm.alpha), ...
    'LabelHorizontalAlignment','right','FontSize',7.5);

xlabel(ax,'dCor (U-stat)','FontSize',10,'FontWeight','bold');
ylabel(ax,'-log_{10}(q)','FontSize',10,'FontWeight','bold');
title(ax,'Set 1','FontSize',10,'FontWeight','bold','Color',[0.18 0.45 0.78]*0.75);
legend(ax, {'Neither sig','FDR only','CI only','Both sig'}, ...
    'Location','northeast','FontSize',7.5,'Box','on');
xlim(ax,[-0.05 1.05]); grid(ax,'on'); box(ax,'on');
ax.GridAlpha = 0.13;
hold(ax,'off');

annotation(fig,'textbox',[0.01 0.002 0.98 0.025], ...
    'String',['Each point = one edge (feature pair)  |  ' ...
              'Vertical bars = Bernstein CI  |  ' ...
              'Green = FDR AND CI both significant  |  ' ...
              'Blue = FDR only  |  Amber = CI only  |  Red = neither'], ...
    'FontSize',7.5,'EdgeColor','none','HorizontalAlignment','center');
set(findall(fig,'Type','axes'),'FontName','Helvetica');
end


% =========================================================================
%  FIGURE N — NETWORK GRAPHS
% =========================================================================
function fig = plot_networks(net1,  thresh_types, thresh_labels, ...
        fig_colors, pm)

n_thresh  = numel(thresh_types);
sets      = {net1};
set_names = pm.Filename;

fig = figure('Name','Network Graphs','NumberTitle','off', ...
    'Position',[40 40 1800 700],'Color',[0.97 0.97 0.97]);
sgtitle('Pairwise dCor Networks  |  Edge colour = dCor magnitude  |  Node size proportional to degree', ...
    'FontSize',11,'FontWeight','bold');

for s = 1:1
    net = sets{s};
    p   = net.n_features;
    lbl = net.labels;
% only plot selection by qval and CI
    for t = n_thresh:n_thresh
        sp_idx = (s-1)*n_thresh + t;
        ax = subplot(1,1,1); % n_thresh, sp_idx);
        hold(ax,'on');

        adj_field = ['adj_' thresh_types{t}];
        adj = net.(adj_field);
        degall = sum(adj, 2);

        % Build graph
        [ri, ci_idx] = find(triu(adj,1));
        weights = net.dcor(sub2ind([p,p], ri, ci_idx));

        if isempty(ri)
            text(ax, 0.5, 0.5, 'No edges', 'Units','normalized', ...
                'HorizontalAlignment','center','FontSize',12,'Color',[0.5 0.5 0.5]);
            axis(ax,'off');
            title(ax, sprintf('%s  %s', set_names{s}, thresh_labels{t}), ...
                'FontSize',8,'FontWeight','bold');
            continue;
        end

        k = pm.select_label
        Gall = graph(ri, ci_idx, weights, p);
        nbrs = neighbors(Gall, k);
        subnodes = [k; nbrs];

        G = subgraph(Gall, subnodes);

        adjsub = adj(subnodes, subnodes);
        deg    = degall(subnodes);

        % Node positions
        switch pm.layout
            case 'circle'
                theta = linspace(0, 2*pi*(1-1/numel(subnodes)), numel(subnodes))';
                pos   = [cos(theta), sin(theta)];
            otherwise   % force layout
                pos = fruchterman_reingold(adjsub, numel(subnodes));
        end

        % Draw edges coloured by dCor weight
        edge_w = G.Edges.Weight;
        cmap_e = parula(64);
        for e = 1:numedges(G)
            nd  = G.Edges.EndNodes(e,:);
            w   = edge_w(e);
            cidx = max(1, min(64, round(w*63)+1));
            ecol = cmap_e(cidx,:);
            plot(ax, pos(nd,1), pos(nd,2), '-', 'Color',[ecol,0.65], ...
                'LineWidth', 0.5 + w*2.5);
        end

        % Draw nodes sized by degree
        max_deg = max(max(deg),1);
        for k = 1:numel(subnodes)
            mk_sz = 40 + 100*(deg(k)/max_deg);
            scatter(ax, pos(k,1), pos(k,2), mk_sz, fig_colors.(thresh_types{t}), ...
                'filled','MarkerEdgeColor',[1 1 1],'LineWidth',0.8);
            text(ax, pos(k,1), pos(k,2) + 0.08, lbl{subnodes(k)}, ...
                'FontSize',6.5,'HorizontalAlignment','center', ...
                'Color',[0.15 0.15 0.15]);
        end

        xlim(ax,[-1.3 1.3]); ylim(ax,[-1.3 1.3]);
        axis(ax,'equal','off');
        title(ax, sprintf('%s\n%s', set_names{s}, thresh_labels{t}), ...
            'FontSize',8,'FontWeight','bold', ...
            'Color',fig_colors.(thresh_types{t})*0.75);

        % Mini colourbar annotation
cb_ax = axes(fig,'Position',ax.Position.*[1 1 0.02 0.5] + ...
    [ax.Position(1)+ax.Position(3)*0.62, ax.Position(2)+ax.Position(4)*0.25, 0, 0]);
        n_g = 32;
image(cb_ax, 1, linspace(0,1,n_g), permute(parula(n_g),[1 3 2]));
dmin = 0; %min(edge_w);
dmax = 1; %max(edge_w);
tick_vals = linspace(dmin, dmax, 5);
%linspace(0.5, n_g-0.5, 5)
set(cb_ax, 'XTick', [], ...
           'YTick', tick_vals, ...
           'YTickLabel', arrayfun(@(v) sprintf('%.2f',v), tick_vals, 'UniformOutput', false), ...
           'FontSize', 6, 'YDir', 'normal');
h=ylabel(cb_ax, 'dCor', 'FontSize', 7.5, 'Rotation', 0, 'VerticalAlignment', 'bottom');
% Shift label upward
pos = h.Position;
pos(2) = pos(2) + 0.52;   % adjust this value as needed
h.Position = pos;

    end
end
set(findall(fig,'Type','axes'),'FontName','Helvetica');
end


% =========================================================================
%  SAVE_FIGURE
% =========================================================================
function save_figure(fig, out_dir, base_name, net1)

if isempty(fig) || ~ishandle(fig)
    warning('save_figure: invalid figure handle for %s — skipping.', base_name);
    return;
end

png_path = fullfile(out_dir, [base_name '.png']);
try
    exportgraphics(fig, png_path, 'Resolution', 300, 'BackgroundColor', 'white');
    fprintf('    [PNG] %s\n', png_path);
catch ME
    print(fig, png_path, '-dpng', '-r300');
    fprintf('    [PNG fallback] %s  (%s)\n', png_path, ME.message);
end

fig_path = fullfile(out_dir, [base_name '.fig']);
try
    savefig(fig, fig_path);
    fprintf('    [FIG] %s\n', fig_path);
catch ME
    warning('save_figure: could not save .fig for %s: %s', base_name, ME.message);
end

json_path = fullfile(out_dir, [base_name '.json']);
try
    mat2rows = @(M) num2cell(M, 2);
    s = struct();
    s.figure_name    = base_name;
    s.timestamp      = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
    s.n_features     = net1.n_features;
    s.labels         = net1.labels;
    s.net1.n             = net1.n;
    s.net1.alpha         = net1.alpha;
    s.net1.ci_min        = net1.ci_min;
    s.net1.degree_full   = net1.degree_full(:)';
    s.net1.degree_fdr    = net1.degree_fdr(:)';
    s.net1.degree_ci     = net1.degree_ci(:)';
    s.net1.degree_fdr_ci = net1.degree_fdr_ci(:)';
    s.net1.bc_full       = net1.bc_full(:)';
    s.net1.bc_fdr        = net1.bc_fdr(:)';
    s.net1.bc_ci         = net1.bc_ci(:)';
    s.net1.bc_fdr_ci     = net1.bc_fdr_ci(:)';
    s.net1.dcor          = mat2rows(net1.dcor);
    s.net1.pval          = mat2rows(net1.pval);
    s.net1.qval          = mat2rows(net1.qval);
    s.net1.ci_lo         = mat2rows(net1.ci_lo);
    s.net1.ci_hi         = mat2rows(net1.ci_hi);
    s.net1.adj_full      = mat2rows(net1.adj_full);
    s.net1.adj_fdr       = mat2rows(net1.adj_fdr);
    s.net1.adj_ci        = mat2rows(net1.adj_ci);
    s.net1.adj_fdr_ci    = mat2rows(net1.adj_fdr_ci);

    json_str = jsonencode(s, 'PrettyPrint', true);
    fid = fopen(json_path, 'w');
    if fid == -1, error('Cannot open file: %s', json_path); end
    fprintf(fid, '%s', json_str);
    fclose(fid);
    fprintf('    [JSON] %s\n', json_path);
catch ME
    warning('save_figure: JSON export failed for %s: %s', base_name, ME.message);
end
end


% =========================================================================
%  STATISTICS HELPERS
% =========================================================================

function dc = dcor_from_ucentered(Au, Bu, n)
if n < 5, dc = NaN; return; end
mask     = ~eye(n);
dF       = n*(n-3);
dcov2_XY = sum(sum(Au.*Bu.*mask)) / dF;
dcov2_XX = sum(sum(Au.^2.*mask))  / dF;
dcov2_YY = sum(sum(Bu.^2.*mask))  / dF;
denom    = sqrt(dcov2_XX * dcov2_YY);
if denom <= 0, dc = 0; return; end
val = dcov2_XY / denom;
dc  = max(0, sign(val)*sqrt(abs(val)));
end

function Au = u_center(A, n)
rs = sum(A,2); gs = sum(rs);
Au = A - rs*ones(1,n)/(n-2) - ones(n,1)*rs'/(n-2) + gs/((n-1)*(n-2));
Au(1:n+1:end) = 0;
end

function pval = dcor_chi2_test(dcor_val, n)
T    = n * dcor_val^2;
pval = 1 - chi2cdf(T+1, 1);
end

function ci = bernstein_ci_dcor(dc, n, alpha)
if n < 2 || isnan(dc), ci = [NaN NaN]; return; end
K      = n;
lnterm = log(2 / alpha);
eps_v  = (2/3) * lnterm / K;
ci     = max(-1, min(1, [dc-eps_v, dc+eps_v]));
end

function bc = betweenness_centrality(adj)
n  = size(adj, 1);
bc = zeros(n, 1);
for s = 1:n
    dist    = -ones(n,1);  dist(s)  = 0;
    sigma   = zeros(n,1);  sigma(s) = 1;
    pred    = cell(n,1);
    queue   = s;
    delta   = zeros(n,1);

    head = 1;
    while head <= numel(queue)
        v    = queue(head);  head = head+1;
        nbrs = find(adj(v,:));
        for w = nbrs
            if dist(w) < 0
                dist(w) = dist(v)+1;
                queue(end+1) = w; %#ok
            end
            if dist(w) == dist(v)+1
                sigma(w) = sigma(w) + sigma(v);
                pred{w}(end+1) = v;
            end
        end
    end

    for idx = numel(queue):-1:2
        w = queue(idx);
        for v = pred{w}
            delta(v) = delta(v) + sigma(v)/sigma(w) * (1 + delta(w));
        end
        if w ~= s
            bc(w) = bc(w) + delta(w);
        end
    end
end
if n > 2
    bc = bc / ((n-1)*(n-2));
end
end

function pos = fruchterman_reingold(adj, n)
if n == 0, pos = zeros(0,2); return; end
if n == 1, pos = [0 0];      return; end

theta = linspace(0, 2*pi*(1-1/n), n)';
pos   = [cos(theta), sin(theta)] + 0.05*randn(n,2);

area  = 4;
k     = sqrt(area/n);
iters = 80;
temp  = 1.0;
cool  = temp / iters;

for iter = 1:iters
    disp_v = zeros(n,2);
    for v = 1:n
        for u = 1:n
            if u == v, continue; end
            d  = pos(v,:) - pos(u,:);
            dl = max(norm(d), 1e-6);
            disp_v(v,:) = disp_v(v,:) + (d/dl) * k^2/dl;
        end
    end
    for v = 1:n
        nbrs = find(adj(v,:));
        for u = nbrs
            if u <= v, continue; end
            d  = pos(v,:) - pos(u,:);
            dl = max(norm(d), 1e-6);
            f  = dl^2/k;
            disp_v(v,:) = disp_v(v,:) - (d/dl)*f;
            disp_v(u,:) = disp_v(u,:) + (d/dl)*f;
        end
    end
    for v = 1:n
        dl = max(norm(disp_v(v,:)), 1e-6);
        pos(v,:) = pos(v,:) + (disp_v(v,:)/dl) * min(dl, temp);
    end
    temp = temp - cool;
end

pos(:,1) = 2*(pos(:,1)-min(pos(:,1)))/(max(pos(:,1))-min(pos(:,1))+1e-9) - 1;
pos(:,2) = 2*(pos(:,2)-min(pos(:,2)))/(max(pos(:,2))-min(pos(:,2))+1e-9) - 1;
end


% =========================================================================
%  FDR_BH
% =========================================================================
function [h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(pvals, q, method, report)
if nargin < 2, q      = .05;    end
if nargin < 3, method = 'pdep'; end
if nargin < 4, report = 'no';   end

if ~isempty(find(pvals<0,1)), error('Some p-values are less than 0.');   end
if ~isempty(find(pvals>1,1)), error('Some p-values are greater than 1.'); end

s = size(pvals);
if (length(s)>2) || s(1)>1
    [p_sorted, sort_ids] = sort(reshape(pvals,1,prod(s)));
else
    [p_sorted, sort_ids] = sort(pvals);
end
[~, unsort_ids] = sort(sort_ids);
m = length(p_sorted);

if strcmpi(method,'pdep')
    thresh = (1:m)*q/m;
    wtd_p  = m*p_sorted./(1:m);
elseif strcmpi(method,'dep')
    denom  = m*sum(1./(1:m));
    thresh = (1:m)*q/denom;
    wtd_p  = denom*p_sorted./(1:m);
else
    error('Argument ''method'' needs to be ''pdep'' or ''dep''.');
end

if nargout > 3
    adj_p = zeros(1,m)*NaN;
    [wtd_p_sorted, wtd_p_sindex] = sort(wtd_p);
    nextfill = 1;
    for k = 1:m
        if wtd_p_sindex(k) >= nextfill
            adj_p(nextfill:wtd_p_sindex(k)) = wtd_p_sorted(k);
            nextfill = wtd_p_sindex(k)+1;
            if nextfill > m, break; end
        end
    end
    adj_p = reshape(adj_p(unsort_ids), s);
end

rej    = p_sorted <= thresh;
max_id = find(rej, 1, 'last');
if isempty(max_id)
    crit_p      = 0;
    h           = pvals*0;
    adj_ci_cvrg = NaN;
else
    crit_p      = p_sorted(max_id);
    h           = pvals <= crit_p;
    adj_ci_cvrg = 1 - thresh(max_id);
end

if strcmpi(report,'yes')
    n_sig = sum(p_sorted <= crit_p);
    if n_sig == 1
        fprintf('Out of %d tests, %d is significant using a false discovery rate of %f.\n',  m, n_sig, q);
    else
        fprintf('Out of %d tests, %d are significant using a false discovery rate of %f.\n', m, n_sig, q);
    end
end
end