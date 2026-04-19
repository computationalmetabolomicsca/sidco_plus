% =========================================================================================================================
% SINGLE SAMPLE ANALYSIS OF UNBIASED DISTANCE CORRELATION
% MCuperlovic-Culf,Ottawa 2026
%==========================================================================================================================
% User enters:
% - csv file with features (columns)and samples (rows)
% - directory for output files
% - alpha values for stat significance defauls 0.05
% - min CI value for selection of significant size - degault 0.1
% - Normalization option - 1 yes, 0 no
% - Imputation option - 1 yes, 0 no
% - Number of samples for bootstrap subset size
% - Number of bootstrap runs

% Run example
%alpha=0.01;
%minCI=0.1;
%Normalization_box=1;
%Imputationa_box=1;  (1=yes,0=no)
%Nsubsample=50;
%Npermutations=20;
%[net1] =...
%dcor_network_single_sample_first_analysis_permutations("Inut_File.csv","outputdirectory", alpha, minCI, Normalization_box, Imputationa_box,Nsubsample,Npermutations);



function net1 = dcor_network_single_sample_first_analysis_permutations(Filename,outputdirectory,...
    alpha, ci_min, Normalize, Imputer,NSubset,NBoot)

warning('OFF')
opts = detectImportOptions(Filename, 'ReadVariableNames', true, 'PreserveVariableNames', true);
opts.DataLines = [2, Inf];
opts.Delimiter = ",";
inputfile = readtable(Filename, opts);
clear opts


% -------------------------------------------------------------------------
% Defaults
% -------------------------------------------------------------------------

    out_dir = outputdirectory;
end



do_plot   = true;  
layout    = 'force'; 
 seed      = 1;  
 B_perm    = 100; 

rng(seed);

% Pack into struct
if ischar(alpha) || isstring(alpha)
    pm.alpha = str2double(alpha);
else
    pm.alpha = alpha;
end
if ischar(ci_min) || isstring(ci_min)
    pm.ci_min = str2double(ci_min);
else
    pm.ci_min = ci_min;
end
pm.Normalize = Normalize;
pm.Imputer   = Imputer;
pm.plot      = do_plot;
pm.layout    = layout;
pm.seed      = seed;
pm.B_perm    = B_perm;
pm.Filename  =Filename;
if ischar(NSubset) || isstring(NSubset)
    pm.NSubset = str2double(NSubset);
else
    pm.NSubset = NSubset;
end
if ischar(NBoot) || isstring(NBoot)
    pm.NBoot = str2double(NBoot);
else
    pm.NBoot = NBoot;
end

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

fig1 = plot_heatmaps(net1,pm, out_dir, 'fig1_dcor_heatmaps');
%save_figure(fig1, out_dir, 'fig1_dcor_heatmaps', net1);

fig3 = plot_degree_comparison(net1, threshold_types, threshold_labels, fig_colors,pm, out_dir, 'fig3_degree_comparison');
%save_figure(fig3, out_dir, 'fig3_degree_comparison', net1);

fig4 = plot_bc_comparison(net1, threshold_types, threshold_labels, fig_colors,pm, out_dir, 'fig4_betweeenness_comparison');
%save_figure(fig4, out_dir, 'fig4_betweenness_centrality', net1);


%fig5 = plot_degree_comparisonF(net1, fig_colors,pm, out_dir, 'fig5_degree_full');
%save_figure(fig3, out_dir, 'fig3_degree_comparison', net1);

%fig6 = plot_bc_comparisonF(net1, threshold_types, threshold_labels, fig_colors,pm, out_dir, 'fig6_betweeenness_full');
%save_figure(fig4, out_dir, 'fig4_betweenness_centrality', net1);

%fig5 = plot_edge_summary(net1, pm);
%save_figure(fig5, out_dir, 'fig5_edge_summary', net1);

fprintf('  All figures saved.\n');
end


% =========================================================================
%  FIGURE 1 — HEATMAPS
% =========================================================================
function fig = plot_heatmaps(net1,pm, out_dir, base_name)

Filename=pm.Filename;
% Create figure without displaying it with ,'Visible', 'off');

fig = figure('Name','dCor Heatmaps','NumberTitle','off', ...
    'Position',[30 30 1600 780],'Color',[0.97 0.97 0.97],'Visible', 'off');
sgtitle('Pairwise dCor (U-stat)  |  p-value  |  Bernstein CI width', ...
    'FontSize',12,'FontWeight','bold');

net = net1;
p   = net.n_features;
lbl = net.labels;

ax1 = subplot(1, 2, 1);
imagesc(ax1, net.dcor, 'AlphaData', ~isnan(net.dcor));
set(ax1, 'Color', 'white');
colorbar(ax1); clim(ax1,[0 1]);
colormap(ax1, parula(256));
set(ax1,'XTick',1:p,'XTickLabel',lbl,'YTick',1:p,'YTickLabel',lbl, ...
    'XTickLabelRotation',90,'FontSize',6,'TickLength',[0 0]);
title(ax1,[Filename, ' —  dCor'],'FontSize',9,'FontWeight','bold');
axis(ax1,'square');

tempdcor = net.dcor; tempdcor(tempdcor==0 | net.adj_fdr_ci==0) = nan;
ax4 = subplot(1, 2, 2);
imagesc(ax4, tempdcor, 'AlphaData', ~isnan(tempdcor));
set(ax4, 'Color', 'white');
colorbar(ax4); clim(ax4,[0 1]);
colormap(ax4, parula(256));
set(ax4,'XTick',1:p,'XTickLabel',lbl,'YTick',1:p,'YTickLabel',lbl, ...
    'XTickLabelRotation',90,'FontSize',6,'TickLength',[0 0]);
title(ax4,[Filename, ' —  dCor FDR+ Bernstein CI pass'],'FontSize',9,'FontWeight','bold');
axis(ax4,'square');

annotation(fig,'textbox',[0.01 0.002 0.98 0.025], ...
    'String',['Diagonal zeroed  |  q-values: Benjamini-Hochberg FDR  |  ' ...
              'CI: Bernstein (Maurer & Pontil 2009)'], ...
    'FontSize',7.5,'EdgeColor','none','HorizontalAlignment','center');
set(findall(fig,'Type','axes'),'FontName','Helvetica');

% =========================================================================
%  SAVE_FIGURE
% =========================================================================


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
    s.net1.dcor          = mat2rows(net1.dcor);           % subplot 1: full dCor
    % --- subplot 2: dcor filtered to only FDR + Bernstein CI passing edges ---
    tempdcor_json = net1.dcor;
    tempdcor_json(tempdcor_json == 0 | net1.adj_fdr_ci == 0) = nan;
    s.net1.dcor_fdr_ci   = mat2rows(tempdcor_json);       % subplot 2: filtered 

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
%  FIGURE 3 — DEGREE COMPARISON
% =========================================================================
function fig = plot_degree_comparison(net1, thresh_types, thresh_labels, fig_colors,pm, out_dir, base_name)

Filename=pm.Filename;
n_thresh = numel(thresh_types);
n_shared = numel(net1.labels);
x_pos    = 1:n_shared;

fig = figure('Name','Degree Information','NumberTitle','off', ...
    'Position',[50 50 1600 900],'Color',[0.97 0.97 0.97],'Visible', 'off');
sgtitle(['Node Degree: ',Filename] , 'FontSize',11,'FontWeight','bold');

% only show original and corrected by both qval and ci
for t = n_thresh: n_thresh
    col = fig_colors.(thresh_types{1});
    d1  = net1.(['degree_' thresh_types{t}]);

    ax1 = subplot(1,1,1); % n_thresh, t);
    hold(ax1,'on');
    bar(ax1, x_pos, d1, 0.35, 'FaceColor', col, 'EdgeColor','none', ...
        'FaceAlpha',0.85, 'DisplayName',Filename);
    hold(ax1,'off');

    ax1.XLim               = [min(x_pos)-0.5, max(x_pos)+0.5];
    ax1.XTickMode          = 'manual';
    ax1.XTickLabelMode     = 'manual';
    ax1.XTick              = x_pos;
    ax1.XTickLabel         = net1.labels;
    ax1.XTickLabelRotation = 90;
    ax1.FontSize           = 6;
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


% =========================================================================
%  SAVE_FIGURE
% =========================================================================
%function save_figure(fig, out_dir, base_name, net1)

%if isempty(fig) || ~ishandle(fig)
%    warning('save_figure: invalid figure handle for %s — skipping.', base_name);
%    return;
%end

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
    t        = n_thresh;                          % same index used in the plot
    deg_field = ['degree_' thresh_types{t}];      % e.g. 'degree_fdr_ci'

    s = struct();
    s.figure_name    = base_name;
    s.timestamp      = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
    s.n_features     = net1.n_features;
    s.labels         = net1.labels;
    s.net1.n             = net1.n;
    s.net1.alpha         = net1.alpha;
    s.net1.ci_min        = net1.ci_min;
%    s.net1.degree_full   = net1.degree_full(:)';
%    s.net1.degree_fdr    = net1.degree_fdr(:)';
%    s.net1.degree_ci     = net1.degree_ci(:)';
    s.net1.degree_fdr_ci = net1.degree_fdr_ci(:)';
%    s.net1.bc_full       = net1.bc_full(:)';
%    s.net1.bc_fdr        = net1.bc_fdr(:)';
%    s.net1.bc_ci         = net1.bc_ci(:)';
%    s.net1.bc_fdr_ci     = net1.bc_fdr_ci(:)';
%    s.net1.dcor          = mat2rows(net1.dcor);
%    s.net1.pval          = mat2rows(net1.pval);
%    s.net1.qval          = mat2rows(net1.qval);
%    s.net1.ci_lo         = mat2rows(net1.ci_lo);
%    s.net1.ci_hi         = mat2rows(net1.ci_hi);
%    s.net1.adj_full      = mat2rows(net1.adj_full);
%    s.net1.adj_fdr       = mat2rows(net1.adj_fdr);
%    s.net1.adj_ci        = mat2rows(net1.adj_ci);
%    s.net1.adj_fdr_ci    = mat2rows(net1.adj_fdr_ci);

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
%  FIGURE 4 — BETWEENNESS CENTRALITY
% =========================================================================
function fig = plot_bc_comparison(net1, thresh_types, thresh_labels, fig_colors,pm, out_dir, base_name)

Filename=pm.Filename;


n_thresh = numel(thresh_types);
n_shared = numel(net1.labels);
x_pos    = 1:n_shared;

fig = figure('Name','Betweenness Centrality','NumberTitle','off', ...
    'Position',[60 60 1600 900],'Color',[0.97 0.97 0.97],'Visible', 'off');
sgtitle(['Betweenness Centrality:', Filename], 'FontSize',11,'FontWeight','bold');
% only show qval and ci corrected
for t = n_thresh:n_thresh
    col = fig_colors.(thresh_types{t});
    b1  = net1.(['bc_' thresh_types{t}]);

    ax1 = subplot(1,1,1); % n_thresh, t);
    hold(ax1,'on');
    bar(ax1, x_pos, b1, 0.35, 'FaceColor', col, 'EdgeColor','none', ...
        'FaceAlpha',0.85, 'DisplayName',Filename );
    hold(ax1,'off');

    ax1.XLim               = [min(x_pos)-0.5, max(x_pos)+0.5];
    ax1.XTickMode          = 'manual';
    ax1.XTickLabelMode     = 'manual';
    ax1.XTick              = x_pos;
    ax1.XTickLabel         = net1.labels;
    ax1.XTickLabelRotation = 90;
    ax1.FontSize           = 6;
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

% =========================================================================
%  SAVE_FIGURE
% =========================================================================
%function save_figure(fig, out_dir, base_name, net1)

%if isempty(fig) || ~ishandle(fig)
%    warning('save_figure: invalid figure handle for %s — skipping.', base_name);
%    return;
%end

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
    t        = n_thresh;                          % same index used in the plot
    deg_field = ['degree_' thresh_types{t}];      % e.g. 'degree_fdr_ci'

    s = struct();
    s.figure_name    = base_name;
    s.timestamp      = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
    s.n_features     = net1.n_features;
    s.labels         = net1.labels;
    s.net1.n             = net1.n;
    s.net1.alpha         = net1.alpha;
    s.net1.ci_min        = net1.ci_min;
%    s.net1.degree_full   = net1.degree_full(:)';
%    s.net1.degree_fdr    = net1.degree_fdr(:)';
%    s.net1.degree_ci     = net1.degree_ci(:)';
%    s.net1.degree_fdr_ci = net1.degree_fdr_ci(:)';
%    s.net1.bc_full       = net1.bc_full(:)';
%    s.net1.bc_fdr        = net1.bc_fdr(:)';
%    s.net1.bc_ci         = net1.bc_ci(:)';
    s.net1.bc_fdr_ci     = net1.bc_fdr_ci(:)';
%    s.net1.dcor          = mat2rows(net1.dcor);
%    s.net1.pval          = mat2rows(net1.pval);
%    s.net1.qval          = mat2rows(net1.qval);
%    s.net1.ci_lo         = mat2rows(net1.ci_lo);
%    s.net1.ci_hi         = mat2rows(net1.ci_hi);
%    s.net1.adj_full      = mat2rows(net1.adj_full);
%    s.net1.adj_fdr       = mat2rows(net1.adj_fdr);
%    s.net1.adj_ci        = mat2rows(net1.adj_ci);
%    s.net1.adj_fdr_ci    = mat2rows(net1.adj_fdr_ci);

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
%  FIGURE 5 — DEGREE COMPARISON
% =========================================================================
function fig = plot_degree_comparisonF(net1,  fig_colors,pm, out_dir, base_name)

Filename=pm.Filename;
%n_thresh = numel(thresh_types);
n_shared = numel(net1.labels);
x_pos    = 1:n_shared;

fig = figure('Name','Degree Information','NumberTitle','off', ...
    'Position',[50 50 1600 900],'Color',[0.97 0.97 0.97],'Visible', 'off');
sgtitle(['Node Degree: ',Filename] , 'FontSize',11,'FontWeight','bold');

% only show original and corrected by both qval and ci
for t = 1: 1
    col = [0 0.5 0]; % dark green
    d1  = net1.degree_full;

    ax1 = subplot(1,1,1); % n_thresh, t);
    hold(ax1,'on');
    bar(ax1, x_pos, d1, 0.35, 'FaceColor', col, 'EdgeColor','none', ...
        'FaceAlpha',0.85, 'DisplayName',Filename);
    hold(ax1,'off');

    ax1.XLim               = [min(x_pos)-0.5, max(x_pos)+0.5];
    ax1.XTickMode          = 'manual';
    ax1.XTickLabelMode     = 'manual';
    ax1.XTick              = x_pos;
    ax1.XTickLabel         = net1.labels;
    ax1.XTickLabelRotation = 90;
    ax1.FontSize           = 6;
    ax1.TickLabelInterpreter = 'none';

    ylabel(ax1, 'Degree', 'FontSize', 9);
    title(ax1, "Degree Weighted", 'FontSize',9, 'FontWeight','bold', 'Color',col*0.75);
    legend(ax1, 'Degree', 'Location','northwest', 'FontSize',7, 'Box','off');
    grid(ax1,'on'); box(ax1,'on'); ax1.GridAlpha = 0.13;
end

annotation(fig,'textbox',[0.01 0.002 0.98 0.025], ...
    'String','Bars = Node Degree', ...
    'FontSize',7.5,'EdgeColor','none','HorizontalAlignment','center');
set(findall(fig,'Type','axes'),'FontName','Helvetica');


% =========================================================================
%  SAVE_FIGURE
% =========================================================================
%function save_figure(fig, out_dir, base_name, net1)

%if isempty(fig) || ~ishandle(fig)
%    warning('save_figure: invalid figure handle for %s — skipping.', base_name);
%    return;
%end

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
    t        = n_thresh;                          % same index used in the plot
    deg_field = ['degree_' thresh_types];      % e.g. 'degree_fdr_ci'

    s = struct();
    s.figure_name    = base_name;
    s.timestamp      = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
    s.n_features     = net1.n_features;
    s.labels         = net1.labels;
    s.net1.n             = net1.n;
    s.net1.alpha         = net1.alpha;
    s.net1.ci_min        = net1.ci_min;
    s.net1.degree_full   = net1.degree_full(:)';
%    s.net1.degree_fdr    = net1.degree_fdr(:)';
%    s.net1.degree_ci     = net1.degree_ci(:)';
%    s.net1.degree_fdr_ci = net1.degree_fdr_ci(:)';
%    s.net1.bc_full       = net1.bc_full(:)';
%    s.net1.bc_fdr        = net1.bc_fdr(:)';
%    s.net1.bc_ci         = net1.bc_ci(:)';
%    s.net1.bc_fdr_ci     = net1.bc_fdr_ci(:)';
%    s.net1.dcor          = mat2rows(net1.dcor);
%    s.net1.pval          = mat2rows(net1.pval);
%    s.net1.qval          = mat2rows(net1.qval);
%    s.net1.ci_lo         = mat2rows(net1.ci_lo);
%    s.net1.ci_hi         = mat2rows(net1.ci_hi);
%    s.net1.adj_full      = mat2rows(net1.adj_full);
%    s.net1.adj_fdr       = mat2rows(net1.adj_fdr);
%    s.net1.adj_ci        = mat2rows(net1.adj_ci);
%    s.net1.adj_fdr_ci    = mat2rows(net1.adj_fdr_ci);

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
%  FIGURE 4 — BETWEENNESS CENTRALITY
% =========================================================================
function fig = plot_bc_comparisonF(net1, thresh_types, thresh_labels, fig_colors,pm, out_dir, base_name)

Filename=pm.Filename;


n_thresh = numel(thresh_types);
n_shared = numel(net1.labels);
x_pos    = 1:n_shared;

fig = figure('Name','Betweenness Centrality','NumberTitle','off', ...
    'Position',[60 60 1600 900],'Color',[0.97 0.97 0.97],'Visible', 'off');
sgtitle(['Betweenness Centrality:', Filename], 'FontSize',11,'FontWeight','bold');
% only show qval and ci corrected
for t = n_thresh:n_thresh
    col = fig_colors.(thresh_types{t});
    b1  = net1.(['bc_' thresh_types{t}]);

    ax1 = subplot(1,1,1); % n_thresh, t);
    hold(ax1,'on');
    bar(ax1, x_pos, b1, 0.35, 'FaceColor', col, 'EdgeColor','none', ...
        'FaceAlpha',0.85, 'DisplayName',Filename );
    hold(ax1,'off');

    ax1.XLim               = [min(x_pos)-0.5, max(x_pos)+0.5];
    ax1.XTickMode          = 'manual';
    ax1.XTickLabelMode     = 'manual';
    ax1.XTick              = x_pos;
    ax1.XTickLabel         = net1.labels;
    ax1.XTickLabelRotation = 90;
    ax1.FontSize           = 6;
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

% =========================================================================
%  SAVE_FIGURE
% =========================================================================
%function save_figure(fig, out_dir, base_name, net1)

%if isempty(fig) || ~ishandle(fig)
%    warning('save_figure: invalid figure handle for %s — skipping.', base_name);
%    return;
%end

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
    t        = n_thresh;                          % same index used in the plot
    deg_field = ['degree_' thresh_types{t}];      % e.g. 'degree_fdr_ci'

    s = struct();
    s.figure_name    = base_name;
    s.timestamp      = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
    s.n_features     = net1.n_features;
    s.labels         = net1.labels;
    s.net1.n             = net1.n;
    s.net1.alpha         = net1.alpha;
    s.net1.ci_min        = net1.ci_min;
%    s.net1.degree_full   = net1.degree_full(:)';
%    s.net1.degree_fdr    = net1.degree_fdr(:)';
%    s.net1.degree_ci     = net1.degree_ci(:)';
%    s.net1.degree_fdr_ci = net1.degree_fdr_ci(:)';
%    s.net1.bc_full       = net1.bc_full(:)';
%    s.net1.bc_fdr        = net1.bc_fdr(:)';
%    s.net1.bc_ci         = net1.bc_ci(:)';
    s.net1.bc_fdr_ci     = net1.bc_fdr_ci(:)';
%    s.net1.dcor          = mat2rows(net1.dcor);
%    s.net1.pval          = mat2rows(net1.pval);
%    s.net1.qval          = mat2rows(net1.qval);
%    s.net1.ci_lo         = mat2rows(net1.ci_lo);
%    s.net1.ci_hi         = mat2rows(net1.ci_hi);
%    s.net1.adj_full      = mat2rows(net1.adj_full);
%    s.net1.adj_fdr       = mat2rows(net1.adj_fdr);
%    s.net1.adj_ci        = mat2rows(net1.adj_ci);
%    s.net1.adj_fdr_ci    = mat2rows(net1.adj_fdr_ci);

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
function [h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(pvals, q)
method = 'pdep'; 
report = 'no';   

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