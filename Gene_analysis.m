% clc
% clear
load('E:\HCP\WM-Getm-over\Results\Gene_expression\Hurst\Hurst_mat\effect_size_Body_Others.mat');
load('E:\HCP\WM-Getm-over\Results\Gene_expression\expression_schaefer400\expression_data.mat');

% Replace -Inf with NaN
no_nan = ~isnan(expression(:,1));

%% 找到相关系数最大的基因
t_gene_corr = corr(T_data(no_nan), expression(no_nan,:));
[max_corr,max_index] = max(t_gene_corr);
[min_corr,min_index] = min(t_gene_corr);
maxcorr_expression = expression(:, max_index);
maxcorr_genename = gene_names(max_index);

mincorr_expression = expression(:, min_index);
mincorr_genename = gene_names(min_index);

maxcorr_expression_noInf = maxcorr_expression;
maxcorr_expression_noInf(isinf(maxcorr_expression_noInf)) = nan;

mincorr_expression_noInf = mincorr_expression;
mincorr_expression_noInf(isinf(mincorr_expression_noInf)) = nan;

% Spin permutation testing for two cortical maps
[spin_p, spin_d] = spin_test(T_data', maxcorr_expression_noInf, ...
    'surface_name', 'fsa5', 'parcellation_name', 'schaefer_400', ...
    'n_rot', 1000, 'type', 'pearson');

% Store p-value and null distribution
p_and_d = cell2struct({[spin_p; spin_d]}, {'wfdc1'}, 2);
 
% Plot null distribution
f = figure;
 
    set(gcf,'color','w');
    set(gcf,'units','normalized','position',[0 0 1 0.5])
    fns = fieldnames(p_and_d);
 
    % Define plot colors
    col = [0.66 0.13 0.11];
 
    % Plot null distributions
    h = histogram(p_and_d.(fns{1})(2:end), 50, 'Normalization', ...
        'pdf', 'edgecolor', 'w', 'facecolor', col, 'facealpha', ...
        1, 'linewidth', 0.5);
    l = line([max_corr max_corr], get(gca, 'ylim'), ...
        'linestyle', '--', ...
             'color', 'k', 'linewidth', 1.5);
    xlabel(['Null correlations' newline '(' strrep(fns{1}, '_', ' ') ')'])
    ylabel('Density')
    legend(l,['{\it r}=' num2str(round(max_corr, 2)) newline ...
              '{\it p}=' num2str(round(p_and_d.(fns{1})(1), 3))])
    legend boxoff
    set(gca, 'box', 'off')

%% 功能注释
data_lh = load_mgh('E:\HCP\WM-Getm-over\Results\Gene_expression\T_contrast\mgh_path\Group_T\lh.2bk_0bk.mgh');
data_rh = load_mgh('E:\HCP\WM-Getm-over\Results\Gene_expression\T_contrast\mgh_path\Group_T\rh.2bk_0bk.mgh');
data = [data_lh;data_rh].';
% 
[correlation, feature] = meta_analytic_decoder(data,'template', 'fsaverage5');
disp(correlation(1:3));
figure();
wc = wordcloud(feature(correlation>0), correlation(correlation>0));