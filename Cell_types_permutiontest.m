clc
clear
%% Cell_Gene permutation test
% 定义文件夹路径
folder_path = 'E:\HCP\WM-Getm-over\Results\Gene_expression\Regress\2bk_0bk_map_GeneResult\Filtered_Genes0.2';
output_path = 'E:\HCP\WM-Getm-over\Results\Gene_expression\Regress\2bk_0bk_map_GeneResult\CellEnrichment';
% 读取七种基本细胞基因文件
basicCellGeneFile = 'C:\Users\DELL\Desktop\小论文\gene_analysis';
basic_cell_genes = readtable(fullfile(basicCellGeneFile,'七种基本细胞基因.csv'));
class_name = unique(basic_cell_genes.Class);
class_name([4,9]) = [];
% 获取文件夹下所有的csv文件
csv_files = dir(fullfile(folder_path, 'PLS3*.csv'));

% 将每种Cell类别的基因整合到一起
[combined_CellGenes] = combined_Gene_Cell(basic_cell_genes);
combined_CellGenes([4,9],:) = [];
% 遍历每个csv文件
for i = 1:numel(csv_files)
    % task_name, Neg or Pos
    task_name = csv_files(i).name(end-6:end-4);
    % 读取当前csv文件
    current_csv = readtable(fullfile(folder_path, csv_files(i).name));
    % 初始化统计结果
    count_table = cell2table(class_name);
    count_table.('count') = zeros(length(class_name),1);
    % 获取当前csv文件中Symbol列的内容
    current_var1 = current_csv.SYMBOL;
    for j = 1:size(current_var1,1)
        % 找到与基本细胞基因中Var5~Var713列内容相同的行索引
        %[row,col] = find(strcmpi(current_var1{j,1},basic_cell_genes{:,5:end})==1);
        % if ~isempty(row)
        %     for k = 1:size(row,1)
        %         index = find(strcmpi(basic_cell_genes.Class{row(k)},count_table.class_name)==1);
        %         count_table.count(index) = count_table.count(index) + 1;
        %     end
        % end
        for m = 1:size(combined_CellGenes,1)
            if sum(strcmpi(current_var1{j,1},combined_CellGenes.CombinedData{m,1}))~=0
                current_var1{j,m+1} = combined_CellGenes.Class{m,1};
                count_table.count(m) = count_table.count(m)+1;
            end

        end
    end

    % 5000次置换检验
    num_genes = sum(count_table.count);
    num_cell_types = size(class_name,1);
    num_iterations = 5000;
    [FDRq] = analyzeGeneCellTypes(num_genes, num_cell_types, num_iterations, count_table);
    count_table.pvalue = FDRq';


    % 保存置换检验p值和基因在每个cell类别上的数量
    output_file1 = fullfile(output_path,['CellType_',task_name,'_PLS3.csv']);
    writetable(count_table, output_file1);

    % 保存所有类别的Genes的总体文件
    output_file2 = fullfile(output_path,['CellType_Genes_',task_name,'_PLS3.csv']);
    writetable(cell2table(current_var1), output_file2);
    % 保存每个cell类别的Genes的csv文件，用来做富集分析
    extract_GenesEachCellType(current_var1,class_name,task_name);

end