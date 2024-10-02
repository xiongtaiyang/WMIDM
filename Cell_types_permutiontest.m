function cell_gene_permutation_test(basic_cell_genes, combined_CellGenes, csv_data, num_iterations)
% Function to perform cell gene permutation test on given data.
% Inputs:
%   basic_cell_genes - Table containing gene data for basic cell types
%   combined_CellGenes - Combined gene data for different cell types
%   csv_data - Cell array containing csv files data (each element is a table with gene data)
%   num_iterations - Number of iterations for the permutation test

    % Get unique cell classes and remove unwanted classes
    class_name = unique(basic_cell_genes.Class);
    class_name([4,9]) = [];
    combined_CellGenes([4,9],:) = [];

    % Iterate over each CSV data provided
    for i = 1:numel(csv_data)
        % Task name extraction
        task_name = csv_data{i}.name(end-6:end-4);
        
        % Initialize count table for statistics
        count_table = cell2table(class_name);
        count_table.('count') = zeros(length(class_name),1);
        
        % Extract SYMBOL column from the current CSV
        current_var1 = csv_data{i}.SYMBOL;
        
        % For each gene in the current CSV, check against basic cell genes
        for j = 1:size(current_var1,1)
            for m = 1:size(combined_CellGenes,1)
                if sum(strcmpi(current_var1{j,1}, combined_CellGenes.CombinedData{m,1})) ~= 0
                    current_var1{j,m+1} = combined_CellGenes.Class{m,1};
                    count_table.count(m) = count_table.count(m) + 1;
                end
            end
        end
        
        % Perform permutation test (5000 iterations as default)
        num_genes = sum(count_table.count);
        num_cell_types = size(class_name,1);
        [FDRq] = analyzeGeneCellTypes(num_genes, num_cell_types, num_iterations, count_table);
        count_table.pvalue = FDRq';
        
        % Save permutation test results (Replace this with actual variable usage as required)
        disp(['Permutation test results for task: ', task_name]);
        disp(count_table);

        % Save gene counts and perform enrichment analysis (Placeholders, modify as needed)
        disp(['Gene counts and enrichment results for task: ', task_name]);
        % Use enrichment analysis function (you would need to modify this as per your actual analysis)
        extract_GenesEachCellType(current_var1, class_name, task_name);
    end
end
