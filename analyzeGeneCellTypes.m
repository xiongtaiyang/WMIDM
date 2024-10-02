function [FDRq] = analyzeGeneCellTypes(num_genes, num_cell_types, num_iterations, count_table)
    % Function to analyze gene-cell type associations with permutation testing and FDR correction.
    % 
    % Parameters:
    % num_genes: The total number of genes.
    % num_cell_types: The number of cell types to assign genes to.
    % num_iterations: The number of iterations for the resampling process.
    % fdr_threshold_q: The FDR threshold (e.g., 0.05).
    %
    % Returns:
    % observed_counts: Observed counts of genes per cell type.
    % p_values: P-values for each cell type.
    % fdr_threshold: The FDR-adjusted p-value threshold.
    % significant: A logical array indicating which cell types are significant.

    % Placeholder for storing resampling results
    resample_counts = zeros(num_iterations, num_cell_types);

    % Resampling and counting
    for i = 1:num_iterations
        % Randomly assign genes to cell types
        cell_type_assignments = randi(num_cell_types, num_genes, 1);
        for j = 1:num_cell_types
            resample_counts(i, j) = sum(cell_type_assignments == j);
        end
    end

    % Observed counts for each cell type
    observed_counts = count_table.count;


    % Calculate p-values
    p_values = zeros(1, num_cell_types);
    for j = 1:num_cell_types
        p_values(j) = mean(resample_counts(:, j) >= observed_counts(j));
    end

    % FDR correction (Benjamini-Hochberg method)
   FDRq = mafdr(p_values);
end
