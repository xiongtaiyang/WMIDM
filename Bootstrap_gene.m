function run_pls_bootstrap(expression_data, gene_names, WMIDM, bootnum)
% Function to perform PLS and bootstrap on given expression data and T-maps.
% Inputs:
%   expression_data - Predictor variables (matrix X)
%   gene_names - List of gene names
%   WMIDM - Structure containing average pearson coefficient based on shcaefer400 atlas(1*400 vector) with fields:
%                   'name', 'mean_r_feature_zscore', and 'folder'
%   bootnum - Number of bootstrap iterations

    % Initialize gene index and remove NaN rows
    genes = gene_names;
    geneindex = 1:size(genes, 2);
    X = expression_data;
    index_nan = isnan(X(:, 1));
    X(index_nan, :) = [];

    % Iterate over T-map data
    for num = 1:length(WMIDM)
        mean_r_feature_zscore = WMIDM(num).mean_r_feature_zscore;
        output_path = fullfile(WMIDM(num).folder, WMIDM(num).name(1:end-4));
        mkdir(output_path);

        if size(mean_r_feature_zscore, 1) == 1
            mean_r_feature_zscore = mean_r_feature_zscore';  % Response variable
        end
        Y = mean_r_feature_zscore;
        Y(index_nan) = [];

        % Perform PLS with 8 components
        dim = 8;
        [XL, YL, XS, YS, BETA, PCTVAR, MSE, stats] = plsregress(X, Y, dim);

        % Calculate correlation and sort components
        [R1, p1] = corr(XS, Y);
        [~, index] = max(R1);

        % Align PLS components for interpretability
        if R1(1, 1) < 0
            stats.W(:, 1) = -1 * stats.W(:, 1);
            XS(:, 1) = -1 * XS(:, 1);
        end
        if R1(index, 1) < 0
            stats.W(:, index) = -1 * stats.W(:, index);
            XS(:, index) = -1 * XS(:, index);
        end

        % Sort genes by PLS weights
        [PLS1w, x1] = sort(stats.W(:, 1), 'descend');
        PLS1ids = genes(x1);
        geneindex1 = geneindex(x1);
        [PLS3w, x2] = sort(stats.W(:, index), 'descend');
        PLS3ids = genes(x2);
        geneindex2 = geneindex(x2);

        % Save ROI scores
        csvwrite([output_path, '\PLS1_ROIscores.csv'], XS(:, 1));
        csvwrite([output_path, '\PLS3_ROIscores.csv'], XS(:, index));

        % Initialize bootstrap weight storage
        PLS1weights = [];
        PLS3weights = [];

        % Bootstrap process
        for i = 1:bootnum
            myresample = randsample(size(X, 1), size(X, 1), 1);
            Xr = X(myresample, :);
            Yr = Y(myresample, :);
            [~, ~, XS, ~, ~, ~, ~, stats] = plsregress(Xr, Yr, dim);

            % Process PLS1 weights
            temp = stats.W(:, 1);
            newW = temp(x1);
            if corr(PLS1w, newW) < 0
                newW = -1 * newW;
            end
            PLS1weights = [PLS1weights, newW];

            % Process PLS3 weights
            temp = stats.W(:, index);
            newW = temp(x2);
            if corr(PLS3w, newW) < 0
                newW = -1 * newW;
            end
            PLS3weights = [PLS3weights, newW];
        end

        % Standard deviation of bootstrap weights
        PLS1sw = std(PLS1weights');
        PLS3sw = std(PLS3weights');

        % Calculate Z-scores for bootstrap weights
        temp1 = PLS1w ./ PLS1sw';
        temp2 = PLS3w ./ PLS3sw';

        % FDR correction for p-values
        temp_p1 = sum((PLS1weights ./ PLS1w > 1), 2) ./ size(PLS1weights, 2);
        temp_p2 = sum((PLS3weights ./ PLS3w > 1), 2) ./ size(PLS3weights, 2);
        p1_fdr = mafdr(temp_p1);
        p2_fdr = mafdr(temp_p2);

        % Sort and save results
        [Z1, ind1] = sort(temp1, 'descend');
        PLS1 = PLS1ids(ind1);
        geneindex1 = geneindex1(ind1);
        [Z2, ind2] = sort(temp2, 'descend');
        PLS3 = PLS3ids(ind2);
        geneindex2 = geneindex2(ind2);

        % Save gene weights
        fid1 = fopen([output_path, '\PLS1_geneWeights.csv'], 'w');
        fprintf(fid1, '%s,%s,%s\n', "Gene_names", "P_adjusted", "Zscore");
        for i = 1:length(genes)
            fprintf(fid1, '%s, %d, %f\n', PLS1{i}, p1_fdr(i), Z1(i));
        end
        fclose(fid1);

        fid2 = fopen([output_path, '\PLS3_geneWeights.csv'], 'w');
        fprintf(fid2, '%s,%s,%s\n', "Gene_names", "P_adjusted", "Zscore");
        for i = 1:length(genes)
            fprintf(fid2, '%s, %d, %f\n', PLS3{i}, p2_fdr(i), Z2(i));
        end
        fclose(fid2);
    end
end
