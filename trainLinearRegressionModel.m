function [trainedModel, validationRMSE, table_fill] = trainRegressionModel_FS(trainingData, responseData, KFolds)
% Function to train a linear regression model with K-Fold cross-validation.
% Inputs:
%   trainingData - Matrix of predictor variables (features).
%   responseData - Vector of response variable (target).
%   KFolds - Number of folds for cross-validation.
%
% Outputs:
%   trainedModel - A structure containing the trained regression model and information.
%   validationRMSE - Validation Root Mean Squared Error (RMSE) from cross-validation.
%   table_fill - A table containing the actual and predicted responses for validation data.

    % Extract predictors and response
    % Convert input matrix to a table with named columns.
    for i = 1:size(trainingData,2)
        var_name{i} =  ['column_' num2str(i)];  % Create column names for the predictors
        ispredictor(i) = false;  % Initialize categorical predictor flag (no categorical variables)
    end
    inputTable = array2table(trainingData, 'VariableNames', var_name);  % Create table for predictors

    predictorNames = var_name;  % List of predictor variable names
    predictors = inputTable(:, predictorNames);  % Extract predictor columns from table
    response = responseData;  % Define the response variable
    isCategoricalPredictor = ispredictor;  % Flag for categorical predictors (none in this case)

    % Train regression model
    % Combine predictors and response for model fitting.
    concatenatedPredictorsAndResponse = predictors;
    concatenatedPredictorsAndResponse.dp = response;  % Add the response column to the table
    linearModel = fitlm(...
        concatenatedPredictorsAndResponse, ...
        'linear', ...
        'RobustOpts', 'off');  % Train linear model without robust options

    % Create function for predicting using the trained model
    predictorExtractionFcn = @(x) array2table(x, 'VariableNames', predictorNames);
    linearModelPredictFcn = @(x) predict(linearModel, x);  % Define prediction function
    trainedModel.predictFcn = @(x) linearModelPredictFcn(predictorExtractionFcn(x));  % Store prediction function in the model structure

    % Add model details to the result structure
    trainedModel.LinearModel = linearModel;
    trainedModel.About = 'This structure contains a trained regression model exported from the Regression Learner app.';
    trainedModel.HowToPredict = sprintf(['To make predictions on new predictor data X, use: \n yfit = c.predictFcn(X) \n' ...
        'where ''c'' is the name of the model structure (e.g., ''trainedModel'').\n \nX must contain exactly 30 columns '...
        'because this model was trained with 30 predictors. X must only contain the predictor columns in the same order and format as the training data.']);

    % Cross-validation
    cvp = cvpartition(size(response, 1), 'KFold', KFolds);  % Create partition for K-Fold cross-validation
    validationPredictions = response;  % Initialize validation predictions

    % Perform K-Fold cross-validation
    for fold = 1:KFolds
        % Get training predictors and response for the current fold
        trainingPredictors = predictors(cvp.training(fold), :);
        trainingResponse = response(cvp.training(fold), :);
        foldIsCategoricalPredictor = isCategoricalPredictor;

        % Train the regression model on the current fold
        concatenatedPredictorsAndResponse = trainingPredictors;
        concatenatedPredictorsAndResponse.dp = trainingResponse;
        linearModel = fitlm(...
            concatenatedPredictorsAndResponse, ...
            'linear', ...
            'RobustOpts', 'off');  % Train linear model without robust options

        % Define prediction function for the current fold
        linearModelPredictFcn = @(x) predict(linearModel, x);
        validationPredictFcn = @(x) linearModelPredictFcn(x);

        % Predict on the validation set for the current fold
        validationPredictors = predictors(cvp.test(fold), :);
        foldPredictions = validationPredictFcn(validationPredictors);

        % Store predictions in the original response order
        validationPredictions(cvp.test(fold), :) = foldPredictions;
    end

    % Calculate Pearson correlation between actual and predicted responses
    [r, p] = corr(response, validationPredictions, 'Type', 'Pearson');
    
    % Fill a table with actual and predicted responses for validation
    table_fill = [response, validationPredictions];

    % Calculate validation RMSE (Root Mean Squared Error)
    isNotMissing = ~isnan(validationPredictions) & ~isnan(response);  % Check for missing values
    validationRMSE = sqrt(nansum((validationPredictions - response).^2) / numel(response(isNotMissing)));  % Calculate RMSE

end
