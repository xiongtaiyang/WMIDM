function [Poincare_Index] = Poincare(rr_segement)
% Function to calculate the Poincare Index (SD1/SD2 ratio) based on an RR interval segment
% Inputs:
%   rr_segement - A vector or matrix containing RR intervals (row vector expected)
% Outputs:
%   Poincare_Index - The ratio SD1/SD2 for the Poincare plot analysis

% Get the dimensions of rr_segement
[COL, ROW] = size(rr_segement);

% Define the length of the segment for iteration (number of RR pairs)
T = ROW - 1;

% Initialize temporary variables for the sum of squared differences
temp1 = 0;
temp2 = 0;

% Loop through each RR interval pair and calculate temp1 and temp2
for i = 1:T
    % Calculate the squared difference for SD1 (short-term variability)
    temp1 = temp1 + ((rr_segement(i) - rr_segement(i + 1))^2) / 2;
    
    % Calculate the squared difference for SD2 (long-term variability)
    temp2 = temp2 + ((rr_segement(i) + rr_segement(i + 1) - 2 * mean(rr_segement))^2) / 2;
end

% Calculate SD1 (short-term variability)
SD1 = sqrt(temp1 / T);

% Calculate SD2 (long-term variability)
SD2 = sqrt(temp2 / T);

% Calculate the Poincare Index (SD1/SD2 ratio)
Poincare_Index = SD1 / SD2;

end
