function D = FractalDim(y, cellmax)
% This function computes the box-counting fractal dimension of a 1D signal.
%
% Inputs:
%   y       - A 1D input signal.
%   cellmax - The maximum size of the boxes (must be a power of 2 and an even number).
%             It should be larger than the length of the input signal y.
%
% Output:
%   D       - The box-counting fractal dimension (generally D >= 1). 
%             D = lim(log(N(e))/log(k/e)), where N(e) is the number of boxes 
%             and k/e is the box size.
%
% First, check if cellmax is greater than or equal to the length of the input signal.
% If not, display an error message.
if cellmax < length(y)
    error('cellmax must be larger than the input signal length!')
end

L = length(y); % Number of samples in the input signal
y_min = min(y); % Minimum value of the signal

% Shift the signal so that its minimum value is 0
y_shift = y - y_min;

% Resample the signal to have cellmax + 1 points
% interp1 performs 1D linear interpolation
x_ord = [0:L-1] / (L-1); % Normalized x-axis for original signal
xx_ord = [0:cellmax] / cellmax; % Normalized x-axis for resampling
y_interp = interp1(x_ord, y_shift, xx_ord); % Interpolated signal

% Scale the resampled signal such that the maximum value is proportional to cellmax
ys_max = max(y_interp); 
factory = cellmax / ys_max; 
yy = abs(y_interp * factory); % Scaled signal

% Calculate the number of iterations based on the logarithm of cellmax
t = log2(cellmax) + 1; % Number of iterations

% Initialize an array to store the number of segments for each iteration
for e = 1:t
    Ne = 0; % Accumulated number of boxes that cover the signal in this iteration
    cellsize = 2^(e-1); % The size of each box in this iteration
    NumSeg(e) = cellmax / cellsize; % Number of segments along the x-axis

    % Iterate through each segment along the x-axis
    for j = 1:NumSeg(e)
        % Define the segment's start and end points
        begin = cellsize * (j-1) + 1;
        tail = min(cellsize * j + 1, cellmax + 1); % Ensure the tail does not exceed bounds
        seg = begin:tail; % The segment's coordinates

        % Find the max and min values of the signal in this segment
        yy_max = max(yy(seg)); 
        yy_min = min(yy(seg));

        % Compute the number of boxes (grid cells) occupied by this segment
        up = ceil(yy_max / cellsize); % Round up to the nearest integer
        down = floor(yy_min / cellsize); % Round down to the nearest integer
        Ns = up - down; % The number of boxes occupied by this segment
        Ne = Ne + Ns; % Accumulate the number of boxes for this iteration
    end

    N(e) = Ne; % Store the total number of boxes for this iteration
end

% Perform linear regression (least squares) on log(N(e)) vs log(k/e)
% The slope of the fitted line gives the fractal dimension D
r = -diff(log2(N)); % Compute differences between consecutive log(N)
id = find(r <= 2 & r >= 1); % Select valid data points where r is between 1 and 2

% Use the selected points for linear regression
Ne = N(id); 
e = NumSeg(id); 
P = polyfit(log2(e), log2(Ne), 1); % Fit a line and return the slope and intercept
D = P(1); % The slope of the fitted line is the fractal dimension
end
