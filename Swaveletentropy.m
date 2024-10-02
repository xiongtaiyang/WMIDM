function SWT = Swaveletentropy(signal, Fs)
% Swaveletentropy - Computes the wavelet entropy of a signal using the 
%                   Synchrosqueezed Wavelet Transform (SWT).
%
% Inputs:
%   signal - The input signal (1D vector)
%   Fs - The sampling frequency (scalar)
%
% Output:
%   SWT - The computed wavelet entropy (scalar)

% Perform synchrosqueezed wavelet transform (SWT)
[wt, f] = wsst(signal, Fs);  % wt is the wavelet coefficients, f is the frequency vector

% Calculate power distribution (energy)
power = abs(wt).^2;  % Power spectrum (squared magnitude of wavelet coefficients)
total_power = sum(power(:));  % Total power in the signal
p = power / total_power;  % Normalize power to get probabilities

% Calculate wavelet entropy
p(p == 0) = eps;   % Replace zero probabilities with a small value (to avoid log(0))
SWT = -sum(p(:) .* log(p(:)));  % Compute the wavelet entropy

% Display the result
disp(['Wavelet Entropy: ', num2str(SWT)]);
end
