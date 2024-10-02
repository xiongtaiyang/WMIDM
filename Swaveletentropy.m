function SWT = Swaveletentropy(signal,Fs)
% 输入信号和采样率

% 执行同步挤压小波变换
[wt, f] = wsst(signal, Fs);

% 计算能量分布
power = abs(wt).^2;
total_power = sum(power(:));
p = power / total_power;

% 计算小波熵
p(p == 0) = eps;   % 避免 log(0)
SWT = -sum(p(:) .* log(p(:)));

% 显示结果
disp(['Wavelet Entropy: ', num2str(SWT)]);
