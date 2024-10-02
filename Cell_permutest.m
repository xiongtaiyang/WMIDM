clc
clear
% 原始数据的数量
original_counts = [50, 60, 70, 80, 90, 100, 110];

% 新生成的数据的数量
new_counts = [20, 30, 40, 50, 60, 70, 90];

% 设置置换检验的参数
num_permutations = 10000;
p_values = zeros(1, length(original_counts));

% 生成数据
rng('default'); % 为了结果可重复设置随机种子
original_data = arrayfun(@(n) randn(n, 1), original_counts, 'UniformOutput', false);
new_data = arrayfun(@(n) randn(n, 1), new_counts, 'UniformOutput', false);

% 置换检验
for i = 1:length(original_counts)
    % 合并两个数据集
    combined_data = [original_data{i}; new_data{i}];
    n1 = length(original_data{i});
    n2 = length(new_data{i});
    
    % 计算原始数据的均值差
    orig_diff = mean(original_data{i}) - mean(new_data{i});
    
    % 置换
    perm_diffs = zeros(num_permutations, 1);
    for j = 1:num_permutations
        perm_indices = randperm(n1 + n2);
        perm_data1 = combined_data(perm_indices(1:n1));
        perm_data2 = combined_data(perm_indices(n1+1:end));
        perm_diffs(j) = mean(perm_data1) - mean(perm_data2);
    end
    
    % 计算 p 值
    p_values(i) = mean(abs(perm_diffs) >= abs(orig_diff));
end

% 显示 p 值
disp('P-values for each category:');
disp(p_values);
