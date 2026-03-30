function plot_dea_kde_raw_multi(data, group_names)
% 绘制多组 DEA 效率值的核密度图（支持边界0和1）
% 输入：
%   data        : 数值矩阵（每列一组数据）或元胞数组（每个元素为一组数据）
%   group_names : 元胞数组，对应每组数据的名称（可选，默认生成 'Group 1','Group 2',...）

    % 统一将输入转换为元胞数组
    if isnumeric(data)
        % 矩阵输入：每列一组数据
        nGroups = size(data, 2);
        data_cell = cell(1, nGroups);
        for i = 1:nGroups
            col = data(:, i);
            col = col(~isnan(col) & ~isinf(col)); % 清理缺失值
            if isempty(col)
                error('第 %d 列数据为空，无法绘图。', i);
            end
            data_cell{i} = col;
        end
    elseif iscell(data)
        data_cell = data;
        nGroups = length(data_cell);
    else
        error('data 必须是数值矩阵或元胞数组！');
    end

    % 处理组名
    if nargin < 2 || isempty(group_names)
        group_names = cell(1, nGroups);
        for i = 1:nGroups
            group_names{i} = sprintf('Group %d', i);
        end
    elseif length(group_names) ~= nGroups
        error('组名数量与数据组数不一致！');
    end

    % 预处理所有数据（截断到 [0,1]）
    all_data = [];
    for i = 1:nGroups
        data = data_cell{i};
        % 强制截断到 [0,1]
        data = max(min(data, 1), 0);
        data_cell{i} = data;
        all_data = [all_data; data];
    end

    % 输出调试信息（每组的最小最大值）
    fprintf('数据范围（截断后）：\n');
    for i = 1:nGroups
        fprintf('  %s : min = %.10f, max = %.10f\n', ...
                group_names{i}, min(data_cell{i}), max(data_cell{i}));
    end

    % 确定横轴范围（基于所有数据）
    data_max = max(all_data);
    if data_max > 1
        x_max = min(1.1, data_max + 0.05);
    else
        x_max = 1.05;
    end
    x_min = 0;
    x = linspace(x_min, x_max, 1000);

    % 准备线型和颜色
    line_styles = {'-', '--', ':', '-.'};
    colors = lines(nGroups);

    % 绘图
    figure;
    hold on;
    for i = 1:nGroups
        data = data_cell{i};

        % 为避免 ksdensity 边界检查失败，将等于 0 或 1 的值微调
        eps_val = 1e-12;
        data_adj = data;
        data_adj(data_adj == 0) = eps_val;
        data_adj(data_adj == 1) = 1 - eps_val;

        % 核密度估计
        [f, xi] = ksdensity(data_adj, x, 'Support', [0, 1]);

        % 选择线型（循环使用）
        style_idx = mod(i-1, length(line_styles)) + 1;
        line_style = line_styles{style_idx};

        plot(xi, f, 'LineWidth', 2, 'LineStyle', line_style, ...
             'Color', colors(i,:), 'DisplayName', group_names{i});
    end

    % 图形装饰
    xlabel('Efficiency Score', 'FontSize', 12);
    ylabel('Density', 'FontSize', 12);
    legend('show', 'Location', 'best');
    grid on;
    xlim([x_min, x_max]);
    title('Kernel Density of DEA Efficiency Scores', 'FontSize', 14);
    hold off;
end