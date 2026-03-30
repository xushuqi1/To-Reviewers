function plot_dea_kde_raw_multi(data, group_names)
% Plot Kernel Density Estimation (KDE) for multiple groups of DEA efficiency values
% (Supports boundaries at 0 and 1)
% Input:
%    data        : Numerical matrix (one group per column) or cell array (one group per element)
%    group_names : Cell array containing names for each data group (optional, defaults to 'Group 1', 'Group 2', ...)

    % Standardize input to a cell array
    if isnumeric(data)
        % Matrix input: each column is a data group
        nGroups = size(data, 2);
        data_cell = cell(1, nGroups);
        for i = 1:nGroups
            col = data(:, i);
            col = col(~isnan(col) & ~isinf(col)); % Clean missing/infinite values
            if isempty(col)
                error('Data in column %d is empty, cannot plot.', i);
            end
            data_cell{i} = col;
        end
    elseif iscell(data)
        data_cell = data;
        nGroups = length(data_cell);
    else
        error('Input "data" must be a numerical matrix or a cell array!');
    end

    % Handle group names
    if nargin < 2 || isempty(group_names)
        group_names = cell(1, nGroups);
        for i = 1:nGroups
            group_names{i} = sprintf('Group %d', i);
        end
    elseif length(group_names) ~= nGroups
        error('Number of group names does not match the number of data groups!');
    end

    % Pre-process all data (truncate to [0, 1])
    all_data = [];
    for i = 1:nGroups
        group_data = data_cell{i};
        % Force truncation to [0, 1]
        group_data = max(min(group_data, 1), 0);
        data_cell{i} = group_data;
        all_data = [all_data; group_data];
    end

    % Output debug info (min/max values for each group)
    fprintf('Data range (after truncation):\n');
    for i = 1:nGroups
        fprintf('  %s : min = %.10f, max = %.10f\n', ...
                group_names{i}, min(data_cell{i}), max(data_cell{i}));
    end

    % Determine x-axis range based on all data
    data_max = max(all_data);
    if data_max > 1
        x_max = min(1.1, data_max + 0.05);
    else
        x_max = 1.05;
    end
    x_min = 0;
    x = linspace(x_min, x_max, 1000);

    % Prepare line styles and colors
    line_styles = {'-', '--', ':', '-.'};
    colors = lines(nGroups);

    % Plotting
    figure;
    hold on;
    for i = 1:nGroups
        group_data = data_cell{i};

        % Adjust values exactly at 0 or 1 to avoid ksdensity support boundary issues
        eps_val = 1e-12;
        data_adj = group_data;
        data_adj(data_adj == 0) = eps_val;
        data_adj(data_adj == 1) = 1 - eps_val;

        % Kernel Density Estimation
        [f, xi] = ksdensity(data_adj, x, 'Support', [0, 1]);

        % Select line style (loop through available styles)
        style_idx = mod(i-1, length(line_styles)) + 1;
        line_style = line_styles{style_idx};

        plot(xi, f, 'LineWidth', 2, 'LineStyle', line_style, ...
             'Color', colors(i,:), 'DisplayName', group_names{i});
    end

    % Figure decoration
    xlabel('Efficiency Score', 'FontSize', 12);
    ylabel('Density', 'FontSize', 12);
    legend('show', 'Location', 'best');
    grid on;
    xlim([x_min, x_max]);
    title('Kernel Density of DEA Efficiency Scores', 'FontSize', 14);
    hold off;
end
