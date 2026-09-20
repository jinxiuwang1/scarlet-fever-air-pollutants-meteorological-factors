clc; clear;

%% Actual data (taking Gansu as an example)
[AA, BB] = xlsread('Number of scarlet fever cases in 31 provinces and municipalities from 2013 to 2020.xlsx');

Label = {'Hainan','Guangdong','Guangxi','Yunnan','Fujian','Guizhou','Hunan','Jiangxi','Zhejiang','Chongqing','Sichuan','Hubei','Shanghai','Tibet','Anhui','Jiangsu','Henan','Shaanxi','Qinghai','Shandong','Ningxia','Shanxi','Gansu','Tianjin','Hebei','Beijing','Xinjiang','Liaoning','Jilin','Inner Mongolia','Heilongjiang'};
Lambda_new = [10882,94998,56156,49214,38369,38304,75274,49704,45862,25666,66883,53544,16462,4100,64722,62453,96248,31398,6820,92545,7150,32700,26164,10157,79685,15739,29885,22279,12288,18693,21923];
Total_P = [8950000,106440000,47190000,46870000,37740000,35020000,66910000,45220000,54980000,29700000,81070000,57990000,24150000,3120000,60300000,79390000,94130000,37640000,5780000,97330000,8950000,36300000,25820000,14720000,73330000,21150000,22640000,43900000,27510000,24980000,38350000];
Total_P_0_14 = [1677230,17754192,10131693,9383374,6359190,7735918,12224457,8908340,6564612,3979800,13417085,8669505,2260440,725712,11131380,10582687,19550801,5676112,1149064,15066684,1323696,5706360,4366162,1647168,13052740,2098080,4715912,4574380,3232425,3429754,4594330];

% Store calculation results for all provinces
for i = 1:31
    [SF, BETA, new_cases, dt, n, time, S, N] = function_beta(Label, AA, BB, i, Lambda_new, Total_P, Total_P_0_14);
    SF_cases(:,i) = SF;
    BETA_cases(:,i) = BETA;
    New_cases(:,i) = new_cases;
    S_cases(:,i) = S;
    N_cases(:,i) = N;
end

% Smooth monthly beta with different windows 
% BETA_cases: original unsmoothed beta(t)
BETA_smooth_2 = movmean(BETA_cases, 2);
BETA_smooth_3 = movmean(BETA_cases, 3);
BETA_smooth_4 = movmean(BETA_cases, 4);
BETA_smooth_5 = movmean(BETA_cases, 5);

% Keep the 3-month moving average as the main analysis
BETA_smooth = BETA_smooth_3;


writematrix(BETA_cases,      'BETA_cases.xlsx');
writematrix(BETA_smooth_2, 'BETA_smooth_2.xlsx');
writematrix(BETA_smooth_3, 'BETA_smooth_3.xlsx');
writematrix(BETA_smooth_4, 'BETA_smooth_4.xlsx');
writematrix(BETA_smooth_5, 'BETA_smooth_5.xlsx');


% Effective reproduction number parameters
sigma = 30/3;
gamma = 30/6.25;
g = 1/(15*12);
d = 1/(76*12);

% Compute Re and annual average Re
for i = 1:31
    Re(:,i) = (BETA_cases(:,i) .* sigma .* Lambda_new(i)/(d+g)) ./ ((sigma + g + d)*(gamma + g + d)*(Lambda_new(i)/d));
    Re_mean(:,i) = sum(reshape(Re(:,i), 12, [])) / 12;
end

% ========== Reorder provinces ==========
special_idx = [1, 20, 13, 19];   % Hainan, Shandong, Shanghai, Qinghai
special_names = {'Hainan','Shandong','Shanghai','Qinghai'};

all_idx = 1:31;
other_idx = setdiff(all_idx, special_idx);
[~, order] = sort(other_idx);
other_idx = other_idx(order);

% Base names of other provinces (without number prefix)
other_names = Label(other_idx);

% ========== Generate labels for the three types of figures ==========
% Special province labels (Fig case A1~A4, Fig Re B1~B4, Fig Beta C1~C4)
special_labels_case = cellfun(@(name, i) sprintf('($A_{%d}$) %s', i, name), special_names, num2cell(1:4), 'UniformOutput', false);
special_labels_re   = cellfun(@(name, i) sprintf('($B_{%d}$) %s', i, name), special_names, num2cell(1:4), 'UniformOutput', false);
special_labels_beta = cellfun(@(name, i) sprintf('($C_{%d}$) %s', i, name), special_names, num2cell(1:4), 'UniformOutput', false);

% Other provinces labels: numbers from 5 to 31
other_labels_case = cell(1, 27);
other_labels_re   = cell(1, 27);
other_labels_beta = cell(1, 27);
for k = 1:27
    idx_num = k + 4;   % 5,6,...,31
    other_labels_case{k} = sprintf('($A_{%d}$) %s', idx_num, other_names{k});
    other_labels_re{k}   = sprintf('($B_{%d}$) %s', idx_num, other_names{k});
    other_labels_beta{k} = sprintf('($C_{%d}$) %s', idx_num, other_names{k});
end

% Divide 27 provinces into 7 groups (each at most 4)
group_size = 4;
num_groups = ceil(27 / group_size);
groups = cell(num_groups, 1);
for g = 1:num_groups
    start_idx = (g-1)*group_size + 1;
    end_idx = min(g*group_size, 27);
    groups{g} = start_idx:end_idx;
end

% Uniform subplot spacing
gap_vert = 0.12;
gap_horiz = 0.08;
marg_top = 0.12;
marg_bottom = 0.08;
marg_left = 0.06;
marg_right = 0.06;

% ========== Plotting section ==========
% Figure 1_A: Case plots for four special provinces
fig1_A = figure('Name', 'Fig1_Special_Provinces', 'NumberTitle', 'off');
[ha_A, ~] = tight_subplot(2, 2, [gap_vert, gap_horiz], [marg_top, marg_bottom], [marg_left, marg_right]);
for k = 1:4
    idx = special_idx(k);
    ax = ha_A(k);
    add_legend_flag = (k == 2);
    plot_case_subplot(ax, time, SF_cases(:,idx), New_cases(:,idx), BETA_cases(:,idx), dt, n, special_labels_case{k}, 14, 14, add_legend_flag);
end
add_global_ylabels(fig1_A, '$\mathrm{New\ cases}$', '$\mathrm{Transmission\ rate}\ (\beta(t))$', 24, 24);

% Figure 2_A: Annual average Re plots for four special provinces
fig2_A = figure('Name', 'Fig2_Special_Provinces', 'NumberTitle', 'off');
[ha_A2, ~] = tight_subplot(2, 2, [gap_vert, gap_horiz], [marg_top, marg_bottom], [marg_left, marg_right]);
for k = 1:4
    idx = special_idx(k);
    ax = ha_A2(k);
    plot_Re_subplot(ax, Re_mean(:,idx), special_labels_re{k}, 14, 14);
end
add_global_leftlabel(fig2_A, '$\mathrm{The\ annual\ average\ value\ of}\ R_{e}(t)$', 0.02, 'k', 24);

% Figure 3_A: Smoothed beta plots for four special provinces
fig3_A = figure('Name', 'Fig3_Special_Provinces', 'NumberTitle', 'off');
[ha_A3, ~] = tight_subplot(2, 2, [gap_vert, gap_horiz], [marg_top, marg_bottom], [marg_left, marg_right]);
for k = 1:4
    idx = special_idx(k);
    ax = ha_A3(k);
    plot_beta_subplot(ax, time, BETA_smooth(:,idx), special_labels_beta{k}, 14, 14);
end
add_global_leftlabel(fig3_A, '$\mathrm{Smoothed}\ \beta(t)$', 0.02, 'k', 24);

% Figure 1_B: Case plots for remaining provinces (7 groups)
for g = 1:num_groups
    fig = figure('Name', sprintf('Fig1_B_Group%d', g), 'NumberTitle', 'off');
    [ha, ~] = tight_subplot(2, 2, [gap_vert, gap_horiz], [marg_top, marg_bottom], [marg_left, marg_right]);
    group_provinces = groups{g};
    num_in_group = length(group_provinces);
    for sub_i = 1:4
        if sub_i <= num_in_group
            idx_in_other = group_provinces(sub_i);
            idx = other_idx(idx_in_other);
            add_legend_flag = (sub_i == 2);
            plot_case_subplot(ha(sub_i), time, SF_cases(:,idx), New_cases(:,idx), BETA_cases(:,idx), dt, n, other_labels_case{idx_in_other}, 14, 14, add_legend_flag);
        else
            set(ha(sub_i), 'Visible', 'off');
        end
    end
    add_global_ylabels(fig, '$\mathrm{New\ cases}$', '$\mathrm{Transmission\ rate}\ \beta(t)$', 24, 24);
end

% Figure 2_B: Annual average Re plots for remaining provinces (7 groups)
for g = 1:num_groups
    fig = figure('Name', sprintf('Fig2_B_Group%d', g), 'NumberTitle', 'off');
    [ha, ~] = tight_subplot(2, 2, [gap_vert, gap_horiz], [marg_top, marg_bottom], [marg_left, marg_right]);
    group_provinces = groups{g};
    num_in_group = length(group_provinces);
    for sub_i = 1:4
        if sub_i <= num_in_group
            idx_in_other = group_provinces(sub_i);
            idx = other_idx(idx_in_other);
            plot_Re_subplot(ha(sub_i), Re_mean(:,idx), other_labels_re{idx_in_other}, 14, 14);
        else
            set(ha(sub_i), 'Visible', 'off');
        end
    end
    add_global_leftlabel(fig, '$\mathrm{The\ annual\ average\ value\ of}\ R_{e}(t)$', 0.02, 'k', 24);
end

% Figure 3_B: Smoothed beta plots for remaining provinces (7 groups)
for g = 1:num_groups
    fig = figure('Name', sprintf('Fig3_B_Group%d', g), 'NumberTitle', 'off');
    [ha, ~] = tight_subplot(2, 2, [gap_vert, gap_horiz], [marg_top, marg_bottom], [marg_left, marg_right]);
    group_provinces = groups{g};
    num_in_group = length(group_provinces);
    for sub_i = 1:4
        if sub_i <= num_in_group
            idx_in_other = group_provinces(sub_i);
            idx = other_idx(idx_in_other);
            plot_beta_subplot(ha(sub_i), time, BETA_smooth(:,idx), other_labels_beta{idx_in_other}, 14, 14);
        else
            set(ha(sub_i), 'Visible', 'off');
        end
    end
    add_global_leftlabel(fig, '$\mathrm{Smoothed}\ \beta(t)$', 0.02, 'k', 24);
end

% ==================== Local function definitions ====================
function plot_case_subplot(ax, time, SF, New, BETA, dt, n, title_str, fs_axis, fs_title, add_legend_flag)
    axes(ax);
    yyaxis left;
  
    h1 = plot(time, SF, 'ko', 'LineWidth', 1.2, 'MarkerSize', 4, 'MarkerFaceColor', 'none');
    hold on;

    h2 = plot([1:dt:n], New, 'r--', 'LineWidth', 1.2);
    ylabel('');
    set(gca, 'xtick', 1:12:length(time));
    set(gca, 'xticklabel', {'2013','2014','2015','2016','2017','2018','2019','2020'});
    set(gca, 'FontSize', fs_axis);
    xlim([0.5 length(time)+0.5]);
    ax.YColor = 'r';
    
    yyaxis right;

    h3 = plot(1:n, BETA, 'b-', 'LineWidth', 2.0);
    ylabel('');
    xlim([0.5 length(time)+0.5]);
    ax.YColor = 'b';
    title(title_str, 'Interpreter', 'latex', 'FontSize', fs_title, 'FontWeight', 'bold');
    
    if add_legend_flag
        legend([h1, h2, h3], ...
            {'$\mathrm{New\ cases}$', '$\mathrm{Simulated\ cases}$', '$\mathrm{Transmission\ rate}$'}, ...
            'FontSize', 14, 'Orientation', 'horizontal', 'Box', 'off', ...
            'Location', 'northeast', 'Interpreter', 'latex');
    end
end

function plot_Re_subplot(ax, Re_mean, title_str, fs_axis, fs_title)
    axes(ax);
    plot(1:8, Re_mean, 'b-*', 'LineWidth', 1.2, 'MarkerSize', 5);
    ylabel('');
    xlabel('');
    set(gca, 'xtick', 1:8);
    set(gca, 'xticklabel', {'2013','2014','2015','2016','2017','2018','2019','2020'});
    set(gca, 'FontSize', fs_axis);
    xlim([0.8 8.2]);
    title(title_str, 'Interpreter', 'latex', 'FontSize', fs_title, 'FontWeight', 'bold');
end

function plot_beta_subplot(ax, time, BETA_smooth, title_str, fs_axis, fs_title)
    axes(ax);
    plot(time, BETA_smooth, 'b-', 'LineWidth', 1.5);
    xlim([0.5 length(time)+0.5]);
    xticks(1:12:length(time));
    set(gca, 'xticklabel', {'2013','2014','2015','2016','2017','2018','2019','2020'});
    set(gca, 'FontSize', fs_axis);
    title(title_str, 'Interpreter', 'latex', 'FontSize', fs_title, 'FontWeight', 'bold');
end

function add_global_ylabels(fig, left_str, right_str, fs_left, fs_right)
    hidden_ax = axes('Parent', fig, 'Position', [0,0,1,1], 'Visible', 'off', 'HandleVisibility', 'off');
    text(hidden_ax, 0.02, 0.5, left_str, 'Units', 'normalized', ...
        'FontSize', fs_left, 'FontWeight', 'bold', 'Color', 'r', ...
        'Rotation', 90, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'Interpreter', 'latex');
    text(hidden_ax, 0.98, 0.5, right_str, 'Units', 'normalized', ...
        'FontSize', fs_right, 'FontWeight', 'bold', 'Color', 'b', ...
        'Rotation', 90, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'Interpreter', 'latex');
    uistack(hidden_ax, 'bottom');
end

function add_global_leftlabel(fig, label_str, xpos, color, fs)
    hidden_ax = axes('Parent', fig, 'Position', [0,0,1,1], 'Visible', 'off', 'HandleVisibility', 'off');
    text(hidden_ax, xpos, 0.5, label_str, 'Units', 'normalized', ...
        'FontSize', fs, 'FontWeight', 'bold', 'Color', color, ...
        'Rotation', 90, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'Interpreter', 'latex');
    uistack(hidden_ax, 'bottom');
end