% This is script is for calculating pearson correlation between pixels from
% cropped single cell ROIs.
%% Quantify correlation
ROIs_dir = 'C:\Users\\xshan\OneDrive\Desktop\matlab_code\correlation_for_colocalization\data\20240504\OE3\Cropped';
pixel_intensity_threshold = 5;
num_cells = 6;
corr_summary = array2table(zeros(num_cells,6), ...
    'VariableNames',{'SPINK1_Golgi','SPINK1_Golgi_shuffled',...
    'COL18A1_Golgi','COL18A1_Golgi_shuffled',...
    'SPINK1_COL18A1','SPINK1_COL18A1_shuffled'});
for i_cell = 1:num_cells
    channel_1 = imread(fullfile(ROIs_dir,sprintf('OE 3_RAW_ch01_cell%d.tif',i_cell)));
    channel_2 = imread(fullfile(ROIs_dir,sprintf('OE 3_RAW_ch02_cell%d.tif',i_cell)));
    channel_3 = imread(fullfile(ROIs_dir,sprintf('OE 3_RAW_ch03_cell%d.tif',i_cell)));
    pixels_to_calculate = channel_3>pixel_intensity_threshold;
    
    channel_1_filtered = double(channel_1(pixels_to_calculate));
    channel_2_filtered = double(channel_2(pixels_to_calculate));
    channel_3_filtered = double(channel_3(pixels_to_calculate));
    % Prepare randomly permutated pixel intensities for negative control
    rng(1);
    channel_1_shuffled = channel_1_filtered(randperm(length(channel_1_filtered)));
    channel_2_shuffled = channel_2_filtered(randperm(length(channel_2_filtered)));
    channel_3_shuffled = channel_3_filtered(randperm(length(channel_3_filtered)));

    % Calculate correlation coefficient for channel 1 and 2
    corr_coef = corrcoef(channel_1_filtered,channel_2_filtered);
    corr_summary{i_cell,1} = corr_coef(1,2);
    % Calculate correlation coefficient for shuffled channel 1 and 2 as
    % negative control
    corr_coef = corrcoef(channel_1_shuffled,channel_2_shuffled);
    corr_summary{i_cell,2} = corr_coef(1,2);

    % Calculate correlation coefficient for channel 2 and 3
    corr_coef = corrcoef(channel_2_filtered,channel_3_filtered);
    corr_summary{i_cell,3} = corr_coef(1,2);
    % Calculate correlation coefficient for shuffled channel 2 and 3 as
    % negative control
    corr_coef = corrcoef(channel_2_shuffled,channel_3_shuffled);
    corr_summary{i_cell,4} = corr_coef(1,2);

    % Calculate correlation coefficient for channel 1 and 3
    corr_coef = corrcoef(channel_1_filtered,channel_3_filtered);
    corr_summary{i_cell,5} = corr_coef(1,2);
    % Calculate correlation coefficient for shuffled channel 2 and 3 as
    % negative control
    corr_coef = corrcoef(channel_1_shuffled,channel_3_shuffled);
    corr_summary{i_cell,6} = corr_coef(1,2);

    % Scatter plot of pixel intensities for channel 1 and 2
    plot(channel_2(pixels_to_calculate),channel_1(pixels_to_calculate),'k.');
    xlabel('GOLGA2');
    ylabel('SPINK1');
    ax = gca;
    exportgraphics(ax,sprintf('results/PANC1_OE3_SPINK1_Golgi_cell%d.jpg',i_cell));
    % Scatter plot of pixel intensities for channel 2 and 3
    plot(channel_2(pixels_to_calculate),channel_3(pixels_to_calculate),'k.');
    xlabel('GOLGA2');
    ylabel('COL18A1');
    ax = gca;
    exportgraphics(ax,sprintf('results/PANC1_OE3_COL18A1_Golgi_cell%d.jpg',i_cell));
    % Scatter plot of pixel intensities for channel 1 and 3
    plot(channel_1(pixels_to_calculate),channel_3(pixels_to_calculate),'k.');
    xlabel('SPINK1');
    ylabel('COL18A1');
    ax = gca;
    exportgraphics(ax,sprintf('results/PANC1_OE3_SPINK1_COL18A1_cell%d.jpg',i_cell));
end
close all;
writetable(corr_summary,'results/PANC1_OE3_correlation_coefficient_summary.csv');
%% Pool all results
for i_sample = 1:3
    if i_sample==1
        merged_results = readtable(sprintf('results/PANC1_OE%d_correlation_coefficient_summary.csv',i_sample));
    else
        to_merge = readtable(sprintf('results/PANC1_OE%d_correlation_coefficient_summary.csv',i_sample));
        merged_results = [merged_results;to_merge];
    end
end
writetable(merged_results,'results/PANC1_all_corr_coef_merged.csv');
%% Calculate mean and sd.
merged_results = readtable('results/MIAPACA/MIAPACA_all_corr_coef_merged.csv');
mean_of_each_corr = mean(table2array(merged_results,1));
sd_of_each_corr = std(table2array(merged_results,1));
p_values = zeros(3,1);
[~,p(1)] = ttest2(merged_results{:,1},merged_results{:,2});
[~,p(2)] = ttest2(merged_results{:,3},merged_results{:,4});
[~,p(3)] = ttest2(merged_results{:,5},merged_results{:,6});