%% Data Preprocessing

% Read in the data from the excel file as a cell
celldata_raw = readcell('Spreadsheets/raw_counts.csv');

%assign variable for numeric data
row_start=4;
col_start=2;

%% Visualization - Raw Data 
transform_raw_data = log2(replaceZeros(cell2mat(celldata_raw(row_start:end,col_start:end)),"lowval"));

%original raw boxplot
figure;hold on
boxplot(transform_raw_data,'Labels',celldata_raw(1,col_start:end),'LabelOrientation','inline')
title('Boxplot of Raw Data (Log-transformed)')
xlabel('Sample IDs')
ylabel ('Raw Counts')
hold off;

%histogram of raw data to check for skew
histogram(transform_raw_data)
title('Histogram of Raw Data')
xlabel('Bins')
ylabel ('Frequency')
hold off;

%Compute alpha outliers
total_counts_raw = sum(cell2mat(celldata_raw(row_start:end,col_start:end)));
alpha_raw = computeAlphaOutliers(length(total_counts_raw));

% total counts of data
total_counts_raw = sum(cell2mat(celldata_raw(row_start:end,col_start:end)));
scatterplotMarkOutliers(total_counts_raw,...
    'ColumnLabels',celldata_raw(1,col_start:end),'xlabel','SampleID',...
    'ylabel','Total Counts','PlotTitle','TotalRawCounts');

%% Normalization and log transformation of raw data

norm_data = SampleNormalizationRF(cell2mat(celldata_raw(row_start:end,col_start:end)));

%log transformation of normalized data
transform_data_norm_data = log2(replaceZeros(norm_data,"lowval"));

%% Visualization - normalized log transformed raw data

% Boxplots with sample ids as labels.
figure;hold on
boxplot(transform_data_norm_data,'Labels',celldata_raw(1, col_start:end),'LabelOrientation','inline')
title('Boxplot of Normalized Raw Data')
xlabel('Sample IDs')
ylabel ('Normalized Counts')
hold off;

% Histogram
set(gcf,'color','w')
figure; hold on; 
histogram(transform_data_norm_data)
title('Distribution of normalized Raw Data')
xlabel('Read Counts')
ylabel('Frequency')
hold off; 

% mark outliers, total counts of normalized data
total_counts_norm = sum(transform_data_norm_data(row_start:end, col_start:end));
scatterplotMarkOutliers(total_counts_norm,...
    'ColumnLabels',celldata_raw(1,col_start:end),'xlabel','SampleID',...
    'ylabel','Total Counts','PlotTitle','TotalNormalizedCounts');

%%Plot IQR
% Calculate the IQR of the transformed normalized data
iqr_data = iqr(transform_data_norm_data);
[outliers, outlier_ids, low_iqr_outliers] = scatterplotMarkOutliers(iqr_data, ...
    'ColumnLabels', celldata_raw(1, col_start:end), ...
    'PlotTitle', 'IQR of Normalized Raw Data', ...
    'ylabel', 'IQR', ...
    'xlabel', 'Sample IDs');

%remove low outliers
transform_data_norm_data(:,low_iqr_outliers)=[];
iqr_data = iqr(transform_data_norm_data);
[outliers, outlier_ids, low_iqr_outliers] = scatterplotMarkOutliers(iqr_data, ...
    'ColumnLabels', celldata_raw(1, col_start:end), ...
    'PlotTitle', 'IQR of Normalized Raw Data after removing outliers', ...
    'ylabel', 'IQR', ...
    'xlabel', 'Sample IDs');

% Spearman
msc_data = SampleCorrelation(norm_data,'Spearman');

% Plot the mean spearman correlation data with the outlier boundary. The
% vertical lines represent the batches.
[~, ~, low_msc_outliers] = scatterplotMarkOutliers(msc_data, ...
                            'ColumnLabels',celldata_raw(1,col_start:end),...
                            'PlotTitle','Mean Spearman Correlation - Normalized Raw Data','ylabel','Mean Spearman Correlation',...
                            'xlabel','Sample IDs');
%% Filtering data
filtered_idx = MarkLowCounts(norm_data, 0.9);

filtered_data = norm_data(~filtered_idx, :); 


%% Visualization - filtered data

figure;hold on
transform_data_filtered_data = log2(replaceZeros(filtered_data,"lowval"));
boxplot(transform_data_filtered_data,'Labels',celldata_raw(1, col_start:end),'LabelOrientation','inline')
title('Boxplot of Normalized Data')
xlabel('Sample IDs')
ylabel ('Normalized Counts')
hold off;

% IQR
iqr_data_filtered = iqr(transform_data_filtered_data);
[outliers, outlier_ids, low_iqr_outliers] = scatterplotMarkOutliers(iqr_data_filtered, ...
                            'ColumnLabels',celldata_raw(1,col_start:end),...
                            'PlotTitle','IQR - Filtered Data','ylabel','IQR',...
                            'xlabel','Sample IDs');

%remove outliers
transform_data_filtered_data(:,low_iqr_outliers)=[];
iqr_data_removed_out = iqr(transform_data_filtered_data);
[outliers, outlier_ids, low_iqr_outliers] = scatterplotMarkOutliers(iqr_data_removed_out, ...
    'ColumnLabels', celldata_raw(1, col_start:end), ...
    'PlotTitle', 'IQR of Normalized Raw Data after removing outliers', ...
    'ylabel', 'IQR', ...
    'xlabel', 'Sample IDs');

%mean spearman correlation
msc_data_filtered = SampleCorrelation(filtered_data,'Spearman');

% Plot the mean spearman correlation data with the outlier boundary. The
% vertical lines represent the batches.
[~, ~, low_msc_outliers] = scatterplotMarkOutliers(msc_data_filtered, ...
                            'ColumnLabels',celldata_raw(1,col_start:end),...
                            'PlotTitle','Mean Spearman Correlation - Filtered Data','ylabel','Mean Spearman Correlation',...
                            'xlabel','Sample IDs');
%remove outliers and visualize
filtered_data_clean = filtered_data;
filtered_data_clean(:, low_msc_outliers) = [];
filtered_data_clean_log = log2(replaceZeros(filtered_data_clean,"lowval"));

clean_labels = celldata_raw(1, col_start:end); 
clean_labels(:, low_msc_outliers) = []; 
clean_disease = celldata_raw(2, col_start:end); 
clean_disease(:, low_msc_outliers) = []; 

figure; hold on
boxplot(filtered_data_clean_log, 'Labels', clean_labels, 'LabelOrientation', 'inline')
title('Boxplot of Normalized + Filtered Data (after removing low-correlation samples)')
xlabel('Sample IDs')
ylabel('Log2 Normalized Counts')
hold off;


%% Prepare final filtered data for unsupervised and supervised learning
%we will have sample ids and disease state as columns, and gene ids as
%features (rows) 

%convert filtered data matrix into cell array
expression_data = num2cell(filtered_data_clean_log);

%use filtered_idx to filter out low count genes from original gene ids list
all_gene_ids = celldata_raw(row_start:end, 1);  
filtered_gene_ids = all_gene_ids(~filtered_idx);

%get sample ids and associated disease status rows from the original cell
%array
sample_id_row   = [{'GeneID'}, clean_labels];
group_label_row = [{'Group'},  clean_disease];

% add the filtered gene_ids as 1st column
expression_data = [filtered_gene_ids, expression_data]; 

% add the sample ids and disease status as first two rows
filtered_data_cell = [sample_id_row; group_label_row; expression_data];

%save to csv
writecell(filtered_data_cell, 'Spreadsheets/filtered_counts.xlsx')

%save sample_metadata
writecell(filtered_data_cell, 'Spreadsheets/filtered_metadata.xlsx')
