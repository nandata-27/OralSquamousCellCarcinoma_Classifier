%% Unsupervised learning

% Load expression data
data = readcell('Spreadsheets/filtered_counts.xlsx');
metadata = readtable('Spreadsheets/filtered_metadata.xlsx');
row_start = 3;
col_start = 2;

% Convert to numeric matrix
expression_data = cell2mat(data(row_start:end, col_start:end));
sample_data = expression_data';  % samples in rows, genes in columns

%% Run t-SNE
Y = tsne(sample_data, 'NumDimensions', 2, 'Perplexity', 20, 'Standardize', false);

%% Plot t-SNE by disease group 
group_labels = metadata.Group;  

% Custom colors 
colors = [0.3 0.5 1;  
          1.0 0.6 0.2; 
          0.7 0.3 0.9];
% Plot
figure;
gscatter(Y(:,1), Y(:,2), group_labels, colors, [], 25);
xlabel('t-SNE 1');
ylabel('t-SNE 2');
title('t-SNE by Disease Group');
grid on;


%% Unsupervised hierarchical clustering with heatmap 
% Prepare the data for plotting  
%Median center data
medianValue = median(expression_data);
medianCtrData = expression_data-medianValue;


% Define color map
cmap = redbluecmap(9); 
display_range = 10;

% Get sample names and Gene names
row_labels = data(row_start:end, 1);       
col_labels = data(2, col_start:end);       


% Create the hierarchical clustering heatmap
cg = clustergram(medianCtrData, ...
    'RowLabels', row_labels, ...
    'ColumnLabels', col_labels, ...
    'RowPDist', 'euclidean', ...
    'ColumnPDist', 'spearman', ...
    'Colormap', cmap, ...
    'DisplayRange', display_range, ...
    'Linkage', 'average');

%% Feature selection 

%feature selection using fsrstest
%get labels from clinical data
group_label = metadata.("Group");
group_idx = grp2idx(group_label);
feature_names = data(:,1);

[index,score] = fsrftest(medianCtrData', group_idx);
%Infinte values check = 0
idxInf = find(isinf(score));

top_features = feature_names(index);
top_features = cell2table(top_features);

summary(score);

%Plot the most important featues 
hold on
bar(score(index))
xlabel('Predictor rank')
ylabel('Predictor importance score')
title('Most important features discriminating between disease Status - fsrstest')
hold off

%pick top 25 features - fsrftest
features_idx = index(1:25);
top_features_labels = feature_names(index(1:25));


%% hierarchical clustering with heatmap with top 25 features
%get index for top 25 features from the celldata array

gene_names = data(3:end, 1);  
medianCtrWithNames = [gene_names, num2cell(medianCtrData)];

sample_headers = data(1:2, : );
medianCtrdata_clean = [sample_headers; medianCtrWithNames];

%get the expression data for top 25 features
top25_medianCtrData = medianCtrData(features_idx, :);

%get top 25 feature names
top25_feature_names = feature_names(features_idx);

% Define color map
cmap = redbluecmap(9); 

display_range = 10;

%plot clustergram

cg_top25 = clustergram(top25_medianCtrData, ...
    'RowLabels', top25_feature_names, ...
    'ColumnLabels', col_labels, ...
    'RowPDist', 'euclidean', ...
    'ColumnPDist', 'spearman', ...
    'Colormap', cmap, ...
    'DisplayRange', display_range, ...
    'Linkage', 'average');

%prepare data for classification app and statistical analysis
X_numeric = top25_medianCtrData'; 
feature_names_clean = matlab.lang.makeValidName(top25_feature_names);
X_table = array2table(X_numeric, 'VariableNames', feature_names_clean);
X_table.Group = metadata.Group;

%% Reorder features as per ReliefF after supervised learning
% Define the order of gene IDs as per ReliefF
desired_order = {
    'ENSG00000143321'
    'ENSG00000162458'
    'ENSG00000102978'
    'ENSG00000134343'
    'ENSG00000146112'
    'ENSG00000214753'
    'ENSG00000103056'
    'ENSG00000133424'
    'ENSG00000275066'
    'ENSG00000167969'
};

% Convert to valid MATLAB variable names
valid_names = matlab.lang.makeValidName(desired_order);

% Reorder table columns
X_table_reordered = X_table(:, valid_names);
X_table_reordered.Group = metadata.Group;

% Save to CSV
writetable(X_table_reordered, 'Spreadsheets/SPSS_input_top10.csv', 'WriteRowNames', true);

