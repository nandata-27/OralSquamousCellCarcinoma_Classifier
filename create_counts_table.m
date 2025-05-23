%% Script for creating counts matrix from counts files 

%define counts_dir
counts_dir="Spreadsheets/counts";
raw_counts_table = generate_counts_table(counts_dir);

%add metadata information, original sample names
metadata=readtable("Spreadsheets/SraRunTable.csv");

%get column with disease status
disease_row = metadata.Group;

%get column with geo id
sample_id = metadata.SampleName;

%get column with run_info
run = metadata.Run;

%transpose them into rows 
disease_row=disease_row';
sample_id=sample_id';
run = run';

%add lables in first cell 
disease_row=[{'Disease_status'}, disease_row];
sample_id=[{'GEO_ID'}, sample_id];
run = [{'Run_id'}, run];

%convert the raw counts table into a cell array
counts_array = table2cell(raw_counts_table);

%add the metadata information to the counts array
final_rawdata = [sample_id; disease_row; run; counts_array];

%write the rawdata cell array to a csv
writecell(final_rawdata, 'Spreadsheets/raw_counts.csv');

%raw counts data is ready for preprocessing and visualization and quality
%control