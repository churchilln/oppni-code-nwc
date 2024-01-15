close all;
clear;

output_directory = '/home/mydata/newstudy';

% now run full pipeline (Base1 --> no physio correction)
P1_diff_dataProcessing( '/home/mydata/input_file_diff.txt', output_directory );

% now running some quality control on the data
P2_qc_diffusion( '/home/mydata/input_file_diff.txt', output_directory );

