close all;
clear;

output_directory = '/home/mydata/newstudy';

input_file = '/home/mydata/input_file_fmri_undist_blip.txt';
pipe_file  = '/home/mydata/pipe1_file_undist_FL1.txt';
param_file = '/home/mydata/param_file_fmri.txt';

% for AFNI blip correction, use:
% input_file = '/home/mydata/input_file_fmri_undist_blip.txt';
% pipe_file  = '/home/mydata/pipe1_file_undist_AF1.txt';

% for FSL fieldmap correction, use:
% input_file = '/home/mydata/input_file_fmri_undist_fieldmap.txt';
% pipe_file  = '/home/mydata/pipe1_file_undist_FL2.txt';

% populate directory structure, copy over raw data, do some minimal processing
P1_fmri_prepareData( input_file, pipe_file, param_file, output_directory );

% now run full pipeline
P2_fmri_dataProcessing( input_file, pipe_file, param_file, output_directory );

% now running some quality control on the data
P3_fmri_motion_diagnostics( input_file, pipe_file, param_file, output_directory );

return;
