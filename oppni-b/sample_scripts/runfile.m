close all;
clear;

output_directory = '/home/mydata/newstudy';

% populate directory structure, copy over raw data, do some minimal processing
P1_fmri_prepareData( '/home/mydata/input_file_fmri.txt', '/home/mydata/pipe1_file.txt', '/home/mydata/param_file_fmri.txt', output_directory );

% now run full pipeline (Base1 --> no physio correction)
P2_fmri_dataProcessing( '/home/mydata/input_file_fmri.txt', '/home/mydata/pipe1_file.txt', '/home/mydata/param_file_fmri.txt', output_directory );

% now running some quality control on the data
P3_fmri_motion_diagnostics( '/home/mydata/input_file_fmri.txt', '/home/mydata/pipe1_file.txt', '/home/mydata/param_file_fmri.txt', output_directory );

return;
%%

% some scripts for quick visual inspection of output data quality

InputStruct = Read_Input_File_fmri('/home/mydata/input_file_fmri.txt'); 

% 1.
config = 'anat_mask_unwarp';
pipe = 'subpipe_001';
for i=1:20
    imname{i} = InputStruct(i).PREFIX;
    bg{i} = sprintf('%s/fmri_proc/%s/anat_proc/anat1_2std.nii.gz',output_directory,imname{i})
    overlay{i} = sprintf('%s/fmri_proc/%s/anat_proc/%s/anat_procss.nii.gz',output_directory,imname{i},pipe)
end
runGL(bg, overlay, imname, config);

% 2.
config = 'func_mask_unwarp';
pipe = 'subpipe_001';
for i=1
    imname{i} = InputStruct(i).PREFIX;
    bg{i} = sprintf('%s/fmri_proc/%s/func_proc_p1/%s/warp/func1_motref.nii.gz',output_directory,imname{i},pipe)
    overlay{i} = sprintf('%s/fmri_proc/%s/func_proc_p1/%s/warp/func1_motref_masked.nii.gz',output_directory,imname{i},pipe)
end
runGL(bg, overlay, imname, config);

% 3.
config = 'anat_warp';
pipe = 'subpipe_001';
pipe2 = 'pipe_Base2';
for i=1
    imname{i} = InputStruct(i).PREFIX;
    bg{i} = sprintf('%s/fmri_proc/%s/anat_proc/%s/anat_warped.nii.gz',output_directory,imname{i},pipe)
    overlay{i} = sprintf('%s/fmri_proc/_group_level/masks/%s/anat_sulcal_mask_grp.nii',output_directory,pipe2)
end
runGL(bg, overlay, imname, config);

% 4.
config = 'func_warp';
pipe = 'subpipe_001';
pipe2 = 'pipe_Base2';
for i=1
    imname{i} = InputStruct(i).PREFIX;
    bg{i} = sprintf('%s/fmri_proc/%s/func_proc_p1/%s/postwarp/func1_warped_tav.nii.gz',output_directory,imname{i},pipe)
    overlay{i} = sprintf('%s/fmri_proc/_group_level/masks/%s/func_brain_mask_grp.nii',output_directory,pipe2)
end
runGL(bg, overlay, imname, config);
