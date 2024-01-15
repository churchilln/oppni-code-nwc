close all;
clear;
% sample file for plotting QC figures, using "runGL" program

e = dir('/Users/tomschweizer/Documents/MATLAB/Neurocovid/fmri_proc/sub-NCOV1F*');
subpipe = 1;
fulpipe = 'Base';

clear ulay olay filn;
for i=1:5 
    ulay{i}=sprintf('/Users/tomschweizer/Documents/MATLAB/Neurocovid/fmri_proc/%s/anat_proc/subpipe_%03u/anat_proc.nii',e(i).name,subpipe);
    olay{i}=sprintf('/Users/tomschweizer/Documents/MATLAB/Neurocovid/fmri_proc/%s/anat_proc/subpipe_%03u/anat_procss.nii',e(i).name,subpipe);
    filn{i}=e(i).name;
end
runGL(ulay,olay,filn,'anat_mask_unwarp');

clear ulay olay filn;
for i=1:6   
    ulay{i}=sprintf('/Users/tomschweizer/Documents/MATLAB/Neurocovid/fmri_proc/%s/anat_proc/subpipe_%03u/anat_warped.nii',e(i).name,subpipe);
    olay{i}=sprintf('/Users/tomschweizer/Documents/MATLAB/Neurocovid/fmri_proc/_group_level/masks/pipe_%s/anat_sulcal_mask_grp.nii',fulpipe);
    filn{i}=e(i).name;
end
runGL(ulay,olay,filn,'anat_warp');

clear ulay olay filn;
for i=1:5 
    ulay{i}=sprintf('/Users/tomschweizer/Documents/MATLAB/Neurocovid/fmri_proc/%s/func_proc_p1/subpipe_%03u/warp/func1_motref.nii',e(i).name,subpipe);
    olay{i}=sprintf('/Users/tomschweizer/Documents/MATLAB/Neurocovid/fmri_proc/%s/func_proc_p1/subpipe_%03u/warp/func1_motref_masked.nii',e(i).name,subpipe);
    filn{i}=e(i).name;
end
runGL(ulay,olay,filn,'func_mask_unwarp');

clear ulay olay filn;
for i=1:5 %numel(e)
    ulay{i}=sprintf('/Users/tomschweizer/Documents/MATLAB/Neurocovid/fmri_proc/%s/func_proc_p1/subpipe_%03u/postwarp/func1_warped_tav.nii',e(i).name,subpipe);
    olay{i}=sprintf('/Users/tomschweizer/Documents/MATLAB/Neurocovid/fmri_proc/_group_level/masks/pipe_%s/func_brain_mask_grp.nii',fulpipe);
    filn{i}=e(i).name;
end
runGL(ulay,olay,filn,'func_warp');


