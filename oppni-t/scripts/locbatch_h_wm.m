close all;
clear;
disp('go onn'); % batching for fixed BW values

e = dir('texturial_data/sub-*t1.nii');

for i=1:numel(e)
    if ~strcmpi(e(i).name,'sub-NCOV1F006_ses-V01')
        % actual command to execute
        prfix = e(i).name(1:21);    
        t1_file = sprintf('texturial_data/%s_t1.nii',prfix);
        mskfile = sprintf('texturial_data/%s_seg.nii',prfix);
        outfile = sprintf('texturial_matfiles_vox/%s_wm',prfix);
        texmap_wrapper3( t1_file, mskfile, [2 41], 1, outfile, 7, 7, [], [], 'HDE', [0.10] );
    else
        disp('not on list!');
    end
end

