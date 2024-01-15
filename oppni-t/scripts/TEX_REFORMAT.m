close all;
clear;
%% This script converts Freesurfer parcellations into nifti format to facilitate texture analysis

%list_s = [9:13 17:20, 48:56, 16];
list_s = [10:13 17:18, 49:54, 16];
%list_c = [11101:11175, 12101:12175];
list_c = [11101:11141 11143:11175, 12101:12141 12143:12175];
masklab = [list_s list_c];
RAD=5;

% hardcode absolute path
pathloc_fs = 'Neurocovid/NCOV_SUBS';
pathloc_af = 'Neurocovid/fmri_proc';
e = dir(sprintf('%s/sub-*-V01',pathloc_fs));
if isempty(e)
    error('no folders found!');
end

for i=1:numel(e)

    if ~strcmpi(e(i).name,'sub-NCOV1F006_ses-V01')
    
    prfix = e(i).name;
    
        if ~exist( sprintf('texturial_data/%s_seg.nii.gz',prfix), 'file')
    
            prfix,
        
            ranii = sprintf('%s/%s/rawdata/anat1_zclip.nii.gz',pathloc_af,e(i).name);
            ramgz = sprintf('Neurocovid/NCOV_SUBS/%s/mri/orig/001.mgz',e(i).name);
            comgz = sprintf('Neurocovid/NCOV_SUBS/%s/mri/orig.mgz',e(i).name);
            avmgz = sprintf('Neurocovid/NCOV_SUBS/%s/mri/rawavg.mgz',e(i).name);
            asmgz = sprintf('Neurocovid/NCOV_SUBS/%s/mri/aparc.a2009s+aseg.mgz',e(i).name);
            
            % copy over nii
            unix(sprintf('cp %s texturial_data/%s_t1.nii.gz',ranii,prfix));
            % % copy over converte
            % unix(sprintf('mri_convert --in_type mgz --out_type nii %s texturial/%s_1.nii.gz',ramgz,prfix));
            % copy over segge
            unix(sprintf('mri_label2vol --seg %s --temp %s --regheader --o texturial_data/%s_seg.mgz', asmgz, ramgz, prfix));
            unix(sprintf('mri_convert --in_type mgz --out_type nii texturial_data/%s_seg.mgz texturial_data/%s_seg.nii',prfix,prfix));
            unix(sprintf('gzip texturial_data/%s_seg.nii',prfix));
            unix(sprintf('rm texturial_data/%s_seg.mgz',prfix));
        end

    else
        disp('not on list!');
    end

end
