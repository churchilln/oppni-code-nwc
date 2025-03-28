function amask_FS1(Adataset, Basedset, odir, ParamCell)
%
% Freesurfer-based anatomical brain masking

% Adataset = 'inputdataset.nii';  %# req/ input dataset
% SubID    = 'subjID';            %# req/ the subject ID
% Basedset = 'basedset.nii';       %# req/ reference dset- must have 4 bricks
% odir     = 'output';             %# opt/ output dir

doclean    = 1;       %# def=1; clean out junk files after terminated; can also set =0 to not clean
watershed_threshold  = 25;        %# def=25; "Increasing the threshold will increase the likelihood that both brain and skull will remain" vals [0 50]

% prefix for temp files
pref = [odir,'/__opptmp_p2anat_mask'];

if ~exist( sprintf('%s/anatBrainMask.nii.gz',odir) ,'file')


    % check for path  to base file, create tpath
    if exist(Basedset,'file')
        [tpath,~,~] = fileparts(Basedset);
    else
        error('cannot find path to Based file');
    end
    
    % Require it to have enough bricks to really be a ref
    [~,nvolbase] = unix(sprintf('3dinfo -nv %s',Basedset)); nvolbase = str2num(nvolbase);
    if nvolbase < 5
        error('needs 5 volumes to serve as base!?!?!');
    end

    % run freesurfer autorecon1; note that the dataset will be cropped to
    % 256x256x256 voxel resolution
    disp('these are the FS1 commands:\n');
    fprintf('recon-all -autorecon1 -cw256 -wsthresh %d -clean-bm -sd %s -s __opptmp_p2anat_mask -i %s\n',watershed_threshold,odir,Adataset);
    fprintf('mri_convert --in_type mgz --out_type nii %s/mri/brainmask.mgz %s/mri/brainmask.nii\n',pref,pref);

    unix(sprintf('recon-all -autorecon1 -cw256 -wsthresh %d -clean-bm -sd %s -s __opptmp_p2anat_mask -i %s',watershed_threshold,odir,Adataset));

    % generate mask from map
    unix(sprintf('mri_convert --in_type mgz --out_type nii %s/mri/brainmask.mgz %s/mri/brainmask.nii',pref,pref));
    unix(sprintf('3dAutomask -prefix %s/mri/anatBrainMask %s/mri/brainmask.nii',pref,pref));

    % resample back into original space and convert to nifti
    unix(sprintf('3dresample -input %s/mri/anatBrainMask+orig.HEAD -master %s -prefix %s/mri/anatBrainMask_resampled',pref,Adataset,pref));
    unix(sprintf('3dAFNItoNIFTI -prefix %s/mri/anatBrainMask %s/mri/anatBrainMask_resampled+orig.HEAD',pref,pref));
    unix(sprintf('gzip %s/mri/anatBrainMask.nii',pref));
    unix(sprintf('cp %s/mri/anatBrainMask.nii.gz %s/anatBrainMask.nii.gz',pref,odir));
    
error('Deliberately Stopped Here...\n');
    
    %# ------------------------------------------------------------------
    %## MODIFIED=> Clean up the junk, keeps only anatBrainMask
    
    if exist(sprintf('%s/anatBrainMask.nii.gz',odir),'file') 
        if doclean>0
            unix(sprintf('rm -rf %s',pref));
        end
    else
        error('freesurfer failure to create anatomical mask')
    end

else
    disp('freesurfermask already exists!')
end
% % The Template dataset ~2~
% % 
% %   Any reference base template dataset, such as
% %   MNI152_2009_template_SSW.nii.gz, must have the first *4* volumes here
% %   (and can have the optional 5th for later uses, as described):
% %     [0] = skull-stripped template brain volume
% %     [1] = skull-on template brain volume
% %     [2] = weight mask for nonlinear registration, with the
% %           brain given greater weight than the skull
% %     [3] = binary mask for the brain
% %     [4] = binary mask for gray matter plus some CSF (slightly dilated)
% %           ++ this volume is not used in this script
% %           ++ it is intended for use in restricting FMRI analyses
% %              to the 'interesting' parts of the brain
% %           ++ this mask should be resampled to your EPI spatial
% %              resolution (see program 3dfractionize), and then
% %              combined with a mask from your experiment reflecting
% %              your EPI brain coverage (see program 3dmask_tool).
