function despike_OP1( Funcfile, prefix, odir, type, base )

% odir => sprintf('%s/prewarp',opath2f)
% prefix => sprintf('func%u',nr)
% Funcfile => sprintf('%s/func%u%s.nii',opath0,nr,drop_tag)

% prefix for temp files
pref = [odir,'/__opptmp_p2func_despike'];

if ~exist(sprintf('%s/%s_despike.nii.gz',odir,prefix),'file')

    % build directory struct recursively
    unix(sprintf('mkdir -p %s',pref));

    % preparing data for first estimates of displacement...smooth 'n' mask
    unix(sprintf('3dmerge -prefix %s/func_smo6.nii.gz -doall -1blur_fwhm 6 %s',pref,Funcfile));

    %--- b. displacement estimation, for qc and despiking purposes only, since we will be processing a bit before MC'ing afterwards...

    % motion correction for gettting mpes, and then smooth, then re-get (better?) brain mask
    unix(sprintf('3dvolreg -prefix %s/func_mc.nii.gz -1Dfile %s/func_mpe -base %u %s',pref,pref,base,Funcfile));
    unix(sprintf('3dmerge -prefix %s/func_mc_smo6.nii.gz -doall -1blur_fwhm 6 %s/func_mc.nii.gz',pref,pref));
    unix(sprintf('3dAutomask -prefix %s/func_mc_smo6_mask.nii.gz %s/func_mc_smo6.nii.gz',pref,pref)); % ** gets a decent functional mask
    % dilate the new brain mask for identifying outliers in despiking
    unix(sprintf('fslmaths %s/func_mc_smo6_mask.nii.gz -dilD %s/func_mc_smo6_mask_dil.nii.gz',pref,pref));

    %-- a. collecting outlier estimates

    % done on non-motcorred data; want to find the *really* extreme bad cases
    outlier_dat = spike_estimator( sprintf('%s/func_smo6.nii.gz',pref), sprintf('%s/func_mc_smo6_mask_dil.nii.gz',pref), sprintf('%s/func_mpe',pref), ['--unused--'], 7,1 );
    save(sprintf('%s/%s_despike_dat.mat',odir,prefix),'outlier_dat');

    % replacing volumes that are both extreme BOLD values and correspond to high-estimated-motion timepoints
    fmri_interpolator( Funcfile, outlier_dat, type, sprintf('%s/%s_despike.nii.gz',odir,prefix) )

    unix(sprintf('rm -rf %s',pref));
else
    disp('afni-despike already exists!')
end
