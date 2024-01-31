function tcfilt_OP1( Funcfile, prefix, odir, base, ParamCell )

type = ParamCell{1};

% odir => sprintf('%s/prewarp',opath2f)
% prefix => sprintf('func%u',nr)
% Funcfile => sprintf('%s/func%u%s.nii',opath0,nr,drop_tag)

% prefix for temp files
pref = [odir,'/__opptmp_p2func_tcfilt'];

if ~exist(sprintf('%s/%s_tcfilt.nii.gz',odir,prefix),'file')

    % build directory struct recursively
    unix(sprintf('mkdir -p %s',pref));

    % preparing data for first estimates of displacement...smooth 'n' mask
    unix(sprintf('3dmerge -prefix %s/func_smo6.nii.gz -doall -1blur_fwhm 6 %s',pref,Funcfile));

    %--- b. displacement estimation, for qc and despiking purposes only, since we will be processing a bit before MC'ing afterwards...

    % motion correction for gettting mpes, and then smooth, then re-get (better?) brain mask
    unix(sprintf('3dvolreg -prefix %s/func_mc.nii.gz -1Dfile %s/func_mpe -base %u %s',pref,pref,base,Funcfile));
    unix(sprintf('3dmerge -prefix %s/func_mc_smo6.nii.gz -doall -1blur_fwhm 6 %s/func_mc.nii.gz',pref,pref));
    unix(sprintf('3dAutomask -prefix %s/func_mc_smo6_mask.nii.gz %s/func_mc_smo6.nii.gz',pref,pref)); % ** gets a decent functional mask
%     % dilate the new brain mask for identifying outliers in despiking (NOT FOR NOW?) 
%     unix(sprintf('fslmaths %s/func_mc_smo6_mask.nii.gz -dilD %s/func_mc_smo6_mask_dil.nii.gz',pref,pref));

    if strcmpi(type,'FF')
    
        M = load_untouch_niiz(sprintf('%s/func_mc_smo6_mask.nii.gz',pref));
        V = load_untouch_niiz(sprintf('%s/func_mc_smo6.nii.gz',pref));
    
        volmat = nifti_to_mat(V,M);
        difmat = volmat(:,2:2:end) - volmat(:,1:2:end);
        Ntime = size(difmat,2);
    
        % volume-wise mean and std
        mu_set = mean(difmat,1);
        sg_set = std(difmat,0,1);
        % significance testing --> abnormal mean value, abnormally elevated dispersion
        throut(:,1) =  ( abs(mu_set-mean(mu_set))./std(mu_set) > 2.5);
        throut(:,2) =  (    (sg_set-mean(sg_set))./std(sg_set) > 1.5);
    
        % flag outliers
        outlier_vol =  double( sum(throut,2)>0 );
        % convert to censor files, where 1=non-outlier, and 0=significant outlier
        tmp=  1-outlier_vol;
        outlier_dat.censor_vol = zeros(2*Ntime,1);
        outlier_dat.censor_vol(1:2:end) = tmp; % each pair of scans gets a censor flag
        outlier_dat.censor_vol(2:2:end) = tmp; % if identified
        outlier_dat.metric = type;
        %
        save(sprintf('%s/%s_tcfilt_dat.mat',odir,prefix),'outlier_dat');
    
        % now discarding outlier volumes from the dataset, resaving output
        remove_nii_files( Funcfile, sprintf('%s/%s_tcfilt.nii.gz',odir,prefix), outlier_dat.censor_vol, 1 );
    end

    unix(sprintf('rm -rf %s',pref));
else
    disp('oppni-tcfilt already exists!')
end
