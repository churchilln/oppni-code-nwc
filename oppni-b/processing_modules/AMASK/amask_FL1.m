function amask_FL1(Adataset, Basedset, odir, ParamCell)


% prefix for temp files
pref = [odir,'/__opptmp_p2anat_mask'];

if ~exist( sprintf('%s/anatBrainMask.nii.gz',odir) ,'file')


    % build directory struct recursively
    unix(sprintf('mkdir -p %s',pref));
    % check for path  to base file, create tpath
    if exist(Basedset,'file')
        [tpath,~,~] = fileparts(Basedset);
    else
        error('cannot find path to Based file');
    end



    % make a copy
    if contains(Adataset,'.nii.gz')
        unix(sprintf('cp %s %s/anat.nii.gz',Adataset,pref));
    elseif contains(Adataset,'.nii')
        unix(sprintf('cp %s %s/anat.nii',Adataset,pref));
        unix(sprintf('gzip %s/anat.nii',pref));
    else
        error('unregocnized Adataset format');
    end



    % Require it to have enough bricks to really be a ref
    [~,nvolbase] = unix(sprintf('3dinfo -nv %s',Basedset)); nvolbase = str2num(nvolbase);
    if nvolbase < 4
        error('needs 4 volumes to serve as base!?!?!');
    end

    [~,refp, refe]=fileparts( Basedset );
    if strcmp(refe,'.nii')
        unix(sprintf('cp %s %s/Base.nii',Basedset,pref))
        unix(sprintf('gzip %s/Base.nii',pref));
    elseif strcmp(refe,'.gz')
        unix(sprintf('cp %s %s/Base.nii.gz',Basedset,pref))
    else
        error('template has unrecognized file ending!')
    end
    unix(sprintf('3dTcat -prefix %s/Base_skullon.nii.gz ''%s/Base.nii.gz[1]''',pref,pref));
    unix(sprintf('3dTcat -prefix %s/Base_mask.nii.gz ''%s/Base.nii.gz[3]''',pref,pref));



    % trim FOV and also get the inverse transform
    % deviates from FSL script, since we already to fslreorient outside of script 
    unix(sprintf('robustfov -i %s/anat.nii.gz -r %s/T1.nii.gz -m %s/T1_roi2orig.mat',pref,pref,pref));

    unix(sprintf('fslmaths %s/T1 -mul 0 %s/lesionmask',pref,pref))
    unix(sprintf('fslmaths %s/lesionmask -bin %s/lesionmask',pref,pref))
    unix(sprintf('fslmaths %s/lesionmask -binv %s/lesionmaskinv',pref,pref))

    % initial masking - very rough / over-inclusive
    unix(sprintf('bet %s/T1 %s/T1_initfast2_brain -m -f 0.1',pref,pref));
    % copying... [trim later?]
    unix(sprintf('fslmaths %s/T1_initfast2_brain %s/T1_initfast2_restore',pref,pref));
    % masking? > with all-ones volume! [can trim later] / tissue segment
    unix(sprintf('fslmaths %s/T1_initfast2_restore -mas %s/lesionmaskinv %s/T1_initfast2_maskedrestore',pref,pref,pref)); %*
    unix(sprintf('fast -o %s/T1_fast -l 20 -b -B -t 1 --iter=10 --nopve --fixed=0 -v %s/T1_initfast2_maskedrestore',pref,pref));

    % extrapolating bias field from central region
    unix(sprintf('fslmaths %s/T1 -div %s/T1_fast_restore -mas %s/T1_initfast2_brain_mask %s/T1_fast_totbias',pref,pref,pref,pref));
    unix(sprintf('fslmaths %s/T1_initfast2_brain_mask -ero -ero -ero -ero -mas %s/lesionmaskinv %s/T1_initfast2_brain_mask2',pref,pref,pref));
    unix(sprintf('fslmaths %s/T1_fast_totbias -sub 1 %s/T1_fast_totbias',pref,pref));
    unix(sprintf('fslsmoothfill -i %s/T1_fast_totbias -m %s/T1_initfast2_brain_mask2 -o %s/T1_fast_bias',pref,pref,pref));
    unix(sprintf('fslmaths %s/T1_fast_bias -add 1 %s/T1_fast_bias',pref,pref));
    unix(sprintf('fslmaths %s/T1_fast_totbias -add 1 %s/T1_fast_totbias',pref,pref));
    unix(sprintf('fslmaths %s/T1 -div %s/T1_fast_bias %s/T1_biascorr',pref,pref,pref));
    
    % registration/brain extraction
%     unix(sprintf('flirt -interp spline -dof 12 -in %s/T1_biascorr -ref /Users/tomschweizer/FSL/data/standard/MNI152_T1_2mm -dof 12 -omat %s/T1_to_MNI_lin.mat -out %s/T1_to_MNI_lin',pref,pref,pref));
    unix(sprintf('flirt -interp spline -dof 12 -in %s/T1_biascorr -ref %s/Base_skullon.nii.gz -dof 12 -omat %s/T1_to_MNI_lin.mat -out %s/T1_to_MNI_lin',pref,pref,pref,pref));
    
    % nonlinear step
    unix(sprintf('fslmaths %s/Base_mask.nii.gz -fillh -dilF %s/MNI152_T1_2mm_brain_mask_dil1',pref,pref));
    unix(sprintf('fnirt --in=%s/T1_biascorr --ref=%s/Base_skullon.nii.gz --fout=%s/T1_to_MNI_nonlin_field --jout=%s/T1_to_MNI_nonlin_jac --iout=%s/T1_to_MNI_nonlin --logout=%s/T1_to_MNI_nonlin.txt --cout=%s/T1_to_MNI_nonlin_coeff --config=/Users/tomschweizer/FSL/etc/flirtsch/T1_2_MNI152_2mm.cnf --aff=%s/T1_to_MNI_lin.mat --refmask=%s/MNI152_T1_2mm_brain_mask_dil1',...
        pref,pref,pref,pref,pref,pref,pref,pref,pref))
    % doing the brain extraction...
    unix(sprintf('invwarp --ref=%s/T1_biascorr -w %s/T1_to_MNI_nonlin_coeff -o %s/MNI_to_T1_nonlin_field',pref,pref,pref))
    unix(sprintf('applywarp --interp=nn --in=%s/Base_mask.nii.gz --ref=%s/T1_biascorr -w %s/MNI_to_T1_nonlin_field -o %s/T1_biascorr_brain_mask', ...
        pref,pref,pref,pref))
    unix(sprintf('fslmaths %s/T1_biascorr_brain_mask -fillh %s/T1_biascorr_brain_mask',pref,pref))
    unix(sprintf('fslmaths %s/T1_biascorr -mas %s/T1_biascorr_brain_mask %s/T1_biascorr_brain',pref,pref,pref))

    unix(sprintf('flirt -interp nearestneighbour -in %s/T1_biascorr_brain_mask -ref %s/anat.nii.gz -applyxfm -init %s/T1_roi2orig.mat -out %s/anatBrainMask',pref,pref,pref,odir))

    if exist(sprintf('%s/anatBrainMask.nii.gz',odir),'file') 
        if doclean>0
            unix(sprintf('rm -rf %s',pref));
        end
    else
        error('fsl_mask failure to create anatomical mask')
    end


else
    disp('fsl-mask already exists!')
end