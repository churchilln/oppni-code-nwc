function aseg_AF1( Anatmasked, Mask, odir, ParamCell )
%
% .aseg_AF1:
% .anatomical tissue segmentation using AFNI utilities
% .implements 3dSeg algorithm

% prefix for temp files
pref = [odir,'/__opptmp_p2anat_aseg'];

if ~exist( sprintf('%s/anat_seg_GM.nii.gz',odir) ,'file') || ...
   ~exist( sprintf('%s/anat_seg_WM.nii.gz',odir) ,'file') || ...
   ~exist( sprintf('%s/anat_seg_CSF.nii.gz',odir) ,'file')

    % build directory struct recursively
    unix(sprintf('mkdir -p %s',pref));

    if contains(Anatmasked,'.nii.gz')
        unix(sprintf('cp %s %s/anatss.nii.gz',Anatmasked,pref));
    elseif contains(Anatmasked,'.nii')
        unix(sprintf('cp %s %s/anatss.nii',Anatmasked,pref));
        unix(sprintf('gzip %s/anatss.nii',pref));
    else
        error('unrecognized Anatmasked format')
    end

    unix(sprintf('3dSeg -anat %s/anatss.nii.gz -mask %s -classes ''CSF ; GM ; WM'' -bias_classes ''GM ; WM'' -bias_fwhm 25 -mixfrac UNI -main_N 5 -blur_meth BFT -prefix %s',...
        pref,Mask,pref))
    unix(sprintf('3dAFNItoNIFTI -prefix %s/Posterior.nii.gz %s/Posterior+orig',pref,pref)); 
    VX = load_untouch_niiz(sprintf('%s/anatss.nii.gz',pref));
    VY = load_untouch_niiz(sprintf('%s/Posterior.nii.gz',pref));

    pvol = double(VY.img);
    pvol = pvol./sum(pvl,3);
    pvol(~isfinite(pvol)) = 0;


    tisslist = {'CSF','GM','WM'}; % tissues in increasing order of T1 intensity
    for i=1:3
        VN = VX;
        VN.img = pvol(:,:,:,i);
        save_untouch_niiz(VN,sprintf('%s/anat_seg_%s.nii.gz', odir, tisslist{i}))
    end

    unix(sprintf('rm -rf %s',pref));
else
    disp('afni-aseg already exists!')
end
