function aseg_AN1( Anatmasked, Mask, odir, ParamCell )
%
% .aseg_AN1:
% .anatomical tissue segmentation using ANTs utilities
% .implements antsAtroposN4.sh script

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

    unix(sprintf('antsAtroposN4.sh -d 3 -a %s/anatss.nii.gz -x %s -c 3 -o %s/x_',pref,Mask,pref));

    % renaming 'em
    VX = load_untouch_niiz(sprintf('%s/anatss.nii.gz',pref));
    for i=1:3
        VS = load_untouch_niiz(sprintf('%s/x_SegmentationPosteriors%u.nii.gz',pref,i));
        sval(i,1) = mean(double(VX.img(VS.img>0.5)));
    end
    isort = sortrows([(1:3)',sval],2,'ascend');
    isort = isort(:,1);           % pves indexed by increasing mean intensity
    tisslist = {'CSF','GM','WM'}; % tissues in increasing order of T1 intensity
    for i=1:3
        unix(sprintf('mv %s/x_SegmentationPosteriors%u.nii.gz %s/anat_seg_%s.nii.gz',pref, isort(i), odir, tisslist{i}));
    end

    unix(sprintf('rm -rf %s',pref));
else
    disp('ants-aseg already exists!')
end
