function amask_AN1(Adataset, Basedset, odir, ParamCell)
%
% .ANTS-based anatomical brain masking

doclean = 1;

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

    unix(sprintf('antsBrainExtraction.sh -d 3 -a %s -e %s/Base_skullon.nii.gz -m %s/Base_mask.nii.gz -o %s/amsk',Adataset, pref, pref,pref));

    unix(sprintf('cp %s/amskBrainExtractionMask.nii.gz %s/anatBrainMask.nii.gz',pref,odir));

    %# ------------------------------------------------------------------
    %## MODIFIED=> Clean up the junk, keeps only anatBrainMask
    
    if exist(sprintf('%s/anatBrainMask.nii.gz',odir),'file') 
        if doclean>0
            unix(sprintf('rm -rf %s',pref));
        end
    else
        error('ants_mask failure to create anatomical mask')
    end

else
    disp('ants-mask already exists!')
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
