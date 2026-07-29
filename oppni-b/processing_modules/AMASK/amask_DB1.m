function amask_DB1(Adataset, Basedset, odir, ParamCell)
%
% .amask_DB1:
% .anatomical masking using DeepBET

if isempty(ParamCell) || isempty(ParamCell{1})
    doTouchup = false;
else
    doTouchup = ParamCell{1};
    if strcmpi(doTouchup,'True')
        doTouchup=true;
    elseif strcmpi(doTouchup,'False')
        doTouchup=false;
    else
        error('doTouchup option must be True or False!')
    end
end

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

    module_path = fileparts(mfilename('fullpath'));
    helper = fullfile(module_path,'run_deepbet_amask.py');
    if ~exist(helper,'file')
        error('cannot find DeepBET helper');
    end

    [dstat,dmsg] = unix(sprintf('python3 %s --input %s --brain %s/deepbet_brain.nii.gz --mask %s/deepbet_mask.nii.gz',helper,Adataset,pref,pref));
    if dstat~=0
        error('deepbet_mask failure:\n%s',dmsg)
    end

    %--> --> --> Extra "Touchup" step if requested
    if doTouchup
        unix(sprintf('fslmaths %s -mul %s/deepbet_mask %s/anat_midmsk.nii.gz',Adataset,pref,pref))
        unix(sprintf('bet %s/anat_midmsk.nii.gz %s/remsk_bet -f 0.3 -m',pref,pref));
        unix(sprintf('cp %s/remsk_bet_mask.nii.gz %s/anatBrainMask.nii.gz',pref,odir))
    else
        unix(sprintf('cp %s/deepbet_mask.nii.gz %s/anatBrainMask.nii.gz',pref,odir))
    end

    %# ------------------------------------------------------------------
    %## MODIFIED=> Clean up the junk, keeps only anatBrainMask

    if exist(sprintf('%s/anatBrainMask.nii.gz',odir),'file')
        if doclean>0
            unix(sprintf('rm -rf %s',pref));
        end
    else
        error('deepbet_mask failure to create anatomical mask')
    end

else
    disp('deepbet-mask already exists!')
end
