function amask_FS2(Adataset, Basedset, odir, ParamCell)
%
% .amask_FS2:
% .anatomical masking using FreeSurfer SynthStrip

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

    [sstat,smsg] = unix(sprintf('mri_synthstrip -i %s -o %s/synthstrip_brain.nii.gz -m %s/synthstrip_mask.nii.gz',Adataset,pref,pref));
    if sstat~=0
        error('synthstrip_mask failure:\n%s',smsg)
    end

    %--> --> --> Extra "Touchup" step if requested
    if doTouchup
        unix(sprintf('fslmaths %s -mul %s/synthstrip_mask %s/anat_midmsk.nii.gz',Adataset,pref,pref))
        unix(sprintf('bet %s/anat_midmsk.nii.gz %s/remsk_bet -f 0.3 -m',pref,pref));
        unix(sprintf('cp %s/remsk_bet_mask.nii.gz %s/anatBrainMask.nii.gz',pref,odir))
    else
        unix(sprintf('cp %s/synthstrip_mask.nii.gz %s/anatBrainMask.nii.gz',pref,odir))
    end

    %# ------------------------------------------------------------------
    %## MODIFIED=> Clean up the junk, keeps only anatBrainMask

    if exist(sprintf('%s/anatBrainMask.nii.gz',odir),'file')
        if doclean>0
            unix(sprintf('rm -rf %s',pref));
        end
    else
        error('synthstrip_mask failure to create anatomical mask')
    end

else
    disp('synthstrip-mask already exists!')
end
