function amask_AF1(Adataset, Basedset, odir, ParamCell)
%
% .amask_AF1:
% .anatomical masking using AFNI utilities
% .adapted from @SSwarper script, uses only initial warp+mask part

if isempty(ParamCell) || isempty(ParamCell{1})
    doTouchup = False;
else
    doTouchup = ParamCell{1};
    if strcmpi(doTouchup,'True')
        doTouchup=True;
    elseif strcmpi(doTouchup,'False')
        doTouchup=False;
    else
        error('doTouchup option must be True or False!')
    end
end

% Adataset = 'inputdataset.nii';  %# req/ input dataset
% SubID    = 'subjID';            %# req/ the subject ID
% Basedset = 'basedset.nii';       %# req/ reference dset- must have 4 bricks
% odir     = 'output';             %# opt/ output dir

warpscale  = '1';        %# def=1; lower warpscale -> less flexible warps may be useful if odd bumps occur.  vals: [0.1, 1.0] 
doclean    = 1;       %# def=1; clean out junk files after terminated; can also set =0 to not clean
DO_GIANT   = 0;     %# def=0; use "giant-move" alignment setting, if =1
DO_UNIFIZE = 1; %# def=1; adjustment for intensity nonuniformity
DO_SKULLST = 1; %# def=1; do skullstripping before alignment
DO_ANISO   = 1; %# def=1; do aniso smoothing processing step
DO_CEIL    = 1;       %# def=1; clipping step - have a ceiling value on anat

cost_nli   = '-lpa';     %# used in:  3dQwarp, initial rounds
gaus_wt    = '4.5';        %# def for 3dQwarp

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
    if nvolbase < 5
        error('needs 5 volumes to serve as base!?!?!');
    end
    
    %## start the work. MODIFIED=> all outputs now go to "pref", not just a select few
    
    %% ## Step #1: Unifize the input T1
    
    if ~exist(sprintf('%s/anatUAC.nii.gz',pref),'file')
        if DO_UNIFIZE>0
            unix(sprintf('3dUnifize -GM -clfrac 0.4 -Urad 30 -prefix %s/anatU.nii.gz -input %s',pref,Adataset));
        else
            if contains(Adataset,'.nii.gz')
                unix(sprintf('cp %s %s/anatU.nii.gz',Adataset,pref));
            elseif contains(Adataset,'.nii')
                unix(sprintf('cp %s %s/anatU.nii',Adataset,pref));
                unix(sprintf('gzip %s/anatU.nii',pref));
            else
                error('unregocnized Adataset format');
            end
        end
        if DO_ANISO>0
            unix(sprintf('3danisosmooth -iters 1 -3D -automask -noneg -prefix %s/anatUA.nii.gz %s/anatU.nii.gz',pref,pref));
        else
            unix(sprintf('cp %s/anatU.nii.gz %s/anatUA.nii.gz',pref,pref));
        end
        if DO_CEIL>0
            [~,vvv] = unix(sprintf('3dBrickStat -percentile 98 1 98 -non-zero %s/anatUA.nii.gz',pref)); vvv=str2num(vvv);
            unix(sprintf('3dcalc -a %s/anatUA.nii.gz -expr ''maxbelow(%f,a)'' -prefix %s/anatUAC.nii.gz',pref,vvv(2),pref))
        else
            unix(sprintf('cp %s/anatUA.nii.gz %s/anatUAC.nii.gz',pref,pref));
        end
    else
        disp('masking - UAC''d file already exists, skipping');
    end

    %% ## Step #2: Strip Skull (Ziad's way)
    if ~exist(sprintf('%s/anatS.nii.gz',pref),'file')
        if DO_SKULLST>0 
            unix(sprintf('3dSkullStrip -input %s/anatUAC.nii.gz -prefix %s/anatS.nii.gz -debug 1 -ld 33 -niter 777 -shrink_fac_bot_lim 0.777 -exp_frac 0.0666 -orig_vol',pref,pref));
        else
            unix(sprintf('cp %s/anatUAC.nii.gz %s/anatS.nii.gz',pref,pref));
        end
    else
        disp('masking - first-strip file already exists, skipping');
    end
    
    %% ## Step #3: run 3dQwarp first time to a moderate level (skull on)
    if ~exist(sprintf('%s/TAL5.nii.gz',pref),'file')
        if DO_GIANT>0
            unix(sprintf('3dQwarp -echo_edu -lite -base "%s[1]" %s -warpscale %s -source %s/anatUAC.nii.gz -weight "%s[2]" -allineate -noneg -maxlev 5 -iwarp -awarp -wtgaus %s -inedge -workhard:3:5 -nopenalty -prefix "%s/TAL5.nii.gz" -allopt ''-twobest 11 -twopass -maxrot 45 -maxshf 40 -source_automask+2 -cmass''',...
                Basedset,cost_nli,warpscale,pref,Basedset,gaus_wt,pref));
        else
            unix(sprintf('3dQwarp -echo_edu -lite -base "%s[1]" %s -warpscale %s -source %s/anatUAC.nii.gz -weight "%s[2]" -allineate -noneg -maxlev 5 -iwarp -awarp -wtgaus %s -inedge -workhard:3:5 -nopenalty -prefix "%s/TAL5.nii.gz"',...
                Basedset,cost_nli,warpscale,pref,Basedset,gaus_wt,pref));
        end
    else
        disp('masking - skull-on warp already exists, skipping')
    end
    %# -------- check early NL alignment to ref vol space ---------> (PAS INLCUS)
    
    if ~exist(sprintf('%s/anatSS.nii.gz',pref),'file')

        %# clearing any intermediate products
        unix(sprintf('rm %s/anatSSd.nii.gz %s/anatSSd.nii.gz %s/de3.nii.gz %s/MASK.nii.gz %s/MASKO.nii.gz %s/TAL5mm.nii.gz %s/TAL5ww.nii.gz',pref,pref,pref,pref,pref,pref,pref))

        %## Step #4: mask off the skull using the template (second skull-strip)
        unix(sprintf('3dmask_tool -input "%s[3]" -dilate_input 2 -prefix %s/MASK.nii.gz',Basedset,pref));
        unix(sprintf('3dcalc -a %s/MASK.nii.gz -b %s/TAL5.nii.gz -expr ''step(a)*b'' -prefix %s/TAL5mm.nii.gz',pref,pref,pref));
        
        %## Step #5: warp this masked dataset back to original space
        unix(sprintf('3dNwarpApply -echo_edu -nwarp %s/TAL5_WARPINV.nii.gz -master %s/anatS.nii.gz -source %s/TAL5mm.nii.gz -prefix %s/TAL5ww.nii.gz',pref,pref,pref,pref));
        
        %## warp the mask itself (dilated) back to orig space
        unix(sprintf('rm -f %s/MASK.nii.gz',pref));
        unix(sprintf('3dmask_tool -echo_edu -input "%s[3]" -dilate_input 3 -prefix %s/MASK.nii.gz',Basedset,pref));
        unix(sprintf('3dNwarpApply -nwarp %s/TAL5_WARPINV.nii.gz -master %s/anatS.nii.gz -source %s/MASK.nii.gz -prefix %s/MASKO.nii.gz -ainterp NN',pref,pref,pref,pref));
        
        %## merge these backward warped datasets with the 3dSkullStrip output to get a better original skull-stripped result
        unix(sprintf('3dcalc -a %s/anatS.nii.gz -b %s/TAL5ww.nii.gz -c %s/MASKO.nii.gz -expr ''step(c)*max(a,b)'' -prefix %s/anatSS.nii.gz',pref,pref,pref,pref));
        
        %# DRG's erode-dilate trick for cleanup of little crap [16 Jan 2019]
        unix(sprintf('3dmask_tool -dilate_inputs -3 3 -prefix %s/de3.nii.gz -input  %s/anatSS.nii.gz',pref,pref));
        unix(sprintf('3dcalc -a %s/anatSS.nii.gz -b %s/de3.nii.gz -expr ''a*step(b)'' -prefix %s/anatSSc.nii.gz',pref,pref,pref));
        
        %# Throw in an automask step to clean up the outer edges [18 Jun 2019]
        unix(sprintf('rm -f %s/de3.nii.gz %s/anatSS.nii.gz',pref,pref));
        unix(sprintf('3dAutomask -apply_prefix %s/anatSSd.nii.gz -prefix %s/anatBrainMask_unt.nii.gz  %s/anatSSc.nii.gz',pref,pref,pref)); % MODIFIED=> create anatBrainMask.nii!
        unix(sprintf('mv -f %s/anatSSd.nii.gz %s/anatSS.nii.gz',pref,pref));
        unix(sprintf('rm -f %s/anatSSc.nii.gz',pref));

        %--> --> --> Extra "Touchup" step if requested
        if doTouchup
            unix(sprintf('fslmaths %s -mul %s/anatBrainMask_unt %s/anat_midmsk.nii.gz',Adataset,pref,pref))
            unix(sprintf('bet %s/anat_midmsk.nii.gz %s/remsk_bet -f 0.3 -m',pref,pref));
            unix(sprintf('cp %s/remsk_bet_mask.nii.gz %s/anatBrainMask.nii.gz',pref,odir))
        else
            unix(sprintf('cp %s/anatBrainMask_unt.nii.gz %s/anatBrainMask.nii.gz',pref,odir))
        end

    else
        disp('masking - tidied up mask already exists?? skipping??')
    end
    
    %# ------------------------------------------------------------------
    %## MODIFIED=> Clean up the junk, keeps only anatBrainMask
    
    if exist(sprintf('%s/anatBrainMask.nii.gz',odir),'file') 
        if doclean>0
            unix(sprintf('rm -rf %s',pref));
        end
    else
        error('afni_mask failure to create anatomical mask')
    end

else
    disp('afni-mask already exists!')
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
