function sswarper_script(Adataset, SubID, Basedset, odir)

% set these variables in-code for now

% Adataset = 'inputdataset.nii';  %# req/ input dataset
% SubID    = 'subjID';            %# req/ the subject ID
% Basedset = 'basedset.nii';       %# req/ reference dset- must have 4 bricks
% odir     = 'output';             %# opt/ output dir

warpscale  = '1';        %# def=1; lower warpscale -> less flexible warps may be useful if odd bumps occur.  vals: [0.1, 1.0] 
minp       = 11;       %# def="11"; the minimum warp patch size
doclean    = 1;       %# def=1; clean out junk files after terminated; can also set =0 to not clean
DO_GIANT   = 0;     %# def=0; use "giant-move" alignment setting, if =1
skipwarp   = 0;      %# def=0; just do masking, skip actual warping part
DO_UNIFIZE = 1; %# def=1; adjustment for intensity nonuniformity
DO_SKULLST = 1; %# def=1; do skullstripping before alignment
DO_ANISO   = 1; %# def=1; do aniso smoothing processing step
DO_CEIL    = 1;       %# def=1; clipping step - have a ceiling value on anat
mask_ss    = []; %# def=[]; initial skullstrip mask --> sets DO_SKULLST=0
if ~isempty(mask_ss)
    DO_SKULLST=0;
end

cost_aff   = 'lpa+ZZ';   %# used in:  3dAllineate -cost ...
cost_nli   = '-lpa';     %# used in:  3dQwarp, initial rounds
cost_nlf   = '-pcl';     %# used in:  3dQwarp, final rounds
gaus_wt    = '4.5';        %# def for 3dQwarp

% build directory struct recursively
unix(sprintf('mkdir -p %s',odir));
% prefix for temp files
pref = [odir,'/junk_ssw_'];

% check for path  to base file, create tpath
if exist(Basedset,'file')
    [tpath,~,~] = fileparts(Basedset);
else
    error('cannae find tha path to tha file');
end

% Require it to have enough bricks to really be a ref
[~,nvolbase] = unix(sprintf('3dinfo -nv %s',Basedset)); nvolbase = str2num(nvolbase);
if nvolbase < 5
    error('needs 5 volumes to serve as base!?!?!');
end

%## Step #-1: copy the raw anat (and perhaps mask_ss) into output
%## dir---often useful to compare later and/or check header info

dset_cp = [ odir, '/anat_cp.',SubID,'.nii'];
unix(sprintf('cp %s %s',Adataset,dset_cp));
if ~isempty(mask_ss)
    mask_ss_cp = [odir, '/mask_ss_cp.',SubID,'.nii'];
    unix(sprintf('cp %s %s',mask_ss,mask_ss_cp));
end

%## start the work

%% ## Step #1: Unifize the input T1

if DO_UNIFIZE>0
    unix(sprintf('3dUnifize -GM -clfrac 0.4 -Urad 30 -prefix %s/anatU.%s.nii -input %s',odir,SubID,Adataset));
else
    unix(sprintf('cp %s %s/anatU.%s.nii',Adataset,odir,SubID));
end

if DO_ANISO>0
    unix(sprintf('3danisosmooth -iters 1 -3D -automask -noneg -prefix %s/anatUA.%s.nii %s/anatU.%s.nii',odir,SubID,odir,SubID));
else
    unix(sprintf('cp %s/anatU.%s.nii %s/anatUA.%s.nii',odir,SubI,odir,SubID));
end

if DO_CEIL>0
    [~,vvv] = unix(sprintf('3dBrickStat -percentile 98 1 98 -non-zero %s/anatUA.%s.nii',odir,SubID)); vvv=str2num(vvv);
    unix(sprintf('3dcalc -a %s/anatUA.%s.nii -expr ''maxbelow(%f,a)'' -prefix %s/anatUAC.%s.nii',odir,SubID,vvv(2),odir,SubID))
else
    unix(sprintf('cp %s/anatUA.%s.nii %s/anatUAC.%s.nii',odir,SubI,odir,SubID));
end

%% ## Step #2: Strip Skull (Ziad's way)
if DO_SKULLST>0 
    unix(sprintf('3dSkullStrip -input %s/anatUAC.%s.nii -prefix %s/anatS.%s.nii -debug 1 -ld 33 -niter 777 -shrink_fac_bot_lim 0.777 -exp_frac 0.0666 -orig_vol',odir,SubID,odir,SubID));
elseif ~isempty(mask_ss)
    unix(sprintf('3dcalc -a %s/anatUAC.%s.nii -b %s -expr ''a*step(b)'' -prefix %s/anatS.%s.nii',odir,SubID,mask_ss,odir,SubID));
else
    unix(sprintf('cp %s/anatUAC.%s.nii %s/anatS.%s.nii',odir,SubI,odir,SubID));
end

%% ## Step #3: run 3dQwarp first time to a moderate level (skull on)
if DO_GIANT>0
    unix(sprintf('3dQwarp -echo_edu -lite -base "%s[1]" %s -warpscale %s -source %s/anatUAC.%s.nii -weight "%s[2]" -allineate -noneg -maxlev 5 -iwarp -awarp -wtgaus %s -inedge -workhard:3:5 -nopenalty -prefix "%sTAL5.nii" -allopt ''-twobest 11 -twopass -maxrot 45 -maxshf 40 -source_automask+2 -cmass''',...
        Basedset,cost_nli,warpscale,odir,SubID,Basedset,gaus_wt,pref));
else
    unix(sprintf('3dQwarp -echo_edu -lite -base "%s[1]" %s -warpscale %s -source %s/anatUAC.%s.nii -weight "%s[2]" -allineate -noneg -maxlev 5 -iwarp -awarp -wtgaus %s -inedge -workhard:3:5 -nopenalty -prefix "%sTAL5.nii"',...
        Basedset,cost_nli,warpscale,odir,SubID,Basedset,gaus_wt,pref));
end

%# -------- check early NL alignment to ref vol space ---------> (PAS INLCUS)

%## Step #4: mask off the skull using the template (second skull-strip)
unix(sprintf('3dmask_tool -input "%s[3]" -dilate_input 2 -prefix %sMASK.nii',Basedset,pref));
unix(sprintf('3dcalc -a %sMASK.nii -b %sTAL5.nii -expr ''step(a)*b'' -prefix %sTAL5mm.nii',pref,pref,pref));

%## Step #5: warp this masked dataset back to original space
unix(sprintf('3dNwarpApply -echo_edu -nwarp %sTAL5_WARPINV.nii -master %s/anatS.%s.nii -source %sTAL5mm.nii -prefix %sTAL5ww.nii',pref,odir,SubID,pref,pref));

%## warp the mask itself (dilated) back to orig space
unix(sprintf('rm -f %sMASK.nii',pref));
unix(sprintf('3dmask_tool -echo_edu -input "%s[3]" -dilate_input 3 -prefix %sMASK.nii',Basedset,pref));
unix(sprintf('3dNwarpApply -nwarp %sTAL5_WARPINV.nii -master %s/anatS.%s.nii -source %sMASK.nii -prefix %sMASKO.nii -ainterp NN',pref,odir,SubID,pref,pref));

%## merge these backward warped datasets with the 3dSkullStrip output to get a better original skull-stripped result
unix(sprintf('3dcalc -a %s/anatS.%s.nii -b %sTAL5ww.nii -c %sMASKO.nii -expr ''step(c)*max(a,b)'' -prefix %s/anatSS.%s.nii',odir,SubID,pref,pref,odir,SubID));

%# DRG's erode-dilate trick for cleanup of little crap [16 Jan 2019]
unix(sprintf('3dmask_tool -dilate_inputs -3 3 -prefix %sde3.nii -input  %s/anatSS.%s.nii',pref,odir,SubID));
unix(sprintf('3dcalc -a %s/anatSS.%s.nii -b %sde3.nii -expr ''a*step(b)'' -prefix %s/anatSSc.%s.nii',odir,SubID,pref,odir,SubID));

%# Throw in an automask step to clean up the outer edges [18 Jun 2019]
unix(sprintf('rm -f %sde3.nii %s/anatSS.%s.nii',pref,odir,SubID));
unix(sprintf('3dAutomask -apply_prefix %s/anatSSd.%s.nii %s/anatSSc.%s.nii',odir,SubID,odir,SubID));
unix(sprintf('mv -f %s/anatSSd.%s.nii %s/anatSS.%s.nii',odir,SubID,odir,SubID));
unix(sprintf('rm -f %s/anatSSc.%s.nii',odir,SubID));

%% %# don't do any more if skipping the final (precision) warp

if skipwarp==0
    % takes: anatSS.[SubID].nii --> UAC'd anatomical, masked
    %        [pref]TAL5_Allin.aff12.1D --> affine align to template (viz skjull-on warp)

    %## Step #6: affine transform that result to template space
    unix(sprintf('3dAllineate -1Dmatrix_apply %sTAL5_Allin.aff12.1D -source %s/anatSS.%s.nii -master %sTAL5mm.nii -cost %s -twopass -autoweight -source_automask -final wsinc5 -prefix %sAffSS.nii',...
        pref,odir,SubID,pref,cost_aff,pref));

    %## warp to template space (skull off), initializing using the previous 3dQwarp -awarp output
    %# Run 3dQwarp in several segments, to avoid gcc OpenMP bug
    %#  where it freezes sometimes with inter-thread conflicts;

    %# Piece number 1
    unix(sprintf('3dQwarp -echo_edu -lite -base "%s[0]" -source %sAffSS.nii -iniwarp %sTAL5_AWARP.nii -warpscale %s %s -inilev 1 -maxlev 5 -wtgaus %s -inedge -pblur -workhard:5:5 -nodset -prefix %sQQ5.nii',...
        Basedset,pref,pref,warpscale,cost_nli,gaus_wt,pref));

    %# Piece number 2
    unix(sprintf('3dQwarp -echo_edu -lite -base "%s[0]" -source %sAffSS.nii -iniwarp %sQQ5_WARP.nii -warpscale %s %s -inilev 6 -maxlev 7 -workhard:6:7 -wtgaus %s -inedge -pblur -nodset -prefix %sQQ7.nii',...
        Basedset,pref,pref,warpscale,cost_nlf,gaus_wt,pref));

    if minp>13
      %# Final piece for coarse final patch size
      unix(sprintf('3dQwarp -echo_edu -lite -base "%s[0]" -source %sAffSS.nii -iniwarp %sQQ7_WARP.nii -warpscale %s %s -inilev 8 -wtgaus %s -inedge -pblur -minpatch %u -Qfinal -prefix %s/anatQQ.%s.nii',...
          Basedset,pref,pref,warpscale,cost_nlf,gaus_wt,minp,odir,SubID));
    else
      %# Penultimate piece for refined final patch size
      mpp = minp + 6;
      unix(sprintf('3dQwarp -echo_edu -lite -base "%s[0]" -source %sAffSS.nii -iniwarp %sQQ7_WARP.nii -warpscale %s %s -inilev 8 -minpatch %u -wtgaus %s -inedge -pblur -nodset -prefix %sQQ9.nii',...
          Basedset,pref,pref,warpscale,cost_nlf,mpp,gaus_wt,pref));

      %# Ultimate piece for refined final patch size
      unix(sprintf('3dQwarp -echo_edu -lite -base "%s[0]" -source %sAffSS.nii -iniwarp %sQQ9_WARP.nii -warpscale %s %s -inilev 10 -wtgaus %s -inedge -pblur -minpatch %u -prefix %s/anatQQ.%s.nii',...
          Basedset,pref,pref,warpscale,cost_nlf,gaus_wt,minp,odir,SubID));
    end

    %# DRG's erode-dilate trick for cleanup of little stuff [16 Jan 2019]
    unix(sprintf('3dmask_tool -dilate_inputs -1 1 -prefix %sde3.nii -input  %s/anatQQ.%s.nii',pref,odir,SubID));
    unix(sprintf('3dcalc -a %s/anatQQ.%s.nii -b %sde3.nii -expr ''a*step(b)'' -prefix %s/anatQQc.%s.nii',odir,SubID,pref,odir,SubID));

    unix(sprintf('rm -f %sde3.nii anatQQ.%s.nii',pref,SubID))
    unix(sprintf('mv -f %s/anatQQc.%s.nii %s/anatQQ.%s.nii',odir,SubID,odir,SubID))
end

%# ------------------------------------------------------------------
%## Clean up the junk

%## Rename affine warp matrix

%# [PT: Feb 19, 2020] Put the mv cmd in an IF condition, because with
%# re-running extra QC images, might not have the first file name;
%# now, avoid an unnecessary error msg.
if exist(sprintf('%sTAL5_Allin.aff12.1D',pref),'file')
    unix(sprintf('mv -f %sTAL5_Allin.aff12.1D %s/anatQQ.%s.aff12.1D',pref,odir,SubID));
end

if doclean>0
    unix(sprintf('rm -f %s*',pref));
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
