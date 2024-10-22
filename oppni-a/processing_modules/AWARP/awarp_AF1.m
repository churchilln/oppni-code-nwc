function awarp_AF1(Adataset, Maskdset, Basedset, odir, ParamCell)
%
% .awarp_AF1:
% .anatomical warping using AFNI utilities
% .adapted from @SSwarper script, uses only latter warping part with predefined mask 

% Adataset = 'inputdataset.nii';  %# req/ input dataset
% SubID    = 'subjID';            %# req/ the subject ID
% Basedset = 'basedset.nii';       %# req/ reference dset- must have 4 bricks
% odir     = 'output';             %# opt/ output dir

% based on modification of afni @SSwarper --> edited to take in arbitrary mask
%                                             subsequent warping done from scratch 

warpscale  = '1';        %# def=1; lower warpscale -> less flexible warps may be useful if odd bumps occur.  vals: [0.1, 1.0] 
minp       = 11;       %# def="11"; the minimum warp patch size
doclean    = 1;       %# def=1; clean out junk files after terminated; can also set =0 to not clean
DO_GIANT   = 0;     %# def=0; use "giant-move" alignment setting, if =1
DO_UNIFIZE = 1; %# def=1; adjustment for intensity nonuniformity
DO_ANISO   = 1; %# def=1; do aniso smoothing processing step
DO_CEIL    = 1;       %# def=1; clipping step - have a ceiling value on anat

cost_aff   = 'lpa+ZZ';   %# used in:  3dAllineate -cost ...
cost_nli   = '-lpa';     %# used in:  3dQwarp, initial rounds
cost_nlf   = '-pcl';     %# used in:  3dQwarp, final rounds
gaus_wt    = '4.5';        %# def for 3dQwarp

% prefix for temp files
pref = [odir,'/__opptmp_p2anat_warp'];

if ~exist( sprintf('%s/anat_warped.nii.gz',odir) ,'file')
    
    % build directory struct recursively
    unix(sprintf('mkdir -p %s',pref));
    % check for path  to base file, create tpath
    if exist(Basedset,'file')
        [tpath,~,~] = fileparts(Basedset);
    else
        error('cannae find the path to the file');
    end
    
    % Require it to have enough bricks to really be a ref
    [~,nvolbase] = unix(sprintf('3dinfo -nv %s',Basedset)); nvolbase = str2num(nvolbase);
    if nvolbase < 4
        error('needs 4 volumes to serve as base!?!?!');
    end
    
    %## start the work
    
    %% ## Step #1: Unifize the input T1
    
    if ~exist(sprintf('%s/anatUAC.nii.gz',odir),'file')
        if DO_UNIFIZE>0
            unix(sprintf('3dUnifize -GM -clfrac 0.4 -Urad 30 -prefix %s/anatU.nii.gz -input %s',odir,Adataset));
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
            unix(sprintf('3danisosmooth -iters 1 -3D -automask -noneg -prefix %s/anatUA.nii.gz %s/anatU.nii.gz',odir,odir));
        else
            unix(sprintf('cp %s/anatU.nii.gz %s/anatUA.nii.gz',odir,odir));
        end
        if DO_CEIL>0
            [~,vvv] = unix(sprintf('3dBrickStat -percentile 98 1 98 -non-zero %s/anatUA.nii.gz',odir)); vvv=str2num(vvv);
            unix(sprintf('3dcalc -a %s/anatUA.nii.gz -expr ''maxbelow(%f,a)'' -prefix %s/anatUAC.nii.gz',odir,vvv(2),odir))
        else
            unix(sprintf('cp %s/anatUA.nii.gz %s/anatUAC.nii.gz',odir,odir));
        end
    else
        disp('warping - UAC''d file already exists, skipping');
    end

    %% ## Step #3: run 3dQwarp first time to a moderate level (skull on)

    if ~exist(sprintf('%s/TAL5.nii.gz',pref),'file')
        if DO_GIANT>0
            unix(sprintf('3dQwarp -echo_edu -lite -base "%s[1]" %s -warpscale %s -source %s/anatUAC.nii.gz -weight "%s[2]" -allineate -noneg -maxlev 5 -iwarp -awarp -wtgaus %s -inedge -workhard:3:5 -nopenalty -prefix "%s/TAL5.nii.gz" -allopt ''-twobest 11 -twopass -maxrot 45 -maxshf 40 -source_automask+2 -cmass''',...
                Basedset,cost_nli,warpscale,odir,Basedset,gaus_wt,pref));
        else
            unix(sprintf('3dQwarp -echo_edu -lite -base "%s[1]" %s -warpscale %s -source %s/anatUAC.nii.gz -weight "%s[2]" -allineate -noneg -maxlev 5 -iwarp -awarp -wtgaus %s -inedge -workhard:3:5 -nopenalty -prefix "%s/TAL5.nii.gz"',...
                Basedset,cost_nli,warpscale,odir,Basedset,gaus_wt,pref));
        end
    end
    
    %# Add in a skull strippin' step --> applied to the UAC'd data
    unix(sprintf('3dcalc -prefix %s/anatSS.nii.gz -a %s/anatUAC.nii.gz -b %s -expr ''a*b''',odir,odir,Maskdset));
    
    %% %# don't do any more if skipping the final (precision) warp
    
    % takes: anatSS.[SubID].nii --> UAC'd anatomical, masked
    %        [pref]TAL5_Allin.aff12.1D --> affine align to template (viz skjull-on warp)
    
    if ~exist(sprintf('%s/AffSS.nii.gz',pref),'file')
        %## Step #6: affine transform that result to template space
        unix(sprintf('3dAllineate -1Dmatrix_apply %s/TAL5_Allin.aff12.1D -source %s/anatSS.nii.gz -master "%s[1]" -cost %s -twopass -autoweight -source_automask -final wsinc5 -prefix %s/AffSS.nii.gz',...
            pref,odir,Basedset,cost_aff,pref)); % MODIFIED=> master is now the refbrick, not the TAL5mm file
    else
        disp('warping - aff-skullstrip found, skipping ahead')
    end
    
    %## warp to template space (skull off), initializing using the previous 3dQwarp -awarp output
    %# Run 3dQwarp in several segments, to avoid gcc OpenMP bug
    %#  where it freezes sometimes with inter-thread conflicts;
    
    if ~exist(sprintf('%s/QQ5.nii.gz',pref),'file')
        %# Piece number 1
        unix(sprintf('3dQwarp -echo_edu -lite -base "%s[0]" -source %s/AffSS.nii.gz -iniwarp %s/TAL5_AWARP.nii.gz -warpscale %s %s -inilev 1 -maxlev 5 -wtgaus %s -inedge -pblur -workhard:5:5 -nodset -prefix %s/QQ5.nii.gz',...
            Basedset,pref,pref,warpscale,cost_nli,gaus_wt,pref));
    else
        disp('warping - nlin piece 1 found, skipping ahead')
    end
    
    if ~exist(sprintf('%s/QQ7.nii.gz',pref),'file')
        %# Piece number 2
        unix(sprintf('3dQwarp -echo_edu -lite -base "%s[0]" -source %s/AffSS.nii.gz -iniwarp %s/QQ5_WARP.nii.gz -warpscale %s %s -inilev 6 -maxlev 7 -workhard:6:7 -wtgaus %s -inedge -pblur -nodset -prefix %s/QQ7.nii.gz',...
            Basedset,pref,pref,warpscale,cost_nlf,gaus_wt,pref));
    else
        disp('warping - nlin piece 2 found, skipping ahead')
    end

    if ~exist(sprintf('%s/anatQQ.nii.gz',pref),'file')
        if minp>13
          %# Final piece for coarse final patch size
          unix(sprintf('3dQwarp -echo_edu -lite -base "%s[0]" -source %s/AffSS.nii.gz -iniwarp %s/QQ7_WARP.nii.gz -warpscale %s %s -inilev 8 -wtgaus %s -inedge -pblur -minpatch %u -Qfinal -prefix %s/anatQQ.nii.gz',...
              Basedset,pref,pref,warpscale,cost_nlf,gaus_wt,minp,odir));
        else
          %# Penultimate piece for refined final patch size
          mpp = minp + 6;
          unix(sprintf('3dQwarp -echo_edu -lite -base "%s[0]" -source %s/AffSS.nii.gz -iniwarp %s/QQ7_WARP.nii.gz -warpscale %s %s -inilev 8 -minpatch %u -wtgaus %s -inedge -pblur -nodset -prefix %s/QQ9.nii.gz',...
              Basedset,pref,pref,warpscale,cost_nlf,mpp,gaus_wt,pref));
        
          %# Ultimate piece for refined final patch size
          unix(sprintf('3dQwarp -echo_edu -lite -base "%s[0]" -source %s/AffSS.nii.gz -iniwarp %s/QQ9_WARP.nii.gz -warpscale %s %s -inilev 10 -wtgaus %s -inedge -pblur -minpatch %u -prefix %s/anatQQ.nii.gz',...
              Basedset,pref,pref,warpscale,cost_nlf,gaus_wt,minp,odir));
        end
    else
        disp('warping - nlin piece 3 found, skipping ahead')
    end
    
    %# DRG's erode-dilate trick for cleanup of little stuff [16 Jan 2019]
    unix(sprintf('3dmask_tool -dilate_inputs -1 1 -prefix %s/de3.nii.gz -input  %s/anatQQ.nii.gz',pref,odir));
    unix(sprintf('3dcalc -a %s/anatQQ.nii.gz -b %s/de3.nii.gz -expr ''a*step(b)'' -prefix %s/anatQQc.nii.gz',odir,pref,odir));
    
    unix(sprintf('rm -f %s/de3.nii.gz anatQQ.nii.gz',pref))
    unix(sprintf('mv -f %s/anatQQc.nii.gz %s/anatQQ.nii.gz',odir,odir))
    
    %# ------------------------------------------------------------------
    %## Clean up the junk (afni-style)
    
    %## Rename affine warp matrix
    
    %# [PT: Feb 19, 2020] Put the mv cmd in an IF condition, because with
    %# re-running extra QC images, might not have the first file name;
    %# now, avoid an unnecessary error msg.
    if exist(sprintf('%s/TAL5_Allin.aff12.1D',pref),'file')
        unix(sprintf('mv -f %s/TAL5_Allin.aff12.1D %s/anatQQ.aff12.1D',pref,odir));
    end
    
    if doclean>0
        unix(sprintf('rm -rf %s',pref));
    end

    %-- 'nother round of cleanup + renaming to fit conventions

    % deleted
    unix(sprintf('rm %s/anatU.nii.gz %s/anatUA.nii.gz',odir,odir))
    % renamed
    unix( sprintf('mv %s/anatUAC.nii.gz %s/anat_proc.nii.gz',odir,odir) ); % processed volume
    unix( sprintf('mv %s/anatSS.nii.gz %s/anat_procss.nii.gz',odir,odir) ); % processed + skull-stripped
    unix( sprintf('mv %s/anatQQ.nii.gz %s/anat_warped.nii.gz',odir,odir) ); % warped volume
    % remask unwarped n warped, for qc etc
    unix(sprintf('3dmask_tool -dilate_input 5 -5 -fill_holes -input %s/anat_warped.nii.gz -prefix %s/anatBrainMask_warped.nii.gz',odir,odir))

else
    disp('afni-warp already exists!')
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
