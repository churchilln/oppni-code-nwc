function P2_perf_dataProcessing( inputfile, pipelinefile, paramlist, outpath, big_skip, mask_subj_idxes, subj_subset )
%
% =========================================================================
% P2_PERF_DATAPROCESSING: this script should be run after the "P0" and "P1"
% steps of the bold pipeline. It does the actual data processing!
% =========================================================================
%
% Syntax:
%
%     P2_perf_dataProcessing( inputfile, pipelinefile, paramlist, outpath, big_skip, mask_subj_idxes )
%
% Input:
%      inputfile : string, giving name of input textfile listing data to process
%      pipelinefiles : string, giving name of pipeline textfile specifying which pipeline steps to apply to the data
%      paramlist : string, giving name of parameter file specifying special arguments to use in pipelines
%      outpath : string, specifying destination directory for processed outputs
%      big_skip : completely bypasses anatomical processing checks / most of functional processing
%                 0=do not skip, 1=do skip. ONLY USE IF YOU KNOW WHAT YOU ARE DOING!
%      mask_subj_idxes : 
%     
%          if mask_subj_idxes=[], it will take all subjects in your input file
%          if mask_subj_idxes=numeric vector, this is will pull corresponding subjects from lines of the input file
%          if mask_subj_idxes=matfile, it will load the file and search for a "mask_subj_idxes" numeric vector, treated as above
%                                      if no such field, will instead seek "importpath" structure
%                                      organized as:
%                                         importpath.brain_maps
%                                         importpath.masks
%                                         importpath.parcellations
%                                      *this will set to IMPORT_MASKFILES --> copies over contents of brain_maps/masks/parcellations
%

% declaring path
CODE_PATH = fileparts(which('P0_perf_populateDirectories.m'));
if CODE_PATH(end)~='/'
    CODE_PATH = [CODE_PATH '/'];
end

if nargin<5
    big_skip=0;
end
if nargin<6
    mask_subj_idxes=[];
end
if nargin<7
    subj_subset=0;
end

% initializing structure and running checkes
[subject_list, InputStruct_aug, PipeStruct_aug, ParamStruct_aug] = P0_perf_populateDirectories( inputfile, pipelinefile, paramlist, outpath );

% now augmenting outpath... do this after P0!
outpath = fullfile(outpath,'perf_proc'); % subdir should be fmri_proc
% check for missing files from previous step too
File_Existence_Checker_perf(InputStruct_aug,outpath,1); 

% check validity of perfusion model --> carry forward
perf_model_struct = check_perf_model( ParamStruct_aug.PERF_MODEL );
% check validity of analysis model etc. --> carry forward information about the analysis too!
analysis_struct = check_perf_analysis_model( ParamStruct_aug.ANALYSIS );
% also, modify "number of components" field if more than one contrast is specified
if( ~isempty(strfind(ParamStruct_aug.CONTRAST,',' )) )
    analysis_struct.num_comp = 'multi_component';
end

% list of subjects for constructing group masks
IMPORT_MASKFILES=0;
if isempty(mask_subj_idxes)
    disp('using all subj in current input file for mask construction!')
    mask_subj_idxes = 1:numel(subject_list);
elseif ischar(mask_subj_idxes) && exist(mask_subj_idxes,'file')
    
    x = load(mask_subj_idxes);
    if isfield(x,'mask_subj_idxes')
        disp('loading list of subject rows in input file for mask construction!');
        if isempty(x.mask_subj_idxes)
            error('list of subject rows for mask-making is empty')
        elseif numel(x.mask_subj_idxes)<10
            warning('not a lot of subjects for group mask-making ... might be unstable')
        end
        mask_subj_idxes = x.mask_subj_idxes; clear x;
    elseif isfield(x,'importpath')
        IMPORT_MASKFILES=1; % case where we import files
        mask_subj_paths.brain_maps    = x.importpath.brain_maps;
        mask_subj_paths.masks         = x.importpath.masks;
        mask_subj_paths.parcellations = x.importpath.parcellations;
        clear mask_subj_idxes x;
    else
        error('mask_subj_idxes matfile contains information in unrecognized format...')
    end
elseif ~isnumeric(mask_subj_idxes)
    error('unrecognized mask id format?')
else
    disp('using numeric list of subj values for mask construction!')
end

if IMPORT_MASKFILES==0
    % store information about masking sublist...
    subject_list_formask = subject_list(mask_subj_idxes);
    % consistency checking
    if exist([outpath,'/_group_level/pipe_',PipeStruct_aug.PNAME{1},'_mask_subj_idxes.mat'],'file')
        x=load([outpath,'/_group_level/pipe_',PipeStruct_aug.PNAME{1},'_mask_subj_idxes.mat']);
        if     ~isempty( setdiff(subject_list_formask,x.subject_list_formask) ) 
            error('custom subject list for masking :: subjects in new list not present in old! delete group level folders if you want to update!')
        elseif ~isempty( setdiff(x.subject_list_formask,subject_list_formask) )
            error('custom subject list for masking :: subjects not in new list that are present in old! delete group level folders if you want to update!')
        else
            disp('custom subject list for masking :: list is consistent with old one ... continuing without modification!')
        end        
    end
    % re-saving, in case indexes need updating ( they may change, as long as subject prefix list doesnt )
    save([outpath,'/_group_level/pipe_',PipeStruct_aug.PNAME{1},'_mask_subj_idxes.mat'],'mask_subj_idxes','subject_list_formask');
else
    % store information about importing paths + consistency checking...
    if exist([outpath,'/_group_level/pipe_',PipeStruct_aug.PNAME{1},'_mask_subj_paths.mat'],'file')
        x=load([outpath,'/_group_level/pipe_',PipeStruct_aug.PNAME{1},'_mask_subj_paths.mat']);
        mismatchpth=0; mistag=[];
        if ~strcmp(mask_subj_paths.brain_maps,x.mask_subj_paths.brain_maps)
            mismatchpth=mismatchpth+1;
            mistag = [mistag, ', brain_maps'];
        elseif ~strcmp(mask_subj_paths.masks,x.mask_subj_paths.masks)
            mismatchpth=mismatchpth+1;
            mistag = [mistag, ', brain_maps'];
        elseif ~strcmp(mask_subj_paths.parcellations,x.mask_subj_paths.parcellations)
            mismatchpth=mismatchpth+1;
            mistag = [mistag, ', brain_maps'];
        end
        if mismatchpth>0
            error('custom imported files for masking :: %s path(s) differ from old! delete group level folders if you want to update!', mistag(3:end));
        else
            disp('custom imported files for masking :: paths are consistent with old ones ... continuing without modification!')
        end
    else
        save([outpath,'/_group_level/pipe_',PipeStruct_aug.PNAME{1},'_mask_subj_paths.mat'],'mask_subj_paths');
    end
    % --> import everything immediately!
    unix(sprintf('cp %s/*.nii* %s/_group_level/brain_maps/pipe_%s',mask_subj_paths.brain_maps,outpath,PipeStruct_aug.PNAME{1}))
    unix(sprintf('cp %s/*.nii* %s/_group_level/masks/pipe_%s',mask_subj_paths.masks,outpath,PipeStruct_aug.PNAME{1}))
    unix(sprintf('cp %s/*.nii* %s/_group_level/parcellations/pipe_%s',mask_subj_paths.parcellations,outpath,PipeStruct_aug.PNAME{1}))
end
% .to reset this part, need to delete group_level folder + strip out P2 processed data! 

% subject listing
if isempty(subj_subset) || subj_subset==0 || subj_subset<0
    subj_list_for_proc = 1:numel(subject_list);
else
    subj_list_for_proc = subj_subset;
end

for ns=subj_list_for_proc % step through anat-proc, func-proc (block-1)

    % check existence of subject specific struct file
    if exist(fullfile( outpath,subject_list{ns},'InputStruct_ssa.mat'),'file')  
        load( fullfile( outpath,subject_list{ns},'InputStruct_ssa.mat') )
    else
        error('cannot find Input struct file for subject: %s \n');
    end

    % quick formatting stuff, again assuRImes that directory structure was already constructed in "P0" pipeline step 
    opath0 = fullfile(outpath,InputStruct_ssa.PREFIX,'rawdata');
    %
    opath1a = fullfile( outpath, InputStruct_ssa.PREFIX,'anat_proc');
    opath2a = fullfile( opath1a,sprintf('subpipe_%03u',    PipeStruct_aug.pipe_idx.Warp));
    opath3a = fullfile( opath2a,sprintf('seg_subpipe_%03u',PipeStruct_aug.pipe_idx.Seg ));
    %
    opath1p = fullfile(outpath,InputStruct_ssa.PREFIX,'phys_proc');
    %
    opath1f = fullfile(outpath,InputStruct_ssa.PREFIX,'perf_proc_p1');
    opath2f = fullfile(opath1f,sprintf('subpipe_%03u',PipeStruct_aug.pipe_idx.P1));
    opath3f = fullfile(opath2f,sprintf('seg_subpipe_%03u',PipeStruct_aug.pipe_idx.Seg ));
    %
    opath4f = fullfile(outpath,InputStruct_ssa.PREFIX,'perf_proc_p2',['pipe_',PipeStruct_aug.PNAME{1}]); 

    fprintf('\n===> subj %s (%u/%u), anat. processing...\n',subject_list{ns},ns,numel(subject_list)),
    %% =======================================================================
    %%      ANAT Processing ...
    %% =======================================================================
    
    nr=1; %-fixing to run=1 for now

    % tag for rawdata in case z-clipping was performed
    if (ischar(InputStruct_ssa.arun(nr).ZCLIP_thr) && strcmpi(InputStruct_ssa.arun(nr).ZCLIP_thr,'AUTO')) || (isnumeric(InputStruct_ssa.arun(nr).ZCLIP_thr) && isfinite(InputStruct_ssa.arun(nr).ZCLIP_thr))
        zclip_tag = '_zclip';
    else
        zclip_tag = '';
    end

    % fixing orientation stuff -- adjusting for obliquity, switching to standard mni-compatible orientation 
    if ~exist(sprintf('%s/anat%u_2std.nii.gz',opath1a,nr),'file')
        unix(sprintf('3dWarp -oblique2card -prefix %s/anat%u_deob.nii.gz -cubic %s/anat%u%s.nii.gz', opath1a,nr, opath0,nr, zclip_tag)); %-wsinc5
        unix(sprintf('fslreorient2std %s/anat%u_deob.nii.gz %s/anat%u_2std.nii.gz', opath1a,nr, opath1a,nr));        
        if exist(sprintf('%s/anat%u_2std.nii.gz',opath1a,nr),'file')
            unix(sprintf('rm %s/anat%u_deob.nii.gz', opath1a,nr)); % not a really useful intermediate - delete it
        else
            error('failed to create deob/reoriented anatomical file')
        end
    end

    % >> Pipeline Step #1: "AMASK" 
    % >> produces a mask of brain voxels within the T1 image
    % >> takes in raw, reoriented anatomical, outputs binary image anatBrainMask.nii
    if strcmpi(PipeStruct_aug.AMASK{1},'AF0')
        %do 3dskullstrip...
        error('afni-0 currently unavailable')
    elseif strcmpi(PipeStruct_aug.AMASK{1},'AF1')
        amask_AF1( sprintf('%s/anat%u_2std.nii.gz', opath1a,nr), ParamStruct_aug.TEMPLATE_loc, opath2a);
    elseif strcmpi(PipeStruct_aug.AMASK{1},'FL0')
        %do bet...
        error('fsl-0 currently unavailable')
    elseif strcmpi(PipeStruct_aug.AMASK{1},'FL1')
        %do glasser bet (script)
        error('fsl-1 currently unavailable')
    elseif strcmpi(PipeStruct_aug.AMASK{1},'AN1')
        amask_AN1( sprintf('%s/anat%u_2std.nii.gz', opath1a,nr), ParamStruct_aug.TEMPLATE_loc, opath2a);
    elseif strcmpi(PipeStruct_aug.AMASK{1},'CST')
        % just check for custom mask
        if ~exist( sprintf('%s/anatBrainMask.nii.gz',opath2a), 'file')
            error('a custom "anatBrainMask.nii.gz" was expected!')
        end
    else
        error('unrecognized mask optione!')
    end

    % >> Pipeline Step #2: "AWARP"
    % >> creates the warp transforming anatomical images to a template.
    % >> processes & transforms anatomical scans; combined with "ASEG" produces group tissue maps, and combined with "FWARP" transforms the functional scans
    % >> takes in raw, reoriented anatomical, produces "anat_proc.nii" (tidied up scan), "anat_procss.nii (skull-stripped version), "anat_warped.nii" and "anatBrainMask_warped.nii"
    % >> also produces requisite warpfiles in output dir. These specific outputs depend on the software package being used 
    if strcmpi(PipeStruct_aug.AWARP{1},'AF1')
        % do afni-style warping
        awarp_AF1( sprintf('%s/anat%u_2std.nii.gz', opath1a,nr), sprintf('%s/anatBrainMask.nii.gz',opath2a), ParamStruct_aug.TEMPLATE_loc, opath2a);
    elseif strcmpi(PipeStruct_aug.AWARP{1},'AN1')
        %%% *** Need fixing for compatibility with param struct
        error('ANTS warping not yet ready to go!')
        % creating bias-corrected, skull-stripped anats
%         unix(sprintf('antsBrainExtraction.sh -d 3 -a %s/anat%u_2std.nii -e MNI152_2009_template_SSW_0.nii -m MNI152_2009_template_SSW_mask.nii -o %s/amsk -k 1',opath1a,nr, opath2a));
%         % params from "newAntsExample.sh"+
%         fixd = 'MNI152_2009_template_SSW_brain.nii';
%         movn = sprintf('%s/amskBrainExtractionBrain.nii.gz',opath2a);
%         nomen = sprintf('%s/aln',opath2a);
%         unix(sprintf(['antsRegistration ', ...
%         '-d 3 ', ...
%         '-r [ %s, %s , 1 ] ', ...
%         '-m mattes[  %s, %s , 1 , 32, regular, 0.3 ] ', ...
%         '-t translation[ 0.1 ] ', ...
%         '-c [ 10000x111110x11110,1.e-8,20 ] ', ...
%         '-s 4x2x1vox  ', ...
%         '-f 6x4x2 -l 1 ', ...
%         '-m mattes[  %s, %s , 1 , 32, regular, 0.3 ] ', ...
%         '-t rigid[ 0.1 ] ', ...
%         '-c [ 10000x111110x11110,1.e-8,20 ] ', ...
%         '-s 4x2x1vox  ', ...
%         '-f 3x2x1 -l 1 ', ...
%         '-m mattes[  %s, %s , 1 , 32, regular, 0.3 ] ', ...
%         '-t affine[ 0.1 ] ', ...
%         '-c [ 10000x111110x11110,1.e-8,20 ] ', ...
%         '-s 4x2x1vox  ', ...
%         '-f 3x2x1 -l 1 ', ...
%         '-m mattes[  %s, %s , 0.5 , 32 ] ', ...
%         '-m cc[  %s, %s , 0.5 , 4 ] ', ...
%         '-t SyN[ .20, 3, 0 ] ', ...
%         '-c [ 100x100x50,-0.01,5 ] ', ...
%         '-s 1x0.5x0vox ', ...
%         '-f 4x2x1 -l 1 -u 1 -z 1 ', ...
%         '-o [ %s,%s_diff.nii.gz,%s_inv.nii.gz]'],...
%         fixd, movn, fixd, movn, fixd, movn, fixd, movn, fixd, movn, fixd, movn, nomen, nomen, nomen ));
%         %
%         unix(sprintf('antsApplyTransforms -d 3 -i %s -r %s -n linear -t %s1Warp.nii.gz -t %s0GenericAffine.mat -o %s_warped.nii.gz',movn,fixd,nomen,nomen,nomen));
    else
        error('unrecognized warp option!')
    end

    disp('t1 mask n warp done! Now for some tissue segmentation...');
    
    % >> Pipeline Step #3: "ASEG"
    % >> parcellates the processed anatomical image from "AWARP" step into CSF, GM and WM regions; combined with "ROIREG" used to regress out physiological signals
    % >> takes in "anat_procss.nii" image and "anatBrainMask.nii" mask, produces anat_seg_CSF.nii, anat_seg_GM.nii, anat_seg_WM.nii 
    if strcmpi(PipeStruct_aug.ASEG{1},'FL1')
        aseg_FL1( sprintf('%s/anat_procss.nii.gz',opath2a), sprintf('%s/anatBrainMask.nii.gz',opath2a), opath3a );
    else
        error('unrecognized segmentation option')
    end
    
    % extra step: warping the anatomical segmentations into template space
    tisslist = {'CSF','GM','WM'}; % list tissues in increasing order of T1 intensity
    for i=1:3
        if ~exist(sprintf('%s/anat_seg_%s_warped.nii.gz',opath3a,tisslist{i}),'file')
            if contains(PipeStruct_aug.AWARP{1},'AF') %afni-styles warp
                unix(sprintf('3dNwarpApply -master %s/anat_warped.nii.gz -source %s/anat_seg_%s.nii.gz -nwarp "%s/anatQQ_WARP.nii.gz %s/anatQQ.aff12.1D" -prefix %s/anat_seg_%s_warped.nii.gz',opath2a,opath3a,tisslist{i},opath2a,opath2a,opath3a,tisslist{i}));
            else
                error('unrecognized warping style!')
            end
        else
            disp('skipping tissue seg warping...')
        end
    end

    fprintf('\n===> subj %s (%u/%u), physio. processing...\n',subject_list{ns},ns,numel(subject_list)),
    %% =======================================================================
    %%      PHYSIO Processing ...
    %% =======================================================================

    fprintf('\n===> phys-proc. now on subj %u/%u: %s...\n',ns,numel(subject_list),subject_list{ns}),
    
    disp('nothin for physio-proc so far. put in soon!')

    fprintf('\n===> subj %s (%u/%u), perf. processing (part-1)...\n',subject_list{ns},ns,numel(subject_list)),
    %% =======================================================================
    %%      PERF Processing, Block-1 ...
    %% =======================================================================

    % --> NB: block-1 processing gets split into prewarp/warp/postwarp subfolders for ease of debugging

    % motref from calibration file!!

    % ===[ Motref Basefile for subsequent spatial alignments etc ]===
    %
    motref_0rel_m0base = 0; %--- current default is to just take first of the run
    dlmwrite( sprintf('%s/motref_0rel_m0base.txt',opath1f), [motref_0rel_m0base] ); % ** found min-disp brick (zero-relative indexing)
    %
    % exrtact basefile
    if ~exist(sprintf('%s/prewarp/motref.nii.gz',opath2f),'file') || ~exist(sprintf('%s/warp/motref_masked_2std.nii.gz',opath2f),'file')
        % --extract the motref file (+maskit, for use in perf-estim)
        unix(sprintf('3dTcat -prefix %s/prewarp/motref.nii.gz ''%s/m0ref_cat.nii.gz[%u]''',opath2f,opath0,motref_0rel_m0base));
        % temporary, smooth it for a looser mask
        unix(sprintf('3dmerge -prefix %s/prewarp/motref_tmpsmo.nii.gz -doall -1blur_fwhm 6.0  %s/prewarp/motref.nii.gz',opath2f,opath2f));
        unix(sprintf('3dAutomask -prefix %s/prewarp/motref_smo_mask.nii.gz %s/prewarp/motref_tmpsmo.nii.gz',opath2f,opath2f)); % ** gets a (more) decent functional mask
        unix(sprintf('rm %s/prewarp/motref_tmpsmo.nii.gz',opath2f));
        % --deoblique it (+masked brain for spatial warping)
        unix(sprintf('3dWarp -oblique2card -prefix %s/warp/motref_deob.nii.gz -wsinc5 %s/prewarp/motref.nii.gz' , opath2f, opath2f));
        unix(sprintf('fslreorient2std %s/warp/motref_deob.nii.gz %s/warp/motref_2std.nii.gz',opath2f,opath2f));
        if exist(sprintf('%s/warp/motref_2std.nii.gz',opath2f),'file')
            unix(sprintf('rm %s/warp/motref_deob.nii.gz', opath2f)); % not really useful intermediate - delete it
        else
            error('failed to create deob/reoriented "motref" file')
        end
        % creat, then apply the mask
        unix(sprintf('3dAutomask -prefix %s/warp/motref_tmpmask.nii.gz %s/warp/motref_2std.nii.gz',opath2f,opath2f));
        unix(sprintf('3dcalc -prefix %s/warp/motref_masked_2std.nii.gz -a %s/warp/motref_2std.nii.gz -b %s/warp/motref_tmpmask.nii.gz -expr ''a*b''',opath2f,opath2f,opath2f));
        if exist(sprintf('%s/warp/motref_masked_2std.nii.gz',opath2f),'file')
            unix(sprintf('rm %s/warp/motref_tmpmask.nii.gz', opath2f)); % not really useful intermediate - delete it
        else
            error('failed to create masked reoriented "motref" file')
        end
    end

    for nr=1:InputStruct_ssa.N_perf

        % tag for rawdata in case scan dropping was performed
        if strcmpi( InputStruct_ssa.PWDROP, 'FIRST' ) || strcmpi( InputStruct_ssa.PWDROP, 'LAST' )
            drop_tag = '_drop';
        elseif strcmpi( InputStruct_ssa.PWDROP, 'NONE' )
            drop_tag = '';
        else
            error('unrecognized drop scheme! Should be NONE, FIRST or LAST!')
        end

        % INIMOT is constructing some preliminary estimates of displacement -- finds the minimum-displacement volume for motion correction later 
        % **** THIS IS ONLY FOR TCFILT (MAYBE OTHER STUFF?) AT THE MOMENT!
        if strcmpi(ParamStruct_aug.INIMOT,'OP1')
            motref_0rel = inimot_OP1( sprintf('%s/perf%u%s.nii.gz',opath0,nr,drop_tag), sprintf('perf%u',nr), sprintf('%s/init_mot_estim',opath1f) );
        else
            error('unrecognized initial motion estimator?!')
        end

        % >> Pipeline Step #4: "TCFILT"
        if strcmpi(PipeStruct_aug.TCFILT{1},'OP1')
            tcfilt_OP1( sprintf('%s/perf%u%s.nii.gz',opath0,nr,drop_tag), sprintf('perf%u',nr), sprintf('%s/prewarp',opath2f), PipeStruct_aug.TCFILT{2}, motref_0rel );
        elseif strcmpi(PipeStruct_aug.TCFILT{1},'OFF')
            unix(sprintf('cp %s/func%u%s.nii.gz %s/prewarp/func%u_tcfilt.nii.gz',opath0,nr,drop_tag, opath2f,nr));
        else
            error('unrecognized despiking?!')
        end
        % >> Pipeline Step #4: "PWALIGN"
        if strcmpi(PipeStruct_aug.PWALIGN{1},'AF1')
            pwalign_AF1( sprintf('%s/prewarp/perf%u_tcfilt.nii.gz',opath2f,nr), sprintf('%s/prewarp/motref.nii.gz',opath2f), sprintf('perf%u',nr), sprintf('%s/prewarp',opath2f) ); % currently just pulls first brick from each m0run as ref
            unix(sprintf('mv %s/prewarp/perf%u_mpe %s/warp',opath2f,nr,opath2f)); % push mpes into warp folder!
        elseif strcmpi(PipeStruct_aug.PWALIGN{1},'OFF')
            unix(sprintf('cp %s/perf%u_tcfilt.nii.gz %s/prewarp/perf%u_pwalign.nii.gz',opath2f,nr, opath2f,nr));
        else
            error('unrecognized m0-aligning?!')
        end
        % >> Pipeline Step #8: "PRESMO" -- NB done here, but not really used until P2
        if strcmpi(PipeStruct_aug.PRESMO{1},'AF1')
            presmo_AF1( sprintf('%s/prewarp/perf%u_pwalign.nii.gz',opath2f,nr), sprintf('perf%u',nr), sprintf('%s/prewarp',opath2f), PipeStruct_aug.PRESMO{2} );
        else
            error('unrecognized smoothing?!')
        end
    end

    % >> Pipeline Step #4: "PWALIGN"
    if strcmpi(PipeStruct_aug.PWALIGN{1},'AF1')
        pwalign_AF1( sprintf('%s/m0ref_cat.nii.gz',opath0), sprintf('%s/prewarp/motref.nii.gz',opath2f), sprintf('m0ref'), sprintf('%s/prewarp',opath2f) ); % currently just pulls first brick from each m0run as ref
        unix(sprintf('mv %s/prewarp/m0ref_mpe %s/warp',opath2f,opath2f)); % push mpes into warp folder!
    elseif strcmpi(PipeStruct_aug.PWALIGN{1},'OFF')
        unix(sprintf('cp %s/m0ref_cat.nii.gz %s/prewarp/m0ref_pwalign.nii.gz',opath0, opath2f));
    else
        error('unrecognized m0-aligning?!')
    end
    % >> Pipeline Step #8: "PRESMO" -- NB done here, but not really used until P2
    if strcmpi(PipeStruct_aug.PRESMO{1},'AF1')
        presmo_AF1( sprintf('%s/prewarp/m0ref_pwalign.nii.gz',opath2f), sprintf('m0ref'), sprintf('%s/prewarp',opath2f), PipeStruct_aug.PRESMO{2} );
    else
        error('unrecognized smoothing?!')
    end
    
    clear prefix_set Funcfile_set prefix_set_wrp Funcfile_set_wrp;
    
    for nr=1:InputStruct_ssa.N_perf

        %*** now run --> perfusion estimation(!)
        
        if strcmpi(ParamStruct_aug.PERF_MODEL,'ATBX1')
    
            % clear the struct
            dataInfo=[];
            % assign values to fields
            dataInfo.TR_MSEC  = InputStruct_ssa.TR_MSEC; % repetition time
            dataInfo.TE_MSEC  = InputStruct_ssa.TE_MSEC; % echo time
            dataInfo.TI1_MSEC = InputStruct_ssa.LDS_MSEC(1); % label time 
            dataInfo.TI2_MSEC = InputStruct_ssa.LDS_MSEC(2); % delay --> full inversion (label+postlabel) 
            dataInfo.TSL_MSEC = InputStruct_ssa.LDS_MSEC(3); % slice time (=0 if no-delay / -1 if auto-estim) 
            % plug in any prespec'd kinetic modelling values...
            par_list = fields(ParamStruct_aug);
            kmix = find( contains(par_list,'KM_'));
            for i=1:numel(kmix)
                dataInfo.(par_list{kmix(i)}) = ParamStruct_aug.(par_list{kmix(i)});
            end

            Input   = sprintf('%s/prewarp/perf%u_presmo.nii.gz',opath2f,nr);
            Calib   = sprintf('%s/prewarp/m0ref_presmo.nii.gz',opath2f);
            Mskname = sprintf('%s/prewarp/motref_smo_mask.nii.gz',opath2f);
            Oname   = sprintf('%s/prewarp/pwfit_%u',opath2f,nr); % perf data goes into prewarp for now, i guess?

            % run the perfusion estimation
            ATBX1( Input, Calib, Mskname, InputStruct_ssa.LAB_METHOD, InputStruct_ssa.LAB_ORDER, ParamStruct_aug.SUBTRACT_TYPE, dataInfo, Oname );                
        end

        % fixing orientation stuff -- adjusting for obliquity, switching to standard mni-compatible orientation 
        kflist = {sprintf('pwfit_%u_PWI',nr),sprintf('pwfit_%u_CBF',nr),sprintf('pwfit_%u_PWI_avg',nr),sprintf('pwfit_%u_CBF_avg',nr)};
        %
        for kf = 1:numel(kflist)
            if ~exist(sprintf('%s/prewarp/%s_2std.nii.gz',opath2f,kflist{kf}),'file')
                unix(sprintf('3dWarp -oblique2card -prefix %s/prewarp/%s_deob.nii.gz -wsinc5 %s/prewarp/%s.nii.gz' , opath2f,kflist{kf}, opath2f,kflist{kf}));
                unix(sprintf('fslreorient2std %s/prewarp/%s_deob.nii.gz %s/prewarp/%s_2std.nii.gz',opath2f,kflist{kf},opath2f,kflist{kf}));
                if exist(sprintf('%s/prewarp/%s_2std.nii.gz',opath2f,kflist{kf}),'file')
                    unix(sprintf('rm %s/prewarp/%s_deob.nii.gz', opath2f,kflist{kf})); % not really useful intermediate - delete it
                else
                    error('failed to create deob/reoriented %s files',kflist{kf})
                end
            else
                fprintf('skipping (%s) alignment prep...\n',kflist{kf})
            end
            % store fields for FWARP step later on...
            prefix_set{nr + InputStruct_ssa.N_perf*(kf-1)} = kflist{kf};
            Funcfile_set{nr + InputStruct_ssa.N_perf*(kf-1)} = sprintf('%s/prewarp/%s_2std.nii.gz',opath2f,kflist{kf});
            % for postwarp smoothing etc...
            prefix_set_wrp{nr + InputStruct_ssa.N_perf*(kf-1)} = sprintf('%s_warped',kflist{kf});
            Funcfile_set_wrp{nr + InputStruct_ssa.N_perf*(kf-1)} = sprintf('%s/postwarp/%s_warped.nii.gz',opath2f,kflist{kf});
        end
    end
    % ==> deobbing just the aligned m0set
    if ~exist(sprintf('%s/prewarp/m0ref_2std.nii.gz',opath2f),'file')
        unix(sprintf('3dWarp -oblique2card -prefix %s/prewarp/m0ref_deob.nii.gz -wsinc5 %s/prewarp/m0ref_pwalign.nii.gz' , opath2f, opath2f));
        unix(sprintf('fslreorient2std %s/prewarp/m0ref_deob.nii.gz %s/prewarp/m0ref_2std.nii.gz',opath2f,opath2f));
        if exist(sprintf('%s/prewarp/m0ref_2std.nii.gz',opath2f),'file')
            unix(sprintf('rm %s/prewarp/m0ref_deob.nii.gz', opath2f)); % not really useful intermediate - delete it
        else
            error('failed to create deob/reoriented m0ref files')
        end
    else
        fprintf('skipping (m0ref) alignment prep...\n')
    end
    % append to cell arrays --> only for warping step!
    prefix_set   = [prefix_set, 'm0ref'];
    Funcfile_set = [Funcfile_set, sprintf('%s/prewarp/m0ref_2std.nii.gz',opath2f)];


    % pre-specifying some input/output directories
    odir1 = sprintf('%s/warp',opath2f);
    odir2 = sprintf('%s/postwarp',opath2f);
    Anatloc = opath2a;

    % >> Pipeline Step #7: "PWFWARP"
    if strcmpi(PipeStruct_aug.PWWARP{1},'AF1')
        pwwarp_AF1( Funcfile_set, prefix_set, odir1, odir2, sprintf('%s/warp/motref_masked_2std.nii.gz',opath2f), Anatloc )
    elseif strcmpi(PipeStruct_aug.PWWARP{1},'AF2')
        pwwarp_AF2( Funcfile_set, prefix_set, odir1, odir2, sprintf('%s/warp/motref_masked_2std.nii.gz',opath2f), Anatloc )
    else
        error('unrecognized functional warping scheme')
    end

    for ni=1:numel(Funcfile_set_wrp) % exclude the m0ref
        % >> Pipeline Step #8: "POSTSMO" -- NB done here, but not really used until P2
        if strcmpi(PipeStruct_aug.POSTSMO{1},'AF1')
            postsmo_AF1( Funcfile_set_wrp{ni}, prefix_set_wrp{ni}, odir2, PipeStruct_aug.POSTSMO{2} );
        else
            error('unrecognized smoothing?!')
        end
    end

    % extra step: make mask, mean, sd maps --> run-1 mean and sd used for creating group-level maps. For runs >1, just kept for qc purposes
    if ~exist(sprintf('%s/postwarp/m0ref_warped_mask.nii.gz',opath2f),'file') 
        unix(sprintf('3dAutomask -prefix %s/postwarp/m0ref_warped_mask.nii.gz %s/postwarp/m0ref_warped.nii.gz',opath2f,opath2f));
    else
        disp('func mask done - skipping ahead.')
    end
    % extra step: tidied up functional mask using anatomical data - for later group mask construction
    if ~exist( sprintf('%s/postwarp/m0ref_warped_mask_clean.nii.gz',opath2f),'file')
        if ~exist(sprintf('%s/warp/anat_warped_rs.nii.gz',opath2f),'file') || ...
           ~exist(sprintf('%s/warp/anat_warped_rs_mask.nii.gz',opath2f),'file') || ...
           ~exist(sprintf('%s/postwarp/m0ref_warped_mask_clean.nii.gz',opath2f),'file')

            unix(sprintf('3dresample -master %s/postwarp/m0ref_warped_mask.nii.gz -input %s/anat_warped.nii.gz -prefix %s/warp/anat_warped_rs.nii.gz',opath2f,opath2a,opath2f));
            unix(sprintf('3dmask_tool -dilate_input 5 -5 -fill_holes -input %s/warp/anat_warped_rs.nii.gz -prefix %s/warp/anat_warped_rs_mask.nii.gz',opath2f,opath2f))
            unix(sprintf('3dmask_tool -input %s/postwarp/m0ref_warped_mask.nii.gz %s/warp/anat_warped_rs_mask.nii.gz -inter -prefix %s/postwarp/m0ref_warped_mask_clean.nii.gz',opath2f,opath2f,opath2f))
        else
            disp('clean m0ref mask done - skipping ahead.')
        end
    else
        disp('skipping newspace masking...')
    end
    % extra step: resampling tissue segmentations into functional space
    tisslist = {'CSF','GM','WM'}; % tissues in increasing order of T1 intensity
    for i=1:3
        if ~exist( sprintf('%s/anat_seg_%s_resam.nii.gz',opath3f,tisslist{i}),'file')
            unix(sprintf('3dresample -master %s/postwarp/m0ref_warped_mask_clean.nii.gz -input %s/anat_seg_%s_warped.nii.gz -prefix %s/anat_seg_%s_resam.nii.gz',...
                opath2f,opath3a,tisslist{i},opath3f,tisslist{i}));
        else
            disp('skipping tissue seg warping...')
        end
    end

end

%% INTERMEZZO-ANATOMICAL: getting group-level maps

if IMPORT_MASKFILES==0
    
    disp('now constructing Anatomical group level brain maps n masks');
    
    % constructing a group-level consensus mask, along with a probabilistic map of brain voxels
    if ~exist( [outpath,'/_group_level/masks/pipe_',PipeStruct_aug.PNAME{1},'/anat_brain_mask_grp.nii'], 'file') || ... % Brain mask
       ~exist( [outpath,'/_group_level/brain_maps/pipe_',PipeStruct_aug.PNAME{1},'/anat_pBRAIN_grp.nii'], 'file')
    
        for ni=1:numel(mask_subj_idxes)
            % check existence of subject specific struct file
            if exist(fullfile( outpath,subject_list{mask_subj_idxes(ni)},'InputStruct_ssa.mat'),'file')  
                load( fullfile( outpath,subject_list{mask_subj_idxes(ni)},'InputStruct_ssa.mat') )
            else
                error('cannot find Input struct file for subject: %s \n');
            end
            opath1a = fullfile( outpath, InputStruct_ssa.PREFIX,'anat_proc');
            opath2a = fullfile( opath1a,sprintf('subpipe_%03u',    PipeStruct_aug.pipe_idx.Warp));
            Vw = load_untouch_niiz(sprintf('%s/anatBrainMask_warped.nii.gz',opath2a)); %***%
            if ni==1
                volref = double(Vw.img);
            else
                volref = volref + double(Vw.img);
            end
        end
        volref = volref./numel(mask_subj_idxes);
        Vw.img = double(volref>0.50); % included in majority of individuals
        save_untouch_niiz(Vw,[outpath,'/_group_level/masks/pipe_',PipeStruct_aug.PNAME{1},'/anat_brain_mask_grp.nii']);
        Vw.img = volref; % prob
        Vw.hdr.dime.bitpix=32;
        Vw.hdr.dime.datatype=16;
        save_untouch_niiz(Vw,[outpath,'/_group_level/brain_maps/pipe_',PipeStruct_aug.PNAME{1},'/anat_pBRAIN_grp.nii']);
    end
    
    % constructing a group-level mean anatomical image
    if ~exist( [outpath,'/_group_level/brain_maps/pipe_',PipeStruct_aug.PNAME{1},'/anat_brain_grp.nii'], 'file') % Mean image
    
        for ni=1:numel(mask_subj_idxes)
        
            % check existence of subject specific struct file
            if exist(fullfile( outpath,subject_list{mask_subj_idxes(ni)},'InputStruct_ssa.mat'),'file')  
                load( fullfile( outpath,subject_list{mask_subj_idxes(ni)},'InputStruct_ssa.mat') )
            else
                error('cannot find Input struct file for subject: %s \n');
            end
            opath1a = fullfile( outpath, InputStruct_ssa.PREFIX,'anat_proc');
            opath2a = fullfile( opath1a,sprintf('subpipe_%03u',    PipeStruct_aug.pipe_idx.Warp));
            Vw = load_untouch_niiz(sprintf('%s/anat_warped.nii.gz',opath2a)); %***%
            if ni==1
                volref = double(Vw.img);
            else
                volref = volref + double(Vw.img);
            end
        end
        Vw.img = volref./numel(mask_subj_idxes);
        save_untouch_niiz(Vw,[outpath,'/_group_level/brain_maps/pipe_',PipeStruct_aug.PNAME{1},'/anat_brain_grp.nii']);
    end

    % constructing group-level probabilistic tissue maps of CSF, GM, WM
    tisslist = {'CSF','GM','WM'};
    for i=1:3
        if ~exist( [outpath,'/_group_level/brain_maps/pipe_',PipeStruct_aug.PNAME{1},'/anat_p',tisslist{i},'_grp.nii'], 'file') % Mean image
        
            for ni=1:numel(mask_subj_idxes)
            
                % check existence of subject specific struct file
                if exist(fullfile( outpath,subject_list{mask_subj_idxes(ni)},'InputStruct_ssa.mat'),'file')  
                    load( fullfile( outpath,subject_list{mask_subj_idxes(ni)},'InputStruct_ssa.mat') )
                else
                    error('cannot find Input struct file for subject: %s \n');
                end
                opath1a = fullfile( outpath, InputStruct_ssa.PREFIX,'anat_proc');
                opath2a = fullfile( opath1a,sprintf('subpipe_%03u',    PipeStruct_aug.pipe_idx.Warp));
                opath3a = fullfile( opath2a,sprintf('seg_subpipe_%03u',PipeStruct_aug.pipe_idx.Seg ));
                Vw = load_untouch_niiz(sprintf('%s/anat_seg_%s_warped.nii.gz',opath3a,tisslist{i})); %***%
                if ni==1
                    volref = double(Vw.img);
                else
                    volref = volref + double(Vw.img);
                end
            end
            Vw.img = volref./numel(mask_subj_idxes);
            save_untouch_niiz(Vw,[outpath,'/_group_level/brain_maps/pipe_',PipeStruct_aug.PNAME{1},'/anat_p',tisslist{i},'_grp.nii']);
        end
    end
      
    % rough "ventricular/sulcal map" for alignment checking
    if ~exist([outpath,'/_group_level/masks/pipe_',PipeStruct_aug.PNAME{1},'/anat_sulcal_mask_grp.nii'],'file')
        Va  = load_untouch_niiz([outpath,'/_group_level/brain_maps/pipe_',PipeStruct_aug.PNAME{1},'/anat_brain_grp.nii']);
        Vb  = load_untouch_niiz([outpath,'/_group_level/masks/pipe_',PipeStruct_aug.PNAME{1},'/anat_brain_mask_grp.nii']);
        vsm = double(Va.img);
        vsi = 1 - (vsm - min(vsm(:)))./(max(vsm(:))-min(vsm(:)));
        vsi = vsi .* double( Vb.img );
        V   = Va;
        V.img = double( vsi > 0.5);% prctile(vsi(vsi>0),75));
        save_untouch_niiz(V,[outpath,'/_group_level/masks/pipe_',PipeStruct_aug.PNAME{1},'/anat_sulcal_mask_grp.nii']);
        clear Va Vb V;
    end

else
    disp('using premade/imported Anatomical group level brain maps n masks');

    % pull represntative first file
    if exist(fullfile( outpath,subject_list{1},'InputStruct_ssa.mat'),'file')  
        load( fullfile( outpath,subject_list{1},'InputStruct_ssa.mat') )
    else
        error('cannot find Input struct file for subject: %s \n');
    end
    opath1a = fullfile( outpath, InputStruct_ssa.PREFIX,'anat_proc');
    opath2a = fullfile( opath1a,sprintf('subpipe_%03u',    PipeStruct_aug.pipe_idx.Warp));
    Vo = load_untouch_niiz(sprintf('%s/anat_warped.nii.gz',opath2a)); %***%
    
    % here is our checklist: subdirs and files
    chklist_a = {'masks',              'brain_maps',     'brain_maps',    'brain_maps',   'brain_maps',  'brain_maps',  'masks'}
    chklist_b = {'anat_brain_mask_grp','anat_pBRAIN_grp','anat_brain_grp','anat_pCSF_grp','anat_pGM_grp','anat_pWM_grp','anat_sulcal_mask_grp'}
      
    for k=1:numel(chklist_a)
        fileimp = [outpath,'/_group_level/',chklist_a{k},'/pipe_',PipeStruct_aug.PNAME{1},'/',chklist_b{k},'.nii'];
        % quick check to make sure everything is there
        % & quick check for compatibility with anaomical processed datas
        if exist( fileimp,'file' )
            Vi = load_untouch_niiz(fileimp);
        elseif exist( [fileimp,'.gz'],'file' )
            Vi = load_untouch_niiz([fileimp,'.gz']);
        else
            error('failed to find imported anatomical mask file %s!',fileimp);
        end

        if sum( Vo.hdr.dime.pixdim(2:5)==Vi.hdr.dime.pixdim(2:5) ) ~= 4
            error('imported masks do not match on vox-res!');
        end
        if sum( Vo.hdr.dime.dim(2:5)==Vi.hdr.dime.dim(2:5) ) ~= 4
            error('imported masks do not match on matrix-size!');
        end        

        Vo.hdr.dime.pixdim(2:5); %vox res
        Vo.hdr.dime.dim(2:5); %matx size
    end
end

%% INTERMEZZO-FUNCTIONAL: getting group-level maps

if IMPORT_MASKFILES==0

    disp('now constructing Functional group level brain maps n masks')
    
    % creating group-level consensus brain mask for restricting analysis
    if ~exist( [outpath,'/_group_level/masks/pipe_',PipeStruct_aug.PNAME{1},'/perf_brain_mask_grp.nii'], 'file') % Brain mask
        for ni=1:numel(mask_subj_idxes)
            % check existence of subject specific struct file
            if exist(fullfile( outpath,subject_list{mask_subj_idxes(ni)},'InputStruct_ssa.mat'),'file')
                load( fullfile( outpath,subject_list{mask_subj_idxes(ni)},'InputStruct_ssa.mat') )
            else
                error('cannot find Input struct file for subject: %s \n');
            end
            opath1f = fullfile(outpath,InputStruct_ssa.PREFIX,'perf_proc_p1');
            opath2f = fullfile(opath1f,sprintf('subpipe_%03u',PipeStruct_aug.pipe_idx.P1));
            Vw=load_untouch_niiz(sprintf('%s/postwarp/m0ref_warped_mask_clean.nii.gz',opath2f)); %***%
            figure,imagesc( Vw.img(:,:,30) );
            if ni==1
                volref = double(Vw.img);
            else
                volref = volref + double(Vw.img);
            end
        end
        volref = volref./numel(mask_subj_idxes);
        Vw.img = double(volref>0.50); % included in majority of individuals
        save_untouch_niiz(Vw,[outpath,'/_group_level/masks/pipe_',PipeStruct_aug.PNAME{1},'/perf_brain_mask_grp.nii']);
    else
        Vw = load_untouch_niiz([outpath,'/_group_level/masks/pipe_',PipeStruct_aug.PNAME{1},'/perf_brain_mask_grp.nii']);
    end
    mask_brain = double(Vw.img);
    
    % creating conservative group-level consensus CSF mask for ROIREG etc
    if ~exist( [outpath,'/_group_level/masks/pipe_',PipeStruct_aug.PNAME{1},'/perf_CSF_mask_grp.nii'], 'file') % CSF mask
        for ni=1:numel(mask_subj_idxes)
            % check existence of subject specific struct file
            if exist(fullfile( outpath,subject_list{mask_subj_idxes(ni)},'InputStruct_ssa.mat'),'file')
                load( fullfile( outpath,subject_list{mask_subj_idxes(ni)},'InputStruct_ssa.mat') )
            else
                error('cannot find Input struct file for subject: %s \n');
            end
            opath1f = fullfile(outpath,InputStruct_ssa.PREFIX,'perf_proc_p1');
            opath2f = fullfile(opath1f,sprintf('subpipe_%03u',PipeStruct_aug.pipe_idx.P1));
            opath3f = fullfile(opath2f,sprintf('seg_subpipe_%03u',PipeStruct_aug.pipe_idx.Seg ));
            Vw=load_untouch_niiz(sprintf('%s/anat_seg_CSF_resam.nii.gz',opath3f)); %***%
            if ni==1
                volref = double(Vw.img);
            else
                volref = volref + double(Vw.img);
            end
        end
        volref = volref./numel(mask_subj_idxes);
        Vw.img = clust_up( double( volref > 0.90 ) ,20);
        save_untouch_niiz(Vw,[outpath,'/_group_level/masks/pipe_',PipeStruct_aug.PNAME{1},'/perf_CSF_mask_grp.nii']);
    else
        Vw = load_untouch_niiz([outpath,'/_group_level/masks/pipe_',PipeStruct_aug.PNAME{1},'/perf_CSF_mask_grp.nii']);
    end
    mask_csf = double(Vw.img);
    
    % creating conservative group-level consensus WM mask for ROIREG etc
    if ~exist( [outpath,'/_group_level/masks/pipe_',PipeStruct_aug.PNAME{1},'/perf_WM_mask_grp.nii'], 'file') % WM mask
        for ni=1:numel(mask_subj_idxes)
            % check existence of subject specific struct file
            if exist(fullfile( outpath,subject_list{mask_subj_idxes(ni)},'InputStruct_ssa.mat'),'file')
                load( fullfile( outpath,subject_list{mask_subj_idxes(ni)},'InputStruct_ssa.mat') )
            else
                error('cannot find Input struct file for subject: %s \n');
            end
            opath1f = fullfile(outpath,InputStruct_ssa.PREFIX,'perf_proc_p1');
            opath2f = fullfile(opath1f,sprintf('subpipe_%03u',PipeStruct_aug.pipe_idx.P1));
            opath3f = fullfile(opath2f,sprintf('seg_subpipe_%03u',PipeStruct_aug.pipe_idx.Seg ));
            Vw=load_untouch_niiz(sprintf('%s/anat_seg_WM_resam.nii.gz',opath3f)); %***%
            if ni==1
                volref = double(Vw.img);
            else
                volref = volref + double(Vw.img);
            end
        end
        volref = volref./numel(mask_subj_idxes);
        volref = smooth3( volref, 'gaussian',[7 7 7], 0.85);
        Vw.img = clust_up( double( volref > 0.90 ) ,20);
        save_untouch_niiz(Vw,[outpath,'/_group_level/masks/pipe_',PipeStruct_aug.PNAME{1},'/perf_WM_mask_grp.nii']);
    else
        Vw = load_untouch_niiz([outpath,'/_group_level/masks/pipe_',PipeStruct_aug.PNAME{1},'/perf_WM_mask_grp.nii']);
    end
    mask_wm = double(Vw.img);
    
% %     % creating group-level consensus temporal variance mask for ROIREG etc
% %     if ~exist( [outpath,'/_group_level/masks/pipe_',PipeStruct_aug.PNAME{1},'/perf_tSD_mask_grp.nii'], 'file') % t-std mask
% %         for ni=1:numel(mask_subj_idxes)
% %             % check existence of subject specific struct file
% %             if exist(fullfile( outpath,subject_list{mask_subj_idxes(ni)},'InputStruct_ssa.mat'),'file')
% %                 load( fullfile( outpath,subject_list{mask_subj_idxes(ni)},'InputStruct_ssa.mat') )
% %             else
% %                 error('cannot find Input struct file for subject: %s \n');
% %             end
% %             opath1f = fullfile(outpath,InputStruct_ssa.PREFIX,'perf_proc_p1');
% %             opath2f = fullfile(opath1f,sprintf('subpipe_%03u',PipeStruct_aug.pipe_idx.P1));
% %             Vw=load_untouch_niiz(sprintf('%s/postwarp/perf1_warped_tsd.nii.gz',opath2f)); %***%
% %             if ni==1
% %                 volref = double(Vw.img);
% %             else
% %                 volref = volref + double(Vw.img);
% %             end
% %         end
% %         volref = volref./numel(mask_subj_idxes);
% %         Vw.img = clust_up( double( volref > prctile(volref(mask_brain>0),90) ) .* mask_brain ,20);
% %         save_untouch_niiz(Vw,[outpath,'/_group_level/masks/pipe_',PipeStruct_aug.PNAME{1},'/perf_tSD_mask_grp.nii']);
% %     else
% %         Vw = load_untouch_niiz([outpath,'/_group_level/masks/pipe_',PipeStruct_aug.PNAME{1},'/perf_tSD_mask_grp.nii']);
% %     end
% %     mask_tsd = double(Vw.img);
    
    % creating liberal group-level consensus GM mask for restricting analysis
    if ~exist( [outpath,'/_group_level/masks/pipe_',PipeStruct_aug.PNAME{1},'/perf_GM_mask_grp.nii'], 'file') % GM inclusive mask
        tisslist = {'CSF','GM','WM'};
        for i=1:3
            for ni=1:numel(mask_subj_idxes)
                % check existence of subject specific struct file
                if exist(fullfile( outpath,subject_list{mask_subj_idxes(ni)},'InputStruct_ssa.mat'),'file')
                    load( fullfile( outpath,subject_list{mask_subj_idxes(ni)},'InputStruct_ssa.mat') )
                else
                    error('cannot find Input struct file for subject: %s \n');
                end
                opath1f = fullfile(outpath,InputStruct_ssa.PREFIX,'perf_proc_p1');
                opath2f = fullfile(opath1f,sprintf('subpipe_%03u',PipeStruct_aug.pipe_idx.P1));
                opath3f = fullfile(opath2f,sprintf('seg_subpipe_%03u',PipeStruct_aug.pipe_idx.Seg ));
                Vw=load_untouch_niiz(sprintf('%s/anat_seg_%s_resam.nii.gz',opath3f,tisslist{i})); %***%
                if ni==1
                    volref = double(Vw.img);
                else
                    volref = volref + double(Vw.img);
                end
            end
            volref = volref./numel(mask_subj_idxes);
            volref = smooth3( volref, 'gaussian',[7 7 7], 0.85);
            volvec(:,i) = volref(mask_brain>0);
        end
        tmp = mask_brain; tmp(tmp>0)= double( (volvec(:,2) ./ sum(volvec,2)) >=0.33 );
        Vw.img = tmp;
        save_untouch_niiz(Vw,[outpath,'/_group_level/masks/pipe_',PipeStruct_aug.PNAME{1},'/perf_GM_mask_grp.nii']);
    else
        Vw = load_untouch_niiz([outpath,'/_group_level/masks/pipe_',PipeStruct_aug.PNAME{1},'/perf_GM_mask_grp.nii']);
    end
    mask_gm = double(Vw.img);
    
% %     % creating group-level mean average epi image for QC
% %     if ~exist( [outpath,'/_group_level/brain_maps/pipe_',PipeStruct_aug.PNAME{1},'/perf_tAV_grp.nii'], 'file') % group mean tav epi map
% %         for ni=1:numel(mask_subj_idxes)
% %             % check existence of subject specific struct file
% %             if exist(fullfile( outpath,subject_list{mask_subj_idxes(ni)},'InputStruct_ssa.mat'),'file')
% %                 load( fullfile( outpath,subject_list{mask_subj_idxes(ni)},'InputStruct_ssa.mat') )
% %             else
% %                 error('cannot find Input struct file for subject: %s \n');
% %             end
% %             opath1f = fullfile(outpath,InputStruct_ssa.PREFIX,'perf_proc_p1');
% %             opath2f = fullfile(opath1f,sprintf('subpipe_%03u',PipeStruct_aug.pipe_idx.P1));
% %             Vw=load_untouch_niiz(sprintf('%s/postwarp/perf1_warped_tav.nii.gz',opath2f)); %***%
% %             if ni==1
% %                 volref = double(Vw.img);
% %             else
% %                 volref = volref + double(Vw.img);
% %             end
% %         end
% %         Vw.img = volref./numel(mask_subj_idxes);
% %         save_untouch_niiz(Vw,[outpath,'/_group_level/brain_maps/pipe_',PipeStruct_aug.PNAME{1},'/perf_tAV_grp.nii']); % group mean sd epi map
% %     end
    
% %     % creating group-level mean temporal SD epi image for QC
% %     if ~exist( [outpath,'/_group_level/brain_maps/pipe_',PipeStruct_aug.PNAME{1},'/perf_tSD_grp.nii'], 'file')
% %         for ni=1:numel(mask_subj_idxes)
% %             % check existence of subject specific struct file
% %             if exist(fullfile( outpath,subject_list{mask_subj_idxes(ni)},'InputStruct_ssa.mat'),'file')
% %                 load( fullfile( outpath,subject_list{mask_subj_idxes(ni)},'InputStruct_ssa.mat') )
% %             else
% %                 error('cannot find Input struct file for subject: %s \n');
% %             end
% %             opath1f = fullfile(outpath,InputStruct_ssa.PREFIX,'perf_proc_p1');
% %             opath2f = fullfile(opath1f,sprintf('subpipe_%03u',PipeStruct_aug.pipe_idx.P1));
% %             Vw=load_untouch_niiz(sprintf('%s/postwarp/perf1_warped_tsd.nii.gz',opath2f)); %***%
% %             if ni==1
% %                 volref = double(Vw.img);
% %             else
% %                 volref = volref + double(Vw.img);
% %             end
% %         end
% %         Vw.img = volref./numel(mask_subj_idxes);
% %         save_untouch_niiz(Vw,[outpath,'/_group_level/brain_maps/pipe_',PipeStruct_aug.PNAME{1},'/perf_tSD_grp.nii']);
% %     end

else

    disp('using premade/imported Functional group level brain maps n masks');

    % pull represntative first file
    if exist(fullfile( outpath,subject_list{1},'InputStruct_ssa.mat'),'file')  
        load( fullfile( outpath,subject_list{1},'InputStruct_ssa.mat') )
    else
        error('cannot find Input struct file for subject: %s \n');
    end
    opath1f = fullfile(outpath,InputStruct_ssa.PREFIX,'perf_proc_p1');
    opath2f = fullfile(opath1f,sprintf('subpipe_%03u',PipeStruct_aug.pipe_idx.P1));
    Vo=load_untouch_niiz(sprintf('%s/postwarp/pwfit_1_PWI_avg_warped.nii.gz',opath2f)); %***%
    
    % here is our checklist: subdirs and files
%     chklist_a = {'masks','masks','masks','masks','masks','brain_maps','brain_maps'};
%     chklist_b = {'perf_brain_mask_grp','perf_CSF_mask_grp','perf_WM_mask_grp','perf_tSD_mask_grp','perf_GM_mask_grp','perf_tAV_grp','perf_tSD_grp'};
    chklist_a = {'masks','masks','masks','masks'};
    chklist_b = {'perf_brain_mask_grp','perf_CSF_mask_grp','perf_WM_mask_grp','perf_GM_mask_grp'};
      
    for k=1:numel(chklist_a)
        fileimp = [outpath,'/_group_level/',chklist_a{k},'/pipe_',PipeStruct_aug.PNAME{1},'/',chklist_b{k},'.nii'],
        % quick check to make sure everything is there
        % & quick check for compatibility with anaomical processed datas
        if exist( fileimp,'file' )
            Vi = load_untouch_niiz(fileimp);
        elseif exist( [fileimp,'.gz'],'file' )
            Vi = load_untouch_niiz([fileimp,'.gz']);
        else
            error('failed to find imported functional mask file %s!',fileimp);
        end

        if sum( Vo.hdr.dime.pixdim(2:5)==Vi.hdr.dime.pixdim(2:5) ) ~= 4
            error('imported masks do not match on vox-res!');
        end
        if sum( Vo.hdr.dime.dim(2:5)==Vi.hdr.dime.dim(2:5) ) ~= 4
            error('imported masks do not match on matrix-size!');
        end        

        Vo.hdr.dime.pixdim(2:5); %vox res
        Vo.hdr.dime.dim(2:5); %matx size
    end

end

for ns=subj_list_for_proc % step through anat-proc, func-proc (block-1)

    % check existence of subject specific struct file
    if exist(fullfile( outpath,subject_list{ns},'InputStruct_ssa.mat'),'file')  
        load( fullfile( outpath,subject_list{ns},'InputStruct_ssa.mat') )
    else
        error('cannot find Input struct file for subject: %s \n');
    end

    % quick formatting stuff, again assuRImes that directory structure was already constructed in "P0" pipeline step 
    opath0 = fullfile(outpath,InputStruct_ssa.PREFIX,'rawdata');
    %
    opath1a = fullfile( outpath, InputStruct_ssa.PREFIX,'anat_proc');
    opath2a = fullfile( opath1a,sprintf('subpipe_%03u',    PipeStruct_aug.pipe_idx.Warp));
    opath3a = fullfile( opath2a,sprintf('seg_subpipe_%03u',PipeStruct_aug.pipe_idx.Seg ));
    %
    opath1p = fullfile(outpath,InputStruct_ssa.PREFIX,'phys_proc');
    %
    opath1f = fullfile(outpath,InputStruct_ssa.PREFIX,'perf_proc_p1');
    opath2f = fullfile(opath1f,sprintf('subpipe_%03u',PipeStruct_aug.pipe_idx.P1));
    opath3f = fullfile(opath2f,sprintf('seg_subpipe_%03u',PipeStruct_aug.pipe_idx.Seg ));
    %
    opath4f = fullfile(outpath,InputStruct_ssa.PREFIX,'perf_proc_p2',['pipe_',PipeStruct_aug.PNAME{1}]); 

    fprintf('\n===> subj %s (%u/%u), perf. processing (part-2)...\n',subject_list{ns},ns,numel(subject_list)),
    %% =======================================================================
    %%      PERF Processing, Block-2 ...
    %% =======================================================================

    % just a lil bit of house keeping at the moment...

    %--> load the average, smoothed CBF/PWI/m0ref
    %--> copy em over
    %--> masked too
    for nr=1:InputStruct_ssa.N_perf

        if ~exist( [opath4f,'/perf',num2str(nr),'_fullproc_mpc.mat'],'file') %only do runs with missing outputs
            
            MBstr = ([outpath,'/_group_level/masks/pipe_',PipeStruct_aug.PNAME{1},'/perf_brain_mask_grp.nii']);
            MB = load_untouch_niiz(MBstr);   
            %
            filetags = {'m0ref',sprintf('pwfit_%u_CBF_avg',nr),sprintf('pwfit_%u_PWI_avg',nr)};
            smotags = {'','_postsmo','_postsmo'};
            for fi=1:numel(filetags)
                unix(sprintf('cp %s/postwarp/%s_warped%s.nii.gz %s/%s_fullproc.nii.gz',opath2f,filetags{fi},smotags{fi},opath4f,filetags{fi}))
                VSstr = sprintf('%s/%s_fullproc.nii.gz',opath4f,filetags{fi}); %***%
                VS = load_untouch_niiz(VSstr);
                volmat_mpc(:,fi) = mean( nifti_to_mat(VS,MB), 2);
            end

            save([opath4f,'/perf',num2str(nr),'_fullproc_mpc.mat'],'volmat_mpc');
        end
    end
         
end

disp('funxionale block-2 done');

