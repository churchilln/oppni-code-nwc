function P2_fmri_dataProcessing( inputfile, pipelinefile, paramlist, outpath, big_skip, mask_subj_idxes, subj_subset )
%
% =========================================================================
% P2_FMRI_DATAPROCESSING: this script should be run after the "P0" and "P1"
% steps of the bold pipeline. It does the actual data processing!
% =========================================================================
%
% Syntax:
%
%     P2_dataProcessing( inputfile, pipelinefile, paramlist, outpath, big_skip, mask_subj_idxes )
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
CODE_PATH = fileparts(which('P0_fmri_populateDirectories.m'));
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
[subject_list, InputStruct_aug, PipeStruct_aug, ParamStruct_aug] = P0_fmri_populateDirectories( inputfile, pipelinefile, paramlist, outpath );

% now augmenting outpath... do this after P0!
outpath = fullfile(outpath,'fmri_proc'); % subdir should be fmri_proc
if ~exist(outpath,'dir') error('fmri proc directory does not exist!'); end
e=dir([outpath,'*']); % dir to get absolute
outpath = fullfile(e.folder,e.name); % convert to absolute
% check for missing files from previous step too
File_Existence_Checker_fmri(InputStruct_aug,outpath,1); 

% check validity of analysis model etc. --> carry forward information about the analysis too!
analysis_struct = check_fmri_analysis_model( ParamStruct_aug.ANALYSIS );
% also, modify "number of components" field if more than one contrast is specified
if( ~isempty(strfind(ParamStruct_aug.CONTRAST,',' )) )
    analysis_struct.num_comp = 'multi_component';
end
% % % check validity of pipelines --> WILL SOON DO MORE??
pipeline_struct = check_fmri_processing_model( PipeStruct_aug );

% Now configuring param defaults if unspecified
%
if ~isfield(ParamStruct_aug,'INIMOT')
    ParamStruct_aug.INIMOT={'OP1',[]};
end
if ~isfield(ParamStruct_aug,'ROIMASK')
    ParamStruct_aug.ROIMASK={'OP1',[]};
end
if ~isfield(ParamStruct_aug,'TRUNC_ANL')
    ParamStruct_aug.TRUNC_ANL=[];
end
if ~isfield(ParamStruct_aug,'GMMASK_ANL')
    ParamStruct_aug.GMMASK_ANL=[];
end


% list of subjects for constructing group masks
IMPORT_MASKFILES=0;
if isempty(mask_subj_idxes)
    disp('using all subj in current input file for mask construction!')
    mask_subj_idxes = 1:numel(subject_list);
elseif ischar(mask_subj_idxes) 
    
    if ~exist(mask_subj_idxes,'file')
        error('mask id file not found!')
    end
    
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
    opath1f = fullfile(outpath,InputStruct_ssa.PREFIX,'func_proc_p1');
    opath2f = fullfile(opath1f,sprintf('subpipe_%03u',PipeStruct_aug.pipe_idx.P1));
    opath3f = fullfile(opath2f,sprintf('seg_subpipe_%03u',PipeStruct_aug.pipe_idx.Seg ));
    %
    opath4f = fullfile(outpath,InputStruct_ssa.PREFIX,'func_proc_p2',['pipe_',PipeStruct_aug.PNAME{1}]); 

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

    % >>> Spatial Masking
    Step = 'AMASK';
    if strcmpi(PipeStruct_aug.(Step){1},'OFF')
        error('cannot turn this step off!');
    else
        % get function handle for analysis model of interest
        currPath=pwd;                               % get current path
        cd(pipeline_struct.(Step).filepath);               % jump to module directory
        pfun= str2func(pipeline_struct.(Step).model_name); % get function handle
        cd(currPath);                               % jump back to current path
        % execute step:
        pfun( sprintf('%s/anat%u_2std.nii.gz', opath1a,nr), ParamStruct_aug.TEMPLATE_loc, opath2a, PipeStruct_aug.(Step)(2:end) );  
    end

    % >>> Spatial Warping
    Step = 'AWARP';
    if strcmpi(PipeStruct_aug.(Step){1},'OFF')
        error('cannot turn this step off!');
    else
        % get function handle for analysis model of interest
        currPath=pwd;                               % get current path
        cd(pipeline_struct.(Step).filepath);               % jump to module directory
        pfun= str2func(pipeline_struct.(Step).model_name); % get function handle
        cd(currPath);                               % jump back to current path
        % execute step:
        pfun( sprintf('%s/anat%u_2std.nii.gz', opath1a,nr), sprintf('%s/anatBrainMask.nii.gz',opath2a), ParamStruct_aug.TEMPLATE_loc, opath2a, PipeStruct_aug.(Step)(2:end) );  
    end

    % >>> Anatomic segmentation
    Step = 'ASEG';
    if strcmpi(PipeStruct_aug.(Step){1},'OFF')
        error('cannot turn this step off!');
    else
        % get function handle for analysis model of interest
        currPath=pwd;                               % get current path
        cd(pipeline_struct.(Step).filepath);               % jump to module directory
        pfun= str2func(pipeline_struct.(Step).model_name); % get function handle
        cd(currPath);                               % jump back to current path
        % execute step:
        pfun( sprintf('%s/anat_procss.nii.gz',opath2a), sprintf('%s/anatBrainMask.nii.gz',opath2a), opath3a, PipeStruct_aug.(Step)(2:end) );  
    end
    
    % extra step: warping the anatomical segmentations into template space
    tisslist = {'CSF','GM','WM'}; % list tissues in increasing order of T1 intensity
    for i=1:3
        if ~exist(sprintf('%s/anat_seg_%s_warped.nii.gz',opath3a,tisslist{i}),'file')
            if contains(PipeStruct_aug.AWARP{1},'AF') %afni-styles warp
                unix(sprintf('3dNwarpApply -master %s/anat_warped.nii.gz -source %s/anat_seg_%s.nii.gz -nwarp "%s/anatQQ_WARP.nii.gz %s/anatQQ.aff12.1D" -prefix %s/anat_seg_%s_warped.nii.gz',...
                    opath2a,opath3a,tisslist{i},opath2a,opath2a,opath3a,tisslist{i}));
            elseif contains(PipeStruct_aug.AWARP{1},'AN') %ants-styles warp
                unix(sprintf('antsApplyTransforms -d 3 -i %s/anat_seg_%s.nii.gz -r %s/anat_warped.nii.gz -n linear -t %s/anatQQ_WARP.nii.gz -t %s/anatQQ_GenericAffine.mat -o %s/anat_seg_%s_warped.nii.gz',...
                    opath3a,tisslist{i}, opath2a, opath2a,opath2a, opath3a,tisslist{i}));
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

    fprintf('\n===> subj %s (%u/%u), func. processing (part-1)...\n',subject_list{ns},ns,numel(subject_list)),
    %% =======================================================================
    %%      FUNC Processing, Block-1 ...
    %% =======================================================================

    %%%%%========== compatibility adjustments for BLOCK1 ... RICOR only
    missing_physio=0;
    for nr=1:InputStruct_ssa.N_func
        if isempty(InputStruct_ssa.frun(nr).PHYSIO_filename)
            missing_physio=missing_physio+1;
        end
    end
    if missing_physio>0 && ~strcmpi(PipeStruct_aug.RICOR{1},'OFF')
        warning('subject %s has %u/%u runs with missing physio. Turning RICOR off!',InputStruct_ssa.PREFIX, missing_physio,InputStruct_ssa.N_func)
        PipeStruct_aug.RICOR{1}='OFF';
    end
    %%%%%========== compatibility adjustments, done

    % --> NB: block-1 processing gets split into prewarp/warp/postwarp subfolders for ease of debugging
    
    % clear for variable run lengths between subj.
    clear prefix_set Funcfile_set base_set prefix_set_wrp Funcfile_set_wrp;

    for nr=1:InputStruct_ssa.N_func

        % tag for rawdata in case scan dropping was performed
        if InputStruct_ssa.frun(nr).DROP_first>0 || InputStruct_ssa.frun(nr).DROP_last>0
            drop_tag = '_drop';
        else
            drop_tag = '';
        end

        % * the steps below (INIMOT, DESPIKE, RICOR, TSHIFT) do slice-specific processing before 
        % we do any deobliqueing/alignment which destroys slice-based information
        
        % INIMOT is constructing some preliminary estimates of displacement -- finds the minimum-displacement volume for motion correction later 
        if strcmpi(ParamStruct_aug.INIMOT{1},'OP1')
            motref_0rel = inimot_OP1( sprintf('%s/func%u%s.nii.gz',opath0,nr,drop_tag), sprintf('func%u',nr), sprintf('%s/init_mot_estim',opath1f) );
        else
            error('unrecognized initial motion estimator?!')
        end

        % >>> Removing "Spikes" in fMRI data
        Step = 'DESPIKE';
        if strcmpi(PipeStruct_aug.(Step){1},'OFF')
            unix(sprintf('cp %s/func%u%s.nii.gz %s/prewarp/func%u_despike.nii.gz',opath0,nr,drop_tag, opath2f,nr));
        else
            % get function handle for analysis model of interest
            currPath=pwd;                               % get current path
            cd(pipeline_struct.(Step).filepath);               % jump to module directory
            pfun= str2func(pipeline_struct.(Step).model_name); % get function handle
            cd(currPath);                               % jump back to current path
            % execute step:
            pfun( sprintf('%s/func%u%s.nii.gz',opath0,nr,drop_tag), sprintf('func%u',nr), sprintf('%s/prewarp',opath2f), motref_0rel, PipeStruct_aug.(Step)(2:end) );  
        end

        % >> Pipeline Step #5: "RICOR"
        if strcmpi(PipeStruct_aug.RICOR{1},'AF1')
                    %-- c. slice-based physio correction
                    %
                    % get the stripped-down physio data matrix
                    fid = fopen( sprintf('%s.slibase.1D',ostr4) );
                    tline = fgetl(fid);
                    kq=0; ise=0;
                    while ischar(tline) 
                        if    ( contains(tline,'# >') ) ise=1;
                        elseif( contains(tline,'# <') ) ise=2;
                        elseif( ise==1 ) % only if currently flagged on
                            kq=kq+1;
                            ricormat(kq,:) = str2num( tline );
                        end
                        tline = fgetl(fid);
                    end
                    fclose(fid);
                    % take the matrix of slicewise physio covariates and regress slice-by-slice from the despiked data
                    ricor_regress( sprintf('%s/prewarp/func%u_despike.nii.gz',opath2f,nr), ricormat, sprintf('%s/prewarp/func%u_ricor.nii.gz',opath2f,nr) );
        elseif strcmpi(PipeStruct_aug.RICOR{1},'OFF')
            unix(sprintf('cp %s/prewarp/func%u_despike.nii.gz %s/prewarp/func%u_ricor.nii.gz',opath2f,nr, opath2f,nr));
        else
            error('unrecognized ricorring?!')
        end

        % >>> Slice-Timing Correction
        Step = 'TSHIFT';
        if strcmpi(PipeStruct_aug.(Step){1},'OFF')
            unix(sprintf('cp %s/prewarp/func%u_ricor.nii.gz %s/prewarp/func%u_tshift.nii.gz',opath2f,nr, opath2f,nr));
        else
            % get function handle for analysis model of interest
            currPath=pwd;                               % get current path
            cd(pipeline_struct.(Step).filepath);               % jump to module directory
            pfun= str2func(pipeline_struct.(Step).model_name); % get function handle
            cd(currPath);                               % jump back to current path
            % execute step:
            pfun( sprintf('%s/prewarp/func%u_ricor.nii.gz',opath2f,nr), sprintf('func%u',nr), sprintf('%s/prewarp',opath2f), InputStruct_ssa.TPATTERN, PipeStruct_aug.(Step)(2:end) );  
        end

        % fixing orientation stuff -- adjusting for obliquity, switching to standard mni-compatible orientation 
        if ~exist(sprintf('%s/prewarp/func%u_2std.nii.gz',opath2f,nr),'file')
            unix(sprintf('3dWarp -oblique2card -prefix %s/prewarp/func%u_deob.nii.gz -wsinc5 %s/prewarp/func%u_tshift.nii.gz' , opath2f,nr, opath2f,nr));
            unix(sprintf('fslreorient2std %s/prewarp/func%u_deob.nii.gz %s/prewarp/func%u_2std.nii.gz',opath2f,nr,opath2f,nr));
            if exist(sprintf('%s/prewarp/func%u_2std.nii.gz',opath2f,nr),'file')
                unix(sprintf('rm %s/prewarp/func%u_deob.nii.gz', opath2f,nr)); % not really useful intermediate - delete it
            else
                error('failed to create deob/reoriented functional file')
            end
        else
            disp('skipping alignment prep...')
        end

        % store fields for FWARP step later on
        prefix_set{nr} = sprintf('func%u',nr);
        Funcfile_set{nr} = sprintf('%s/prewarp/func%u_2std.nii.gz',opath2f,nr);
        base_set(nr) = motref_0rel;
        % and for SMOTTINGT
        prefix_set_wrp{nr} = sprintf('func%u_warped',nr);
        Funcfile_set_wrp{nr} = sprintf('%s/postwarp/func%u_warped.nii.gz',opath2f,nr);
    end

    % pre-specifying some input/output directories
    odir1 = sprintf('%s/warp',opath2f);
    odir2 = sprintf('%s/postwarp',opath2f);
    Anatloc = opath2a;

    % >>> Functional Warping
    Step = 'FWARP';
    if strcmpi(PipeStruct_aug.(Step){1},'OFF')
        error('cannot turn this step off!');
    else
        % get function handle for analysis model of interest
        currPath=pwd;                               % get current path
        cd(pipeline_struct.(Step).filepath);               % jump to module directory
        pfun= str2func(pipeline_struct.(Step).model_name); % get function handle
        cd(currPath);                               % jump back to current path
        % execute step:
        pfun( Funcfile_set, prefix_set, odir1, odir2, base_set, Anatloc, PipeStruct_aug.(Step)(2:end) );  
    end

    % >>> Spatial Smoothing
    Step = 'SMOOTH';
    if strcmpi(PipeStruct_aug.(Step){1},'OFF')
        error('cannot turn this step off!');
    else
        % get function handle for analysis model of interest
        currPath=pwd;                               % get current path
        cd(pipeline_struct.(Step).filepath);               % jump to module directory
        pfun= str2func(pipeline_struct.(Step).model_name); % get function handle
        cd(currPath);                               % jump back to current path
        % execute step:
        for ni=1:numel(Funcfile_set_wrp)
            pfun( Funcfile_set_wrp{ni}, prefix_set_wrp{ni}, odir2, PipeStruct_aug.(Step)(2:end) );  
        end
    end

    % extra step: make mask, mean, sd maps --> run-1 mean and sd used for creating group-level maps. For runs >1, just kept for qc purposes
    for nr=1:InputStruct_ssa.N_func
        if ~exist(sprintf('%s/postwarp/func%u_warped_mask.nii.gz',opath2f,nr),'file') || ...
           ~exist(sprintf('%s/postwarp/func%u_warped_tav.nii.gz',opath2f,nr),'file')  || ...
           ~exist(sprintf('%s/postwarp/func%u_warped_tsd.nii.gz',opath2f,nr),'file') 

            unix(sprintf('3dAutomask -prefix %s/postwarp/func%u_warped_mask.nii.gz %s/postwarp/func%u_warped.nii.gz',opath2f,nr,opath2f,nr));
            unix(sprintf('3dTstat -mean  -prefix %s/postwarp/func%u_warped_tav.nii.gz %s/postwarp/func%u_warped.nii.gz',opath2f,nr,opath2f,nr));
            unix(sprintf('3dTstat -stdev -prefix %s/postwarp/func%u_warped_tsd.nii.gz %s/postwarp/func%u_warped.nii.gz',opath2f,nr,opath2f,nr));
        else
            disp('func mask done - skipping ahead.')
        end
    end
    % extra step: tidied up functional mask using anatomical data - for later group mask construction
    if ~exist( sprintf('%s/postwarp/func%u_warped_mask_clean.nii.gz',opath2f,nr),'file')
        nr=1; % run-1 only
        if ~exist(sprintf('%s/warp/anat_warped_rs.nii.gz',opath2f),'file') || ...
           ~exist(sprintf('%s/warp/anat_warped_rs_mask.nii.gz',opath2f),'file') || ...
           ~exist(sprintf('%s/postwarp/func%u_warped_mask_clean.nii.gz',opath2f,nr),'file')

            unix(sprintf('3dresample -master %s/postwarp/func%u_warped_mask.nii.gz -input %s/anat_warped.nii.gz -prefix %s/warp/anat_warped_rs.nii.gz',opath2f,nr,opath2a,opath2f));
            unix(sprintf('3dmask_tool -dilate_input 5 -5 -fill_holes -input %s/warp/anat_warped_rs.nii.gz -prefix %s/warp/anat_warped_rs_mask.nii.gz',opath2f,opath2f))
            unix(sprintf('3dmask_tool -input %s/postwarp/func%u_warped_mask.nii.gz %s/warp/anat_warped_rs_mask.nii.gz -inter -prefix %s/postwarp/func%u_warped_mask_clean.nii.gz',opath2f,nr,opath2f,opath2f,nr))
        else
            disp('clean func mask done - skipping ahead.')
        end
    else
        disp('skipping newspace masking...')
    end
    % extra step: resampling tissue segmentations into functional space
    tisslist = {'CSF','GM','WM'}; % tissues in increasing order of T1 intensity
    for i=1:3
        if ~exist( sprintf('%s/anat_seg_%s_resam.nii.gz',opath3f,tisslist{i}),'file')
            unix(sprintf('3dresample -master %s/postwarp/func1_warped_mask_clean.nii.gz -input %s/anat_seg_%s_warped.nii.gz -prefix %s/anat_seg_%s_resam.nii.gz',...
                opath2f,opath3a,tisslist{i},opath3f,tisslist{i}));
        else
            disp('skipping tissue seg warping...')
        end
    end

% %     % SCRATCH -- extra step: effex maps
% %     mkdir(sprintf('%s/effex',opath2f));
% %     if ns==20
% %         return;
% %     end
% %     %
% %     V1 = load_untouch_niiz(sprintf('%s/func%u%s.nii',opath0,nr,drop_tag));
% %     V2 = load_untouch_niiz(sprintf('%s/prewarp/func%u_despike.nii',opath2f,nr));
% %     map = kurtosis(double(V2.img),0,4)-kurtosis(double(V1.img),0,4);
% %     map(~isfinite(map))=eps;
% %     mosaic_viewer( map, 3, [], [], 'jet', 1 )
end

%% Checkpoint

% check for all completed base-proccing of participants before masking can
% start (in case only a subset were run):
nrmax = 0;
for ns=1:numel(subject_list)
    % check existence of subject specific struct file
    if exist(fullfile( outpath,subject_list{ns},'InputStruct_ssa.mat'),'file')  
        load( fullfile( outpath,subject_list{ns},'InputStruct_ssa.mat') )
    else
        error('cannot find Input struct file for subject: %s \n');
    end
    nrmax = max([nrmax InputStruct_ssa.N_func]);
end
complet_mat = NaN*ones(numel(subject_list),nrmax);

for ns=1:numel(subject_list)
    
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
    opath1f = fullfile(outpath,InputStruct_ssa.PREFIX,'func_proc_p1');
    opath2f = fullfile(opath1f,sprintf('subpipe_%03u',PipeStruct_aug.pipe_idx.P1));
    opath3f = fullfile(opath2f,sprintf('seg_subpipe_%03u',PipeStruct_aug.pipe_idx.Seg ));
    %
    opath4f = fullfile(outpath,InputStruct_ssa.PREFIX,'func_proc_p2',['pipe_',PipeStruct_aug.PNAME{1}]); 

    for nr=1:InputStruct_ssa.N_func
        complet_mat(ns,nr)=2; %--fully processed already
        if ~exist( [opath4f,'/func',num2str(nr),'_fullproc.mat'],'file') %only do runs with missing outputs
            complet_mat(ns,nr)=1; %--processed up to masking
            if ~exist(sprintf('%s/postwarp/func%u_warped_smo.nii.gz',opath2f,nr),'file')
                complet_mat(ns,nr)=0; %--not processed
            end
        end
    end
end
if sum(complet_mat(:)==0)>0
    disp('not all subjects processed enough to mask. Halting for now!');
    disp('Unfinished:')
    ix = find(sum(complet_mat==0,2)>0);
    for i=1:numel(ix)
        ix2 = find( complet_mat(ix(i),:)==0);
        for j=1:numel(ix2)
            flagd = complet_mat(ix(i),ix2(j));
            fprintf('%s -- run %u.\n',subject_list{ix(i)},ix2(j));
        end
    end
    return;
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
        Vw2=Vw;
        Vw2.img = volref; % prob
        Vw2.hdr.dime.bitpix=32;
        Vw2.hdr.dime.datatype=16;
        save_untouch_niiz(Vw2,[outpath,'/_group_level/brain_maps/pipe_',PipeStruct_aug.PNAME{1},'/anat_pBRAIN_grp.nii']);
        clear Vw Vw2;
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
        volref = volref./numel(mask_subj_idxes);
        Vw.img = volref;
        save_untouch_niiz(Vw,[outpath,'/_group_level/brain_maps/pipe_',PipeStruct_aug.PNAME{1},'/anat_brain_grp.nii']);
        clear Vw;
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
            volref = volref./numel(mask_subj_idxes);
            Vw.img = volref;
            save_untouch_niiz(Vw,[outpath,'/_group_level/brain_maps/pipe_',PipeStruct_aug.PNAME{1},'/anat_p',tisslist{i},'_grp.nii']);
            clear Vw;
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
    chklist_a = {'masks',              'brain_maps',     'brain_maps',    'brain_maps',   'brain_maps',  'brain_maps',  'masks'};
    chklist_b = {'anat_brain_mask_grp','anat_pBRAIN_grp','anat_brain_grp','anat_pCSF_grp','anat_pGM_grp','anat_pWM_grp','anat_sulcal_mask_grp'};
      
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

maskisnew=0;

if IMPORT_MASKFILES==0

    disp('now constructing Functional group level brain maps n masks')
    
    % creating group-level consensus brain mask for restricting analysis
    if ~exist( [outpath,'/_group_level/masks/pipe_',PipeStruct_aug.PNAME{1},'/func_brain_mask_grp.nii'], 'file') % Brain mask
        for ni=1:numel(mask_subj_idxes)
            % check existence of subject specific struct file
            if exist(fullfile( outpath,subject_list{mask_subj_idxes(ni)},'InputStruct_ssa.mat'),'file')
                load( fullfile( outpath,subject_list{mask_subj_idxes(ni)},'InputStruct_ssa.mat') )
            else
                error('cannot find Input struct file for subject: %s \n');
            end
            opath1f = fullfile(outpath,InputStruct_ssa.PREFIX,'func_proc_p1');
            opath2f = fullfile(opath1f,sprintf('subpipe_%03u',PipeStruct_aug.pipe_idx.P1));
            Vw=load_untouch_niiz(sprintf('%s/postwarp/func1_warped_mask_clean.nii.gz',opath2f)); %***%
            if ni==1
                volref = double(Vw.img);
            else
                volref = volref + double(Vw.img);
            end
        end
        volref = volref./numel(mask_subj_idxes);
        Vw.img = double(volref>0.50); % included in majority of individuals
        save_untouch_niiz(Vw,[outpath,'/_group_level/masks/pipe_',PipeStruct_aug.PNAME{1},'/func_brain_mask_grp.nii']);
        Vw2=Vw;
        Vw2.img = volref; % prob
        Vw2.hdr.dime.bitpix=32;
        Vw2.hdr.dime.datatype=16;
        save_untouch_niiz(Vw2,[outpath,'/_group_level/brain_maps/pipe_',PipeStruct_aug.PNAME{1},'/func_pBRAIN_grp.nii']);
        clear Vw Vw2;

        maskisnew=1;
    end

    % creating group-level mean average epi image for QC
    if ~exist( [outpath,'/_group_level/brain_maps/pipe_',PipeStruct_aug.PNAME{1},'/func_tAV_grp.nii'], 'file') % group mean tav epi map
        for ni=1:numel(mask_subj_idxes)
            % check existence of subject specific struct file
            if exist(fullfile( outpath,subject_list{mask_subj_idxes(ni)},'InputStruct_ssa.mat'),'file')
                load( fullfile( outpath,subject_list{mask_subj_idxes(ni)},'InputStruct_ssa.mat') )
            else
                error('cannot find Input struct file for subject: %s \n');
            end
            opath1f = fullfile(outpath,InputStruct_ssa.PREFIX,'func_proc_p1');
            opath2f = fullfile(opath1f,sprintf('subpipe_%03u',PipeStruct_aug.pipe_idx.P1));
            Vw=load_untouch_niiz(sprintf('%s/postwarp/func1_warped_tav.nii.gz',opath2f)); %***%
            if ni==1
                volref = double(Vw.img);
            else
                volref = volref + double(Vw.img);
            end
        end
        volref = volref./numel(mask_subj_idxes);
        Vw.img = volref;
        save_untouch_niiz(Vw,[outpath,'/_group_level/brain_maps/pipe_',PipeStruct_aug.PNAME{1},'/func_tAV_grp.nii']); % group mean sd epi map
        clear Vw;

        maskisnew=1;
    end
    % creating group-level mean temporal SD map (+and variance mask!)
    if ~exist( [outpath,'/_group_level/brain_maps/pipe_',PipeStruct_aug.PNAME{1},'/func_tSD_grp.nii'], 'file') % t-std mask
        for ni=1:numel(mask_subj_idxes)
            % check existence of subject specific struct file
            if exist(fullfile( outpath,subject_list{mask_subj_idxes(ni)},'InputStruct_ssa.mat'),'file')
                load( fullfile( outpath,subject_list{mask_subj_idxes(ni)},'InputStruct_ssa.mat') )
            else
                error('cannot find Input struct file for subject: %s \n');
            end
            opath1f = fullfile(outpath,InputStruct_ssa.PREFIX,'func_proc_p1');
            opath2f = fullfile(opath1f,sprintf('subpipe_%03u',PipeStruct_aug.pipe_idx.P1));
            Vw=load_untouch_niiz(sprintf('%s/postwarp/func1_warped_tsd.nii.gz',opath2f)); %***%
            if ni==1
                volref = double(Vw.img);
            else
                volref = volref + double(Vw.img);
            end
        end
        volref = volref./numel(mask_subj_idxes);
        Vw.img = volref;
        save_untouch_niiz(Vw,[outpath,'/_group_level/brain_maps/pipe_',PipeStruct_aug.PNAME{1},'/func_tSD_grp.nii']); % group mean sd epi map
        clear Vw;

        maskisnew=1;
    end

    % constructing group-level probabilistic tissue maps of CSF, GM, WM
    tisslist = {'CSF','GM','WM'};
    for i=1:3
        if ~exist( [outpath,'/_group_level/brain_maps/pipe_',PipeStruct_aug.PNAME{1},'/func_p',tisslist{i},'_grp.nii'], 'file') % Mean image
        
            for ni=1:numel(mask_subj_idxes)
            
                % check existence of subject specific struct file
                if exist(fullfile( outpath,subject_list{mask_subj_idxes(ni)},'InputStruct_ssa.mat'),'file')  
                    load( fullfile( outpath,subject_list{mask_subj_idxes(ni)},'InputStruct_ssa.mat') )
                else
                    error('cannot find Input struct file for subject: %s \n');
                end
                opath1f = fullfile(outpath,InputStruct_ssa.PREFIX,'func_proc_p1');
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
            Vw.img = volref;
            save_untouch_niiz(Vw,[outpath,'/_group_level/brain_maps/pipe_',PipeStruct_aug.PNAME{1},'/func_p',tisslist{i},'_grp.nii']);
            clear Vw;

            maskisnew=1;
        end
    end

    % ROIMASK is constructing binary functional masks for roi-based regression
    pcsf_file = [outpath,'/_group_level/brain_maps/pipe_',PipeStruct_aug.PNAME{1},'/func_pCSF_grp.nii'];
    pwm_file  = [outpath,'/_group_level/brain_maps/pipe_',PipeStruct_aug.PNAME{1},'/func_pWM_grp.nii'];
    pgm_file  = [outpath,'/_group_level/brain_maps/pipe_',PipeStruct_aug.PNAME{1},'/func_pGM_grp.nii'];
    tsd_file  = [outpath,'/_group_level/brain_maps/pipe_',PipeStruct_aug.PNAME{1},'/func_tSD_grp.nii'];
    bmsk_file = [outpath,'/_group_level/masks/pipe_',PipeStruct_aug.PNAME{1},'/func_brain_mask_grp.nii'];
    
    % constructing binary tissue masks -- for ROI regression (And other things)
    if strcmpi(ParamStruct_aug.ROIMASK{1},'OP1')
        roimask_OP1( pcsf_file,pwm_file,pgm_file,tsd_file,bmsk_file, [outpath,'/_group_level/masks/pipe_',PipeStruct_aug.PNAME{1}], ParamStruct_aug.ROIMASK(2:end) );
    else
        error('unrecognized roimask estimator?!')
    end

    % --> FOR ROIREG... special masks/parcellations, depending on what variant you are using!!!
    
    if strcmpi(PipeStruct_aug.ROIREG{1},'OP1') || strcmpi(PipeStruct_aug.ROIREG{1},'OP2') %% creating PCA-based spatial maps
    
        masklist = {'CSF','WM','tSD'};
        parcdir  = [outpath,'/_group_level/parcellations/pipe_',PipeStruct_aug.PNAME{1}];
        parcpref = [parcdir,'/__opptmp_inter_spwt'];
        if ~exist([parcdir,'/Ugrp_WM.nii'], 'file') || ~exist([parcdir,'/Ugrp_CSF.nii'], 'file') || ~exist([parcdir,'/Ugrp_tSD.nii'], 'file')
            %
            for i=1:numel(masklist)
                mkdir_r(parcpref);
                maskname = [outpath,'/_group_level/masks/pipe_',PipeStruct_aug.PNAME{1},'/func_',masklist{i},'_mask_grp.nii'];
                Mtmp = load_untouch_niiz(maskname);
                ucat=[];
                unix(sprintf('rm %s/volblur.nii',parcpref))
                for ni=1:numel(mask_subj_idxes)
                    [ni numel(mask_subj_idxes)],
                    % check existence of subject specific struct file
                    if exist(fullfile( outpath,subject_list{mask_subj_idxes(ni)},'InputStruct_ssa.mat'),'file')
                        load( fullfile( outpath,subject_list{mask_subj_idxes(ni)},'InputStruct_ssa.mat') )
                    else
                        error('cannot find Input struct file for subject: %s \n');
                    end
                    % NB: only uses run-1 per subject
                    opath1f = fullfile(outpath,InputStruct_ssa.PREFIX,'func_proc_p1');
                    opath2f = fullfile(opath1f,sprintf('subpipe_%03u',PipeStruct_aug.pipe_idx.P1));
                    unix(sprintf('3dBlurInMask -input %s/postwarp/func1_warped.nii.gz -prefix %s/volblur.nii -mask %s -fwhm 6',opath2f,parcpref,maskname));
                    Vtmp = load_untouch_niiz(sprintf('%s/volblur.nii',parcpref));
                    volmat = nifti_to_mat(Vtmp,Mtmp); npc_10 = floor(0.1*size(volmat,2));
                    volmat = bsxfun(@rdivide, bsxfun(@minus,volmat,mean(volmat,2)), std(volmat,0,2)+eps);
                    [u,~,~]=svd( volmat,'econ');
                    ucat = [ucat, u(:,1:npc_10)];
                    npc_10_list(ni,1) = npc_10;
                    unix(sprintf('rm %s/volblur.nii',parcpref))
                end
                [Ugrp,~,~] = svd( ucat,'econ');
                Ugrp = Ugrp(:,1:round(median(npc_10_list)));
                save([parcdir,'/Ugrp_',masklist{i},'.mat'],'Ugrp','npc_10_list');
            
                clear TMPVOL;
                for(p=1:size(Ugrp,2) )
                    tmp=double(Mtmp.img);
                    tmp(tmp>0)= Ugrp(:,p);
                    TMPVOL(:,:,:,p) = tmp;
                end
                nii=Vtmp;
                nii.img = TMPVOL;
                nii.hdr.dime.datatype = 16;
                nii.hdr.hist = Vtmp.hdr.hist;
                nii.hdr.dime.dim(5) = size(Ugrp,2);
                save_untouch_niiz(nii,[parcdir,'/Ugrp_',masklist{i},'.nii']); 
            
                unix(sprintf('rm %s/*.nii', parcpref))
                unix(sprintf('rmdir %s', parcpref))
            end
        end
    
    elseif strcmpi(PipeStruct_aug.ROIREG{1},'OP3') %% resampling pre-made parcellations
        
        roipath = [CODE_PATH,'/reference_maps/roimask']; %% HARD-CODED for now
        e=dir([roipath,'/*.nii']);
        for i=1:numel(e)
            [~,epr,esf]=fileparts(e(i).name);
            if ~strcmp(esf,'.nii')
                error('roimask file does not seem to be in nifti format');
            else
                prefixe = ['roimask_resam_',epr];
            end
            if ~exist(sprintf('%s/%s.nii',[outpath,'/_group_level/parcellations/pipe_',PipeStruct_aug.PNAME{1}],prefixe),'file')
                unix(sprintf('3dresample -master %s -input %s -prefix %s/%s.nii -rmode NN',[outpath,'/_group_level/masks/pipe_',PipeStruct_aug.PNAME{1},'/func_brain_mask_grp.nii'],fullfile(e(i).folder,e(i).name), [outpath,'/_group_level/parcellations/pipe_',PipeStruct_aug.PNAME{1}],prefixe  ));
            end
        end
    elseif strcmpi(PipeStruct_aug.ROIREG{1},'OFF')
        disp('no ROIREG step - no parcellations imported!')
    else
        error('unrecognized ROIREG setting-- not sure how to configure tissue templates!')
    end

else

    disp('using premade/imported Functional group level brain maps n masks');

    % try to pull a represntative first file >> for qc checking
    if exist(fullfile( outpath,subject_list{1},'InputStruct_ssa.mat'),'file')  
        load( fullfile( outpath,subject_list{1},'InputStruct_ssa.mat') )
    else
        error('cannot find Input struct file for subject: %s \n');
    end
    opath1f = fullfile(outpath,InputStruct_ssa.PREFIX,'func_proc_p1');
    opath2f = fullfile(opath1f,sprintf('subpipe_%03u',PipeStruct_aug.pipe_idx.P1));
    Vo=load_untouch_niiz(sprintf('%s/postwarp/func1_warped_tav.nii.gz',opath2f)); %***%
    
    % here is our checklist: subdirs and files
    chklist_a = {'masks','masks','masks','masks','masks','brain_maps','brain_maps'};
    chklist_b = {'func_brain_mask_grp','func_CSF_mask_grp','func_WM_mask_grp','func_tSD_mask_grp','func_GM_mask_grp','func_tAV_grp','func_tSD_grp'};
      
    for k=1:numel(chklist_a)
        fileimp = [outpath,'/_group_level/',chklist_a{k},'/pipe_',PipeStruct_aug.PNAME{1},'/',chklist_b{k},'.nii'];
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

    maskisnew = 1; % for now, default is to redo if you are using an imported mask set
end

% pre-load mask for blck-2
MBstr = ([outpath,'/_group_level/masks/pipe_',PipeStruct_aug.PNAME{1},'/func_brain_mask_grp.nii']);
MB = load_untouch_niiz(MBstr);
maskS = double(MB.img);
% declare empty struct --> nuisance regressor correlation maps
RegCorr2.det=[];
RegCorr2.glb=[];
RegCorr2.mot=[];
RegCorr2.roi=[];

for ns=subj_list_for_proc % step through func-proc (block-2)

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
    opath1f = fullfile(outpath,InputStruct_ssa.PREFIX,'func_proc_p1');
    opath2f = fullfile(opath1f,sprintf('subpipe_%03u',PipeStruct_aug.pipe_idx.P1));
    opath3f = fullfile(opath2f,sprintf('seg_subpipe_%03u',PipeStruct_aug.pipe_idx.Seg ));
    %
    opath4f = fullfile(outpath,InputStruct_ssa.PREFIX,'func_proc_p2',['pipe_',PipeStruct_aug.PNAME{1}]); 

    fprintf('\n===> subj %s (%u/%u), func. processing (part-2)...\n',subject_list{ns},ns,numel(subject_list)),
    %% =======================================================================
    %%      FUNC Processing, Block-2 ...
    %% =======================================================================



    %%%%%========== compatibility adjustments for BLOCK2 ... TASKREG only
    missing_physio=0;
    for nr=1:InputStruct_ssa.N_func
        if isempty(InputStruct_ssa.frun(nr).PHYSIO_filename)
            missing_physio=missing_physio+1;
        end
    end
    if analysis_struct.uses_taskfile==0  && ~strcmpi(PipeStruct_aug.TASKREG{1},'OFF')
        warning('your analysis model does not specify a task design. Turning TASKREG off!')
        PipeStruct_aug.TASKREG{1}={'OFF'};
    end
    %%%%%========== compatibility adjustments, done

    anymiss=0;
    for nr = 1:InputStruct_ssa.N_func
        if ~exist( [opath4f,'/func',num2str(nr),'_fullproc.mat'],'file')
            anymiss=1;
        end
    end

    % appending info about run-length (for regressor construction and interp-contrast-list)
    for nr = 1:InputStruct_ssa.N_func
        %
        opath0   = fullfile(outpath,InputStruct_ssa.PREFIX,'rawdata');
        unix(sprintf('gunzip %s/func%u.nii.gz',opath0,nr));
        hdr      = load_nii_hdr(sprintf('%s/func%u.nii',opath0,nr)); %***%
        unix(sprintf('gzip %s/func%u.nii',opath0,nr));
        InputStruct_ssa.frun(nr).Nt_raw   = hdr.dime.dim(5);
        InputStruct_ssa.frun(nr).Nt_adj   = hdr.dime.dim(5) - InputStruct_ssa.frun(nr).DROP_first - InputStruct_ssa.frun(nr).DROP_last;
    end

    % extract contrasts specified for analysis --> appended into InputStruct, delete unformatted split-info field 
    InputStruct_ssa = interpret_contrast_list( InputStruct_ssa, analysis_struct, ParamStruct_aug.CONTRAST); % generate contrast list for each subject and run    
    %InputStruct_ssa = rmfield(InputStruct_ssa,'task_unformat');
    % extract seed contrasts SAA
    InputStruct_ssa = interpret_seed_list( InputStruct_ssa, analysis_struct, ParamStruct_aug.CONTRAST); % generate contrast list for each subject and run 
    %InputStruct_ssa = rmfield(InputStruct_ssa,'seed_unformat');

    if anymiss==0 && maskisnew==0 % we can skip all of this, if data are processed and masks havent been re-constructed
        disp('Files found, no update to mask(s). Skipping last proc. stage!')
    else
        %--> previously specified opaths over here?
        
        for nr=1:InputStruct_ssa.N_func
    
            if ~exist( [opath4f,'/func',num2str(nr),'_fullproc.mat'],'file') || maskisnew==1 %only do runs with missing outputs / or if mask was just updated
        
                % loading data into mats:
                VSstr = sprintf('%s/postwarp/func%u_warped_smo.nii.gz',opath2f,nr); %***%
                VS = load_untouch_niiz(VSstr);
                volmatS = nifti_to_mat(VS,MB); 
            
                %% first part of regression proc: constructing regressors from DETREG,GSREG,MOTREG,ROIREG,TASKREG 
            
                % >>> Temporal Detrending
                Step = 'DETREG';
                if strcmpi(PipeStruct_aug.(Step){1},'OFF')
                    xdet = ones( InputStruct_ssa.frun(nr).Nt_adj, 1 );
                    statd = [];
                else
                    % get function handle for analysis model of interest
                    currPath=pwd;                               % get current path
                    cd(pipeline_struct.(Step).filepath);               % jump to module directory
                    pfun= str2func(pipeline_struct.(Step).model_name); % get function handle
                    cd(currPath);                               % jump back to current path
                    % execute step:
                    [xdet,statd] = pfun( InputStruct_ssa.frun(nr).Nt_adj, InputStruct_ssa.TR_MSEC, PipeStruct_aug.(Step)(2:end) );  
                end

                % >>> Global Signal Regression
                Step = 'GSREG';
                if strcmpi(PipeStruct_aug.(Step){1},'OFF')
                    xglb = [];
                    statg = [];
                else
                    % get function handle for analysis model of interest
                    currPath=pwd;                               % get current path
                    cd(pipeline_struct.(Step).filepath);               % jump to module directory
                    pfun= str2func(pipeline_struct.(Step).model_name); % get function handle
                    cd(currPath);                               % jump back to current path
                    % execute step:
                    [xglb,statg] = pfun( volmatS, PipeStruct_aug.(Step)(2:end) );  
                end

                % >>> Motion Parameter Regression
                Step = 'MOTREG';
                if strcmpi(PipeStruct_aug.(Step){1},'OFF')
                    xmot = [];
                    statm = [];
                else
                    % get function handle for analysis model of interest
                    currPath=pwd;                               % get current path
                    cd(pipeline_struct.(Step).filepath);               % jump to module directory
                    pfun= str2func(pipeline_struct.(Step).model_name); % get function handle
                    cd(currPath);                               % jump back to current path
                    % execute step:
                    [xmot,statm] = pfun( sprintf('%s/warp/func%u_mpe',opath2f,nr), PipeStruct_aug.(Step)(2:end) );  
                end

                % >>> Noise ROI Regression
                Step = 'ROIREG';
                if strcmpi(PipeStruct_aug.(Step){1},'OFF')
                    xroi = [];
                    statr = [];
                else
                    maskpath  = [outpath,'/_group_level/masks/pipe_',PipeStruct_aug.PNAME{1}];
                    parcpath  = [outpath,'/_group_level/parcellations/pipe_',PipeStruct_aug.PNAME{1}];
                    % get function handle for analysis model of interest
                    currPath=pwd;                               % get current path
                    cd(pipeline_struct.(Step).filepath);               % jump to module directory
                    pfun= str2func(pipeline_struct.(Step).model_name); % get function handle
                    cd(currPath);                               % jump back to current path
                    % execute step:
                    [xroi,statr] = pfun( sprintf('%s/postwarp/func%u_warped.nii.gz',opath2f,nr), {maskpath,parcpath}, PipeStruct_aug.(Step)(2:end) );  
                end
                
                % >>> Regression of Task Design
                Step = 'TASKREG';
                if strcmpi(PipeStruct_aug.(Step){1},'OFF')
                    Xsignal = [];
                elseif strcmpi(PipeStruct_aug.TASKREG{1},'OP1')
                    Xsignal = InputStruct_ssa.task_contrast{nr}.design_mat;
                    Xsignal = bsxfun(@rdivide, Xsignal, sqrt(sum(Xsignal.^2)));
                else
                    error('unknown taskreg option!');
                end

            
                %% second part of regression proc: concat vectors and apply to data matrix 
        
                % --> some regressional stats
                ztmp = zscore(volmatS')';
                if nr==1 && ~isempty(xdet) && size(xdet,2)>1
                    xtmp = zscore(xdet(:,2:end));
                    RegCorr2.det(:,:,ns) = (ztmp*xtmp)./(InputStruct_ssa.frun(nr).Nt_adj-1);
                    RegStat.det(ns,:) = statd;
                end
                if nr==1 && ~isempty(xglb)
                    xtmp = zscore(xglb);
                    RegCorr2.glb(:,:,ns) = (ztmp*xtmp)./(InputStruct_ssa.frun(nr).Nt_adj-1);
                    RegStat.glb(ns,:) = statg;
                end
                if nr==1 && ~isempty(xmot)
                    xtmp = zscore(xmot);
                    RegCorr2.mot(:,:,ns) = (ztmp*xtmp)./(InputStruct_ssa.frun(nr).Nt_adj-1);
                    RegStat.mot(ns,:) = statm;
                end
                if nr==1 && ~isempty(xroi)
                    xtmp = zscore(xroi);
                    RegCorr2.roi(:,:,ns) = (ztmp*xtmp)./(InputStruct_ssa.frun(nr).Nt_adj-1);
                    RegStat.roi(ns,:) = statr;
                end
                % --> some regressional stats
        
                Xnoise = [xdet(:,2:end),xglb,xmot,xroi]; % full noise matrix, normed to unit length
                Xnoise = bsxfun(@rdivide,Xnoise,sqrt(sum(Xnoise.^2)));
                [ output ] = GLM_model_fmri( volmatS, [0], [Xnoise], [Xsignal], 1, 1 );
                volmatF = output.vol_denoi;
            
                % >>> Low-pass Filtering
                Step = 'LOPASS';
                if strcmpi(PipeStruct_aug.(Step){1},'OFF')
                    disp('skipping lopass');
                else
                    % get function handle for analysis model of interest
                    currPath=pwd;                               % get current path
                    cd(pipeline_struct.(Step).filepath);               % jump to module directory
                    pfun= str2func(pipeline_struct.(Step).model_name); % get function handle
                    cd(currPath);                               % jump back to current path
                    % execute step:
                    volmatF = pfun( volmatF, InputStruct_ssa.TR_MSEC, PipeStruct_aug.(Step)(2:end) );  
                end
                
                % >>> Component-based filtering
                if strcmpi(PipeStruct_aug.COMPFILT{1},'OFF')
                    disp('skipping compfilt');
                else
                    error('no compfilt options enabled yet!!')
                end
            
                %% Export fully processed data as .mat file and full .nii file
                clear TMPVOL;
                for(p=1:size(volmatF,2) )
                    tmp=double(MB.img);
                    tmp(tmp>0)= volmatF(:,p);
                    TMPVOL(:,:,:,p) = tmp;
                end
                nii=VS;
                nii.img = TMPVOL;
                nii.hdr.dime.datatype = 16;
                nii.hdr.hist = VS.hdr.hist;
                nii.hdr.dime.dim(5) = size(volmatF,2);
                save_untouch_niiz(nii,[opath4f,'/func',num2str(nr),'_fullproc.nii.gz']); 
                save([opath4f,'/func',num2str(nr),'_fullproc.mat'],'volmatF');
            end
        end
    end


    if strcmpi( ParamStruct_aug.ANALYSIS, 'NONE')
        disp('no further analysis, just doing one last quick check...')
        % extract and put the runs into cell array
        clear volcel;
        for nr=1:InputStruct_ssa.N_func
           x=load([opath4f,'/func',num2str(nr),'_fullproc.mat']);
           % quick check in case mask somehow doesnt match matfile
           if size(x.volmatF,1) ~= sum(maskS(:))
                error('fully processed matfile in %s does not match functional mask! Try deleting group level folders and rerunning P2!',opath4f)
           end
        end
        disp('ok!');

    else
        disp('now doing subject-level analysis...')

        % construct relevant path, load runfiles
        opath5f = fullfile(opath4f, ParamStruct_aug.Variable_ID);
        mkdir_r(opath5f);

        if ~isempty(ParamStruct_aug.GMMASK_ANL)
            if strcmpi(ParamStruct_aug.GMMASK_ANL,'AUTO')
                gmmask_anl_file = [outpath,'/_group_level/masks/pipe_',PipeStruct_aug.PNAME{1},'/func_GM_mask_grp.nii'];
            else
                if exist(ParamStruct_aug.GMMASK_ANL,'file')
                    gmmask_anl_file = ParamStruct_aug.GMMASK_ANL;
                else
                    error('analysis grey matter mask %s cannot be located?',ParamStruct_aug.GMMASK_ANL)
                end
            end
            if  contains(gmmask_anl_file,'.nii')
                unix(sprintf('cp %s %s/gmmask_anl_pre.nii',gmmask_anl_file,opath5f));
                unix(sprintf('gzip %s/gmmask_anl_pre.nii',opath5f))
            elseif contains(gmmask_anl_file,'.nii.gz')
                unix(sprintf('cp %s %s/gmmask_anl_pre.nii',gmmask_anl_file,opath5f));
            else
                error('non-nifti format of analysis grey matter mask')
            end
            VGA = load_untouch_niiz(gmmask_anl_file);
            VGA.img = VGA.img .* maskS;
            save_untouch_nii(VGA,sprintf('%s/gmmask_anl.nii',opath5f));
            unix(sprintf('rm %s/gmmask_anl_pre.nii',opath5f));
            vga_vec = nifti_to_mat(VGA,MB);
            kepix_vga = find(vga_vec>0);
            maskSsub = VGA.img;
        else
            kepix_vga = 1:sum(maskS(:));
            maskSsub = maskS;
        end
        
        % extract and put the runs into cell array
        clear volcel;
        for nr=1:InputStruct_ssa.N_func
           x=load([opath4f,'/func',num2str(nr),'_fullproc.mat']);

           % quick check in case mask somehow doesnt match matfile
           if size(x.volmatF,1) ~= sum(maskS(:))
                error('fully processed matfile in %s does not match functional mask! Try deleting group level folders and rerunning P2!',opath4f)
           end

           volcel{nr} = x.volmatF(kepix_vga,:);
        end

        % if seed-based analysis, make sure to produce appropriately resampled copies.....
        if analysis_struct.uses_roifile>0
            seedpath = fullfile(outpath,InputStruct_ssa.PREFIX,'func_seeds',sprintf('subpipe_%03u',PipeStruct_aug.pipe_idx.P1));
            for i=1:numel(InputStruct_ssa.seed_info.contrast)
                if ~exist(sprintf('%s/%s.resam.nii.gz',seedpath,InputStruct_ssa.seed_info.contrast(i).label),'file')

                    % copy over base mask and gzip if unzipped
                    if contains(InputStruct_ssa.seed_info.contrast(i).location,'.nii.gz')
                        unix(sprintf('cp %s %s/%s.nii.gz',InputStruct_ssa.seed_info.contrast(i).location, seedpath, InputStruct_ssa.seed_info.contrast(i).label))
                    elseif contains(InputStruct_ssa.seed_info.contrast(i).location,'.nii')
                        unix(sprintf('cp %s %s/%s.nii',InputStruct_ssa.seed_info.contrast(i).location, seedpath, InputStruct_ssa.seed_info.contrast(i).label));
                        unix(sprintf('gzip %s/%s.nii', seedpath, InputStruct_ssa.seed_info.contrast(i).label));
                    end
                    unix(sprintf('3dresample -master %s/postwarp/func1_warped_mask_clean.nii.gz -input %s/%s.nii.gz -prefix %s/%s.resam.nii.gz',opath2f, seedpath,InputStruct_ssa.seed_info.contrast(i).label, seedpath,InputStruct_ssa.seed_info.contrast(i).label));
                else
                    disp('seedmap found! skipping to next');
                end
                % loading data into mats:
                VSstr = sprintf('%s/%s.resam.nii.gz',seedpath, InputStruct_ssa.seed_info.contrast(i).label); %***%
                VS = load_untouch_niiz(VSstr);
                stmp = nifti_to_mat(VS,MB); 
                if size(stmp,2)>1
                    error('seed nifti files must be 3D (single-volume)!')
                end
                InputStruct_ssa.seed_info.seedmat(:,i) = stmp(kepix_vga); 
            end
        end

        % get function handle for analysis model of interest
        currPath=pwd;                               % get current path
        cd(analysis_struct.filepath);               % jump to module directory
        pfun= str2func(analysis_struct.model_name); % get function handle
        cd(currPath);                               % jump back to current path
        
        % call function
        ParamStruct_aug.TR_MSEC = InputStruct_ssa.TR_MSEC;
        ParamStruct_aug.mask    = maskSsub;
        %
        if analysis_struct.uses_taskfile>0
            out_analysis = pfun( volcel, InputStruct_ssa.task_info, ParamStruct_aug );  
        elseif analysis_struct.uses_roifile>0
            out_analysis = pfun( volcel, InputStruct_ssa.seed_info, ParamStruct_aug );  
        else
            out_analysis = pfun( volcel, ParamStruct_aug );  
        end
        % store submask:
        out_analysis.submask = maskSsub;
        %
        save([opath5f,'/out_analysis.mat'],'out_analysis');

        fi = fieldnames( out_analysis.image );

        for i=1:numel(fi)
            clear TMPVOL;
            for p=1:size(out_analysis.image.(fi{i}),2)
                tmp=double(maskS);
                tmp(tmp>0)= out_analysis.image.(fi{i})(:,p);
                TMPVOL(:,:,:,p) = tmp;
            end
            nii=MB;
            nii.img = TMPVOL;
            nii.hdr.dime.datatype = 16;
            nii.hdr.hist = MB.hdr.hist;
            nii.hdr.dime.dim(5) = size(out_analysis.image.(fi{i}),2);
            save_untouch_niiz(nii,[opath5f,'/',fi{i},'.nii.gz']); 
        end
    end
end

save([outpath,'/_group_level/brain_maps/pipe_',PipeStruct_aug.PNAME{1},'/regstat.mat'],'RegCorr2');

disp('funxionale block-2 done');

