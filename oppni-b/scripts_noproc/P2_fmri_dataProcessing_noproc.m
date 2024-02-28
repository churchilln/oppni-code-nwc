function P2_fmri_dataProcessing_noproc( inputfile, paramlist, outpath, big_skip, mask_subj_idxes, subj_subset )
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
[subject_list, InputStruct_aug, PipeStruct_aug, ParamStruct_aug] = P0_fmri_populateDirectories_noproc( inputfile, paramlist, outpath );

% now augmenting outpath... do this after P0!
outpath = fullfile(outpath,'fmri_proc'); % subdir should be fmri_proc
% check for missing files from previous step too
File_Existence_Checker_fmri_noproc(InputStruct_aug,outpath,1); 

% check validity of analysis model etc. --> carry forward information about the analysis too!
analysis_struct = check_fmri_analysis_model( ParamStruct_aug.ANALYSIS );
% also, modify "number of components" field if more than one contrast is specified
if( ~isempty(strfind(ParamStruct_aug.CONTRAST,',' )) )
    analysis_struct.num_comp = 'multi_component';
end

% Now configuring param defaults if unspecified
%
if ~isfield(ParamStruct_aug,'INIMOT')
    ParamStruct_aug.INIMOT='OP1';
end
if ~isfield(ParamStruct_aug,'ROIMASK')
    ParamStruct_aug.ROIMASK='OP1';
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

    fprintf('\n===> subj %s (%u/%u), physio. processing...\n',subject_list{ns},ns,numel(subject_list)),
    %% =======================================================================
    %%      PHYSIO Processing ...
    %% =======================================================================

    fprintf('\n===> phys-proc. now on subj %u/%u: %s...\n',ns,numel(subject_list),subject_list{ns}),
    
    disp('nothin for physio-proc so far. put in soon!')

    fprintf('\n===> subj %s (%u/%u), func. processing (part-1)...\n',subject_list{ns},ns,numel(subject_list)),
    %% =======================================================================
    %%      FUNC Processing ...
    %% =======================================================================

    % port over the functional data
    for nr=1:InputStruct_ssa.N_func
        unix(sprintf('cp %s/func%u.nii.gz %s/func%u_fullproc.nii.gz',opath0,nr,opath4f,nr))
        if ~isempty(InputStruct_ssa.frun(nr).FMASK_filename)
        unix(sprintf('cp %s/func%u_mask.nii.gz %s/postwarp/func%u_warped_mask_clean.nii.gz',opath0,nr,opath2f,nr))
        end
        if ~isempty(InputStruct_ssa.frun(nr).MPE_filename)
        unix(sprintf('cp %s/func%u_mpe %s/func%u_mpe',opath0,nr,opath2f,nr))
        end
    end

    % extra step: make mask, mean, sd maps --> run-1 mean and sd used for creating group-level maps. For runs >1, just kept for qc purposes
    for nr=1:InputStruct_ssa.N_func
        if ~exist(sprintf('%s/postwarp/func%u_warped_mask_clean.nii.gz',opath2f,nr),'file') || ...
           ~exist(sprintf('%s/postwarp/func%u_warped_tav.nii.gz',opath2f,nr),'file')  || ...
           ~exist(sprintf('%s/postwarp/func%u_warped_tsd.nii.gz',opath2f,nr),'file') 

            unix(sprintf('3dAutomask -prefix %s/postwarp/func%u_warped_mask_clean.nii.gz %s/func%u_fullproc.nii.gz',opath2f,nr,opath4f,nr));
            unix(sprintf('3dTstat -mean  -prefix %s/postwarp/func%u_warped_tav.nii.gz %s/func%u_fullproc.nii.gz',opath2f,nr,opath4f,nr));
            unix(sprintf('3dTstat -stdev -prefix %s/postwarp/func%u_warped_tsd.nii.gz %s/func%u_fullproc.nii.gz',opath2f,nr,opath4f,nr));
        else
            disp('func mask done - skipping ahead.')
        end
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
    chklist_a = {'masks','brain_maps','brain_maps'};
    chklist_b = {'func_brain_mask_grp','func_tAV_grp','func_tSD_grp'};
      
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
        InputStruct_ssa.frun(nr).Nt_raw   = hdr.dime.dim(5) + InputStruct_ssa.frun(nr).DROP_first + InputStruct_ssa.frun(nr).DROP_last;
        InputStruct_ssa.frun(nr).Nt_adj   = hdr.dime.dim(5);
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
    
            if ~exist( [opath4f,'/func',num2str(nr),'_fullproc.mat'],'file') %only do runs with missing outputs
        
                % loading data into mats:
                VSstr = sprintf('%s/func%u_fullproc.nii.gz',opath4f,nr); %***%
                VS = load_untouch_niiz(VSstr);
                volmatF = nifti_to_mat(VS,MB); 
                save([opath4f,'/func',num2str(nr),'_fullproc.mat'],'volmatF');
            end
        end
    end


    if ~strcmpi( ParamStruct_aug.ANALYSIS, 'NONE')

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

disp('funxionale block-2 done');
