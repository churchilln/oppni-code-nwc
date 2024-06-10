function [subject_list, InputStruct_aug, PipeStruct_aug, ParamStruct_aug] = P0_perf_populateDirectories( inputfile, pipelinefile, paramlist, outpath, push_overwrite )
%
% =========================================================================
% P0_PERF_POPULATEDIRECTORIES: this script should be run before any other steps
% in bold pipeline. It does the following:
% . (a) checks if input/pipelin/param/template/task/seed files can be read 
% . (b) generates output directory structure, 
% . (c) creates dir-specific InputStructs to protect against overwriting
% =========================================================================
%
% Syntax:
%
%      P0_perf_populateDirectories( inputfile, pipelinefile, paramlist, outpath, push_overwrite )
%
% Inputs:
%      inputfile : string, giving name of input textfile listing data to process
%      pipelinefiles : string, giving name of pipeline textfile specifying which pipeline steps to apply to the data
%      paramlist : string, giving name of parameter file specifying special arguments to use in pipelines
%      outpath : string, specifying destination directory for processed outputs
%      push_overwrite : optional binary argument. If =0, and your input file doesn't match what was already processed, throws error.
%                       If =1, it will overwrite input file information and suppress this error. USE WITH CAUTION!!
%

% declaring path
CODE_PATH = fileparts(which('P0_perf_populateDirectories.m'));
if CODE_PATH(end)~='/'
    CODE_PATH = [CODE_PATH '/'];
end

if nargin<5
    push_overwrite=0;
end

%% Preliminary check of files

outpath = fullfile(outpath,'perf_proc'); % subdir should be fmri_proc
mkdir_r(outpath); % construct path
e=dir([outpath,'*']); % dir to get absolute
outpath = fullfile(e.folder,e.name); % convert to absolute

% basic file checks ... construct input/pipeline/param structure -> should throw error if inputs non-valid
InputStruct = Read_Input_File_perf(inputfile);
PipeStruct  = Read_Pipeline_File_perf(pipelinefile);
ParamStruct = Read_Params_File_perf(paramlist);
if ~exist( ParamStruct.TEMPLATE, 'file')
    error('cannot find TEMPLATE image specified in paramfile:\n\t%s\n');
end

%% Construct group-level directory structure, with template file copying

% creating group-level directory structure
mkdir_r( fullfile(outpath,'_group_level',{'masks','brain_maps','parcellations','templates'}) );
mkdir_r( fullfile(outpath,'_group_level','QC',{'qc.compat','qc.manual','qc.quant'}) );

% template checking
[~,refp, refe]=fileparts( ParamStruct.TEMPLATE );
[~,refp2,~]=fileparts(refp); % make sure any ext is stripped
% parsing and remaking template info.
ParamStruct.TEMPLATE_loc    = fullfile(outpath,'_group_level','templates',[refp,refe]);
ParamStruct.TEMPLATE_prefix = refp2;
ParamStruct.TEMPLATE_orig   = ParamStruct.TEMPLATE;
ParamStruct = rmfield(ParamStruct,'TEMPLATE');

if exist(ParamStruct.TEMPLATE_loc,'file')

    Vtemp_new = load_untouch_niiz(ParamStruct.TEMPLATE_orig);
    Vtemp_old = load_untouch_niiz(ParamStruct.TEMPLATE_loc);
    
    if  max(abs(Vtemp_new.img(:)-Vtemp_old.img(:)))>1E-9
        error('template: %s is same name as pre-existing one, but image data are different!',ParamStruct.TEMPLATE_orig);
    elseif(~isequal(Vtemp_new.hdr,Vtemp_old.hdr)) 
        error('template: %s is same name as pre-existing one, image matched but headers are different!',ParamStruct.TEMPLATE_orig); 
    else
        disp('template: matches existing version')
    end
else 
    %%% otherwise keep a copy in group folder for provenance tracking
    unix(sprintf('cp %s %s',ParamStruct.TEMPLATE_orig,ParamStruct.TEMPLATE_loc));
end

% also constructing output path for contrast+analysis
ParamStruct.Fixed_ID = [ParamStruct.TEMPLATE_prefix,'.',ParamStruct.VOXRES,'mm'];
% and the output path for analysis/contrast combos
pstr = ParamStruct.PERF_MODEL;
astr = ParamStruct.ANALYSIS;
cstr = ParamStruct.CONTRAST;
if( strcmpi(pstr,'NONE') || strcmpi(pstr,'0') ) pstr = 'no_perf_model'; end
if( strcmpi(cstr,'NONE') || strcmpi(cstr,'0') ) cstr = 'no_contrast'; end
if( strcmpi(astr,'NONE') || strcmpi(astr,'0') ) astr = 'no_analysis'; end
ParamStruct.Variable_ID = [pstr,'.',cstr,'.',astr];
% ouput "augmented" param structure
ParamStruct_aug = ParamStruct; 
clear ParamStruct;

%% Pipeline manager step

mkdir_r( fullfile(outpath,'_pipe_manager') );

if ~exist(fullfile(outpath,'_pipe_manager','pipe_key.mat'),'file')

    pipe_key.Warp{1,1} = 1;
    pipe_key.Warp{1,2} = PipeStruct.Warp_ID;
    pipe_key.Seg{1,1}  = 1;
    pipe_key.Seg{1,2}  = PipeStruct.Seg_ID;
    pipe_key.P1{1,1}   = 1;
    pipe_key.P1{1,2}   = PipeStruct.P1_ID;
    pipe_key.P2{1,1} = PipeStruct.PNAME{1};
    pipe_key.P2{1,2} = PipeStruct.P2_ID;

    pipe_idx.Warp = 1;
    pipe_idx.Seg  = 1;
    pipe_idx.P1   = 1;

    save(fullfile(outpath,'_pipe_manager','pipe_key.mat'),'pipe_key');
else
    load( fullfile(outpath,'_pipe_manager','pipe_key.mat') );
    % modular blocks

    % search warp
    ix = find( strcmp(PipeStruct.Warp_ID,pipe_key.Warp(:,2)));
    if ~isempty(ix)
        pipe_idx.Warp = pipe_key.Warp{ix,1}; % allocate correct pipe ID#
    else
        pid = cell2mat(pipe_key.Warp(:,1)); % pull list of pipe ID#s
        oi = setdiff( 1:(max(pid)+1), pid ); % get smallest possible pipe ID# that is unclaimed
        oi = oi(1);
        pipe_idx.Warp = oi(1);  % allocate correct pipe ID#
        iqq = numel(pipe_key.Warp(:,1))+1; % increment size of pipekey array
        pipe_key.Warp{iqq,1} = oi(1);
        pipe_key.Warp{iqq,2} = PipeStruct.Warp_ID;
    end

    % search seg
    ix = find( strcmp(PipeStruct.Seg_ID,pipe_key.Seg(:,2)));
    if ~isempty(ix)
        pipe_idx.Seg = pipe_key.Seg{ix,1};
    else
        pid = cell2mat(pipe_key.Seg(:,1)); % pull list of pipe ID#s
        oi = setdiff( 1:(max(pid)+1), pid ); % get smallest possible pipe ID# that is unclaimed
        oi = oi(1);
        pipe_idx.Seg = oi(1);  % allocate correct pipe ID#
        iqq = numel(pipe_key.Seg(:,1))+1; % increment size of pipekey array
        pipe_key.Seg{iqq,1} = oi(1);
        pipe_key.Seg{iqq,2} = PipeStruct.Seg_ID;
    end

    % search P1
    ix = find( strcmp(PipeStruct.P1_ID,pipe_key.P1(:,2)));
    if ~isempty(ix)
        pipe_idx.P1 = pipe_key.P1{ix,1};
    else
        pid = cell2mat(pipe_key.P1(:,1)); % pull list of pipe ID#s
        oi = setdiff( 1:(max(pid)+1), pid ); % get smallest possible pipe ID# that is unclaimed
        oi = oi(1);
        pipe_idx.P1 = oi(1);  % allocate correct pipe ID#
        iqq = numel(pipe_key.P1(:,1))+1; % increment size of pipekey array
        pipe_key.P1{iqq,1} = oi(1);
        pipe_key.P1{iqq,2} = PipeStruct.P1_ID;
    end

    % search P2
    ix = find( strcmp(PipeStruct.PNAME{1},pipe_key.P2(:,1)));
    if ~isempty(ix)
        % "authenticate" -> match was found for pname, make sure if corresponds to the same actual sequence
        if ~strcmp(PipeStruct.P2_ID,pipe_key.P2{ix,2})
            error('Pipeline with name %s does not match existing one! What changed?\n\tOld: %s\n\tNew: %s\n',PipeStruct.PNAME{1},pipe_key.P2{ix,2},PipeStruct.P2_ID)
        end
    else
        % if found to be a unique pipeline name, double-check it wasnt already constructed under different name!
        ix2 = find( strcmp(PipeStruct.P2_ID,pipe_key.P2(:,2)) );
        if ~isempty(ix2)
            error('This pipeline %s was already specified under a different name of %s. Stopping because it is redundant.\n\tOld: %s\n\tNew: %s\n',PipeStruct.PNAME{1},pipe_key.P2{ix2,1},pipe_key.P2{ix2,2},PipeStruct.P2_ID)
        end

        tmpix= numel(pipe_key.P2(:,1))+1;
        pipe_key.P2{tmpix,1} = PipeStruct.PNAME{1};
        pipe_key.P2{tmpix,2} = PipeStruct.P2_ID;
    end
    
    save(fullfile(outpath,'_pipe_manager','pipe_key.mat'),'pipe_key');
end

% store numeric indexing for pipeline structures
PipeStruct.pipe_idx = pipe_idx;

% just using same naming convention as other structs for now
PipeStruct_aug = PipeStruct; clear PipeStruct;

%% Construct group-level pipeline directory structure

mkdir( fullfile(outpath,'_group_level','masks',sprintf('pipe_%s',PipeStruct_aug.PNAME{1})) )
mkdir( fullfile(outpath,'_group_level','brain_maps',sprintf('pipe_%s',PipeStruct_aug.PNAME{1})) )
mkdir( fullfile(outpath,'_group_level','parcellations',sprintf('pipe_%s',PipeStruct_aug.PNAME{1})) )

%% Construct subject-level pipeline directory structure

% now store all files to output
for(ns=1:numel(InputStruct))

    % take ith subject input structure
    InputStruct_temp = InputStruct(ns);
    % now reconstruct full folder hierarchy for this subject
    mkdir_r( fullfile( outpath,InputStruct_temp.PREFIX, {'rawdata','perf_contrasts','anat_proc','perf_proc_p1','perf_proc_p2'}) );
    mkdir_r( fullfile( outpath,InputStruct_temp.PREFIX,'perf_proc_p1','init_mot_estim') );
    
    % and under specific pipeline combos
    mkdir_r( fullfile( outpath,InputStruct_temp.PREFIX,'anat_proc',sprintf('subpipe_%03u',PipeStruct_aug.pipe_idx.Warp),sprintf('seg_subpipe_%03u',PipeStruct_aug.pipe_idx.Seg) ) );
    mkdir_r( fullfile( outpath,InputStruct_temp.PREFIX,'perf_proc_p1',sprintf('subpipe_%03u',PipeStruct_aug.pipe_idx.P1),{'prewarp','warp','postwarp',sprintf('seg_subpipe_%03u',PipeStruct_aug.pipe_idx.Seg)} ) );
    mkdir_r( fullfile( outpath,InputStruct_temp.PREFIX,'perf_proc_p2',sprintf('pipe_%s',PipeStruct_aug.PNAME{1}) ) );

    % append param. list
    InputStruct_temp.Param_Fixed_ID = ParamStruct_aug.Fixed_ID;
    % append task design information (sans contrast)
    for nr=1:InputStruct_temp.N_perf
        [InputStruct_temp.task_unformat{nr}] = Read_Task_File_perf(  InputStruct_temp.prun(nr).TASK_filename  );
    end
    % this gets put in augmented inputstruct
    InputStruct_aug(ns) = InputStruct_temp;

    if exist(fullfile( outpath,InputStruct_temp.PREFIX,'InputStruct_ssa.mat'),'file') && push_overwrite==0

        %%% if previously existent file - only continue if params are same...and keep the older version
        load(fullfile( outpath,InputStruct_temp.PREFIX,'InputStruct_ssa.mat'));

        if ~isequal(InputStruct_temp,InputStruct_ssa)

            errstr = sprintf('Input struct for "%s" does not match the one you used to previously run these data...',InputStruct_ssa.PREFIX);
            if ~isequal(InputStruct_ssa.Param_Fixed_ID,InputStruct_temp.Param_Fixed_ID)
                errstr=sprintf('%s\n\tYour param file''s fixed values have changed!',errstr);
            end
            for nr=1:InputStruct_temp.N_perf
                if ~isequal(InputStruct_ssa.task_unformat{nr},InputStruct_temp.task_unformat{nr})
                    errstr=sprintf('%s\n\tYour task design file (run %u) has changed!',errstr,nr);
                end
            end
            error(errstr);
        else
            disp('Input struct unchanged from last submission. Proceeding to next step...');
        end
    else
        disp('Input struct being created...');
        %%% redefine InputStruct as "temp" version 
        InputStruct_ssa = InputStruct_temp; clear InputStruct_temp;
        save( fullfile( outpath,InputStruct_ssa.PREFIX,'InputStruct_ssa.mat'),'InputStruct_ssa' );
    end
    subject_list{ns} = InputStruct(ns).PREFIX;
end

% store "augmented" input structre for reference
save( fullfile(outpath,'_group_level','InputStruct_aug.mat'),'InputStruct_aug');

disp('Step-0 Complete!');
