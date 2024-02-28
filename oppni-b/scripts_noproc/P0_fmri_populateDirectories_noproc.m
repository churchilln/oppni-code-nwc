function [subject_list, InputStruct_aug, PipeStruct_aug, ParamStruct_aug] = P0_fmri_populateDirectories_noproc( inputfile, paramlist, outpath, push_overwrite )
%
% =========================================================================
% P0_FMRI_POPULATEDIRECTORIES: this script should be run before any other steps
% in bold pipeline. It does the following:
% . (a) checks if input/pipelin/param/template/task/seed files can be read 
% . (b) generates output directory structure, 
% . (c) creates dir-specific InputStructs to protect against overwriting
% =========================================================================
%
% Syntax:
%
%      P0_fmri_populateDirectories( inputfile, pipelinefile, paramlist, outpath, push_overwrite )
%
% Inputs:
%      inputfile : string, giving name of input textfile listing data to process
%      paramlist : string, giving name of parameter file specifying special arguments to use in pipelines
%      outpath : string, specifying destination directory for processed outputs
%      push_overwrite : optional binary argument. If =0, and your input file doesn't match what was already processed, throws error.
%                       If =1, it will overwrite input file information and suppress this error. USE WITH CAUTION!!
%

% declaring path
CODE_PATH = fileparts(which('P0_fmri_populateDirectories.m'));
if CODE_PATH(end)~='/'
    CODE_PATH = [CODE_PATH '/'];
end

if nargin<5
    push_overwrite=0;
end

%% Preliminary check of files

outpath = fullfile(outpath,'fmri_proc'); % subdir should be fmri_proc

% basic file checks ... construct input/pipeline/param structure -> should throw error if inputs non-valid
InputStruct = Read_Input_File_fmri_noproc(inputfile);
ParamStruct = Read_Params_File_fmri(paramlist);

%% Construct group-level directory structure, with template file copying

% creating group-level directory structure
mkdir_r( fullfile(outpath,'_group_level',{'masks','brain_maps','parcellations'}) );
mkdir_r( fullfile(outpath,'_group_level','QC',{'qc.compat','qc.quant'}) );

% also constructing output path for contrast+analysis
ParamStruct.Fixed_ID = 'NULL';
% and the output path for analysis/contrast combos
astr = ParamStruct.ANALYSIS;
cstr = ParamStruct.CONTRAST;
if( strcmpi(cstr,'NONE') || strcmpi(cstr,'0') ) cstr = 'no_contrast'; end
if( strcmpi(astr,'NONE') || strcmpi(astr,'0') ) astr = 'no_analysis'; end
ParamStruct.Variable_ID = [cstr,'.',astr];
% ouput "augmented" param structure
ParamStruct_aug = ParamStruct; 
clear ParamStruct;

% pipeline defaults
PipeStruct_aug.PNAME{1}      = 'noproc';
PipeStruct_aug.pipe_idx.P1   = 1;
PipeStruct_aug.pipe_idx.Seg  = 1;
PipeStruct_aug.pipe_idx.Warp = 1;

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
    mkdir_r( fullfile( outpath,InputStruct_temp.PREFIX, {'rawdata','func_seeds','func_contrasts','anat_proc','func_proc_p1','func_proc_p2','phys_proc'}) );
    mkdir_r( fullfile( outpath,InputStruct_temp.PREFIX,'func_proc_p1','init_mot_estim') );
    
    % and under specific pipeline combos
    mkdir_r( fullfile( outpath,InputStruct_temp.PREFIX,'anat_proc',sprintf('subpipe_%03u',PipeStruct_aug.pipe_idx.Warp),sprintf('seg_subpipe_%03u',PipeStruct_aug.pipe_idx.Seg) ) );
    mkdir_r( fullfile( outpath,InputStruct_temp.PREFIX,'func_proc_p1',sprintf('subpipe_%03u',PipeStruct_aug.pipe_idx.P1),{'prewarp','warp','postwarp',sprintf('seg_subpipe_%03u',PipeStruct_aug.pipe_idx.Seg)} ) );
    mkdir_r( fullfile( outpath,InputStruct_temp.PREFIX,'func_proc_p2',sprintf('pipe_%s',PipeStruct_aug.PNAME{1}) ) );
    % special subdir for resliced func-seeds ... can depend on alignment steps
    mkdir_r( fullfile( outpath,InputStruct_temp.PREFIX,'func_seeds',sprintf('subpipe_%03u',PipeStruct_aug.pipe_idx.P1)) );

    % append param. list
    InputStruct_temp.Param_Fixed_ID = ParamStruct_aug.Fixed_ID;
    % append task design information (sans contrast)
    for nr=1:InputStruct_temp.N_func
        [InputStruct_temp.task_unformat{nr}] = Read_Task_File_fmri(  InputStruct_temp.frun(nr).TASK_filename  );
    end
    % append seed information 
    [InputStruct_temp.seed_unformat] = Read_Seed_File_fmri(  InputStruct_temp.SEED_filename  );
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
            for nr=1:InputStruct_temp.N_func
                if ~isequal(InputStruct_ssa.task_unformat{nr},InputStruct_temp.task_unformat{nr})
                    errstr=sprintf('%s\n\tYour task design file (run %u) has changed!',errstr,nr);
                end
            end
            if ~isequal(InputStruct_ssa.seed_unformat,InputStruct_temp.seed_unformat)
                errstr=sprintf('%s\n\tYour seed info file has changed!',errstr);
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
