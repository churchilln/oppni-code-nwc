function P0_fmri_removeDirectories( inputfile, pipelinefile, paramlist, outpath, PartToDelete )
%
% . this script: 
% . strips out a specific pipeline

% declaring path
CODE_PATH = fileparts(which('P0_fmri_populateDirectories.m'));
if CODE_PATH(end)~='/'
    CODE_PATH = [CODE_PATH '/'];
end

%% Preliminary check of files

outpath = fullfile(outpath,'fmri_proc'); % subdir should be fmri_proc
if ~exist(outpath,'dir') error('fmri proc directory does not exist!'); end
e=dir([outpath,'*']); % dir to get absolute
outpath = fullfile(e.folder,e.name); % convert to absolute

% basic file checks ... construct input/pipeline/param structure -> should throw error if inputs non-valid
InputStruct = Read_Input_File_fmri(inputfile);
PipeStruct  = Read_Pipeline_File_fmri(pipelinefile);
ParamStruct = Read_Params_File_fmri(paramlist);
if ~exist( ParamStruct.TEMPLATE, 'file')
    error('cannot find TEMPLATE image specified in paramfile:\n\t%s\n');
end

%% Verify group-level directory structure, with template file checking

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
    error('specific template does not exist?')
end

% still constructing output path for contrast+analysis --> ** strip option here to be added ** 
ParamStruct.Fixed_ID = [ParamStruct.TEMPLATE_prefix,'.',ParamStruct.VOXRES,'mm'];
% and the output path for analysis/contrast combos
astr = ParamStruct.ANALYSIS;
cstr = ParamStruct.CONTRAST;
if( strcmpi(cstr,'NONE') || strcmpi(cstr,'0') ) cstr = 'no_contrast'; end
if( strcmpi(astr,'NONE') || strcmpi(astr,'0') ) astr = 'no_analysis'; end
ParamStruct.Variable_ID = [cstr,'.',astr];
% ouput "augmented" param structure
ParamStruct_aug = ParamStruct; 
clear ParamStruct;

%% Pipeline manager step

if ~exist(fullfile(outpath,'_pipe_manager','pipe_key.mat'),'file')
    error('no pipeline key found - no way to safely delete!')
else
    load( fullfile(outpath,'_pipe_manager','pipe_key.mat') );

    % where in "branch" to delete from
    if strcmpi(PartToDelete,'Warp') % search warp
        ix = find( strcmp(PipeStruct.Warp_ID,pipe_key.Warp(:,2)));
        if ~isempty(ix)
            %> all descend depend on Warp
            cell_idx_todel.Warp = ix; % row of cell array to del
            pipe_idx_todel.Warp = pipe_key.Warp(cell_idx_todel.Warp,1); % name of folder todel
            % 1st descend - all segs containing warp
            cell_idx_todel.Seg  = find( contains( pipe_key.Seg(:,2), PipeStruct.Warp_ID ) ); %
            pipe_idx_todel.Seg  = pipe_key.Seg( cell_idx_todel.Seg, 1 );
            % 2nd descend - all p1s containing warp
            cell_idx_todel.P1  = find( contains( pipe_key.P1(:,2), PipeStruct.Warp_ID ) ); %
            pipe_idx_todel.P1  = pipe_key.P1( cell_idx_todel.P1, 1 );
            % 3rd descend - all p2s containing warp
            cell_idx_todel.P2  = find( contains( pipe_key.P2(:,2), PipeStruct.Warp_ID ) ); %
            pipe_idx_todel.P2  = pipe_key.P2( cell_idx_todel.P2, 1 );
        else
            error('pipeline block Warp not found -- does this pipeline exist?')
        end
    elseif strcmpi(PartToDelete,'Seg') % search seg
        ix = find( strcmp(PipeStruct.Seg_ID,pipe_key.Seg(:,2)));
        if ~isempty(ix)
            %> Only P2 descend depend on Seg
            cell_idx_todel.Warp = [];
            pipe_idx_todel.Warp = [];
            %
            cell_idx_todel.P1   = [];
            pipe_idx_todel.P1   = [];
            %
            cell_idx_todel.Seg = ix; % row of cell array to del
            pipe_idx_todel.Seg = pipe_key.Seg(cell_idx_todel.Seg,1); % name of folder todel
            % 1st descend - all p2s containing seg
            cell_idx_todel.P2  = find( contains( pipe_key.P2(:,2), PipeStruct.Seg_ID ) ); %
            pipe_idx_todel.P2  = pipe_key.P2( cell_idx_todel.P2, 1 );
        else
            error('pipeline block Seg not found -- does this pipeline exist?')
        end
    elseif strcmpi(PartToDelete,'P1') % search P1
        ix = find( strcmp(PipeStruct.P1_ID,pipe_key.P1(:,2)));
        if ~isempty(ix)
            %> Only P2 descend depend on P1
            cell_idx_todel.Warp = [];
            pipe_idx_todel.Warp = [];
            %
            cell_idx_todel.Seg  = [];
            pipe_idx_todel.Seg  = [];
            %
            cell_idx_todel.P1 = ix; % row of cell array to del
            pipe_idx_todel.P1 = pipe_key.P1(cell_idx_todel.P1,1); % name of folder todel
            % 1st descend - all p2s containing p1
            ixl = strfind( PipeStruct.P1_ID, '#P1#' );
            p1wpart = [PipeStruct.P1_ID(1:ixl-2),'-#'];
            p1ppart = [PipeStruct.P1_ID(ixl:end),'-#'];
            cell_idx_todel.P2  = find( contains( pipe_key.P2(:,2), p1wpart ) & contains( pipe_key.P2(:,2), p1ppart ) ); %
            pipe_idx_todel.P2  = pipe_key.P2( cell_idx_todel.P2, 1 );
        else
            error('pipeline block P1 not found -- does this pipeline exist?')
        end
    elseif strcmpi(PartToDelete,'P2') % search P2
        ix = find( strcmp(PipeStruct.PNAME{1},pipe_key.P2(:,1))); % slightly different - match on NAME, then check STRING
        if ~isempty(ix)
            %> P2 has no descendents
            % "authenticate" -> match was found for pname, make sure if corresponds to the same actual sequence
            if ~strcmp(PipeStruct.P2_ID,pipe_key.P2{ix,2})
                error('Pipeline with name %s does not match existing one! What changed?\n\tOld: %s\n\tNew: %s\n',PipeStruct.PNAME{1},pipe_key.P2{ix,2},PipeStruct.P2_ID)
            end
            %
            cell_idx_todel.Warp = [];
            pipe_idx_todel.Warp = [];
            %
            cell_idx_todel.Seg  = [];
            pipe_idx_todel.Seg  = [];
            %
            cell_idx_todel.P1  = [];
            pipe_idx_todel.P1  = [];
            %
            cell_idx_todel.P2 = ix; % row of cell array to del
            pipe_idx_todel.P2 = pipe_key.P2(cell_idx_todel.P2,1); % name of folder todel
        else
            error('pipeline block P2 not found -- does this pipeline exist?')
        end
    elseif strcmpi(PartToDelete,'Analysis')
        % special case --> only remove the analysis results!
        % looks like P2, but we do not delete this segment...
        ix = find( strcmp(PipeStruct.PNAME{1},pipe_key.P2(:,1))); % slightly different - match on NAME, then check STRING
        if ~isempty(ix)
            %> P2 has no descendents
            % "authenticate" -> match was found for pname, make sure if corresponds to the same actual sequence
            if ~strcmp(PipeStruct.P2_ID,pipe_key.P2{ix,2})
                error('Pipeline with name %s does not match existing one! What changed?\n\tOld: %s\n\tNew: %s\n',PipeStruct.PNAME{1},pipe_key.P2{ix,2},PipeStruct.P2_ID)
            end
            %
            cell_idx_todel.Warp = [];
            pipe_idx_todel.Warp = [];
            %
            cell_idx_todel.Seg  = [];
            pipe_idx_todel.Seg  = [];
            %
            cell_idx_todel.P1  = [];
            pipe_idx_todel.P1  = [];
            %
            cell_idx_todel.P2 = []; % row of cell array to del
            pipe_idx_todel.P2 = []; % name of folder todel

            pipe_AnalysisLoc  = pipe_key.P2(ix,1);
        else
            error('pipeline block P2 not found -- does this pipeline exist?')
        end
    else
        error('unrecognized PartToDelete: %s',PartToDelete)
    end

    % save a backup copy in case something goes wrong
    save(fullfile(outpath,'_pipe_manager','_opptmp_pipe_key.mat'),'pipe_key');

    % now update -- by deleting the selected field, then re-saving pipe key
    if ~isempty(cell_idx_todel.Warp)
        pipe_key.Warp( cell_idx_todel.Warp,: ) = [];
    end
    if ~isempty(cell_idx_todel.Seg)
        pipe_key.Seg( cell_idx_todel.Seg,: ) = [];
    end
    if ~isempty(cell_idx_todel.P1)
        pipe_key.P1( cell_idx_todel.P1,: ) = [];
    end
    if ~isempty(cell_idx_todel.P2)
        pipe_key.P2( cell_idx_todel.P2,: ) = [];
    end
    
    save(fullfile(outpath,'_pipe_manager','pipe_key.mat'),'pipe_key');
end

%% Delete from group-level pipeline directory structure

if ~isempty(cell_idx_todel.P2)
    for i=1:numel(pipe_idx_todel.P2)
        unix(sprintf('rm -rf %s', fullfile(outpath,'_group_level','masks',sprintf('pipe_%s',pipe_idx_todel.P2{i}))) );
        unix(sprintf('rm -rf %s', fullfile(outpath,'_group_level','brain_maps',sprintf('pipe_%s',pipe_idx_todel.P2{i}))) );
        unix(sprintf('rm -rf %s', fullfile(outpath,'_group_level','parcellations',sprintf('pipe_%s',pipe_idx_todel.P2{i}))) );
    end
end

%% Delete from subject-level pipeline directory structure

% now store all files to output
for(ns=1:numel(InputStruct))

    % take ith subject input structure
    InputStruct_temp = InputStruct(ns);

    if strcmpi(PartToDelete,'Analysis') % special case
        % just delete the analysis folder!!
        unix(sprintf('rm -rf %s', fullfile( outpath,InputStruct_temp.PREFIX,'func_proc_p2',sprintf('pipe_%s',pipe_idx_todel.P2{1}),ParamStruct_aug.Variable_ID )) );
    else
        % .proceeding from "leaf" (P2) to "trunk" (warp)
        
        % delete all specified p2 pipelines
        if ~isempty(cell_idx_todel.P2)
            for i=1:numel(pipe_idx_todel.P2)
                unix(sprintf('rm -rf %s', fullfile( outpath,InputStruct_temp.PREFIX,'func_proc_p2',sprintf('pipe_%s',pipe_idx_todel.P2{i}) )) );
            end
        end
        % then delete all specified p1 pipelines
        if ~isempty(cell_idx_todel.P1)
            for i=1:numel(pipe_idx_todel.P1)
                unix(sprintf('rm -rf %s', fullfile( outpath,InputStruct_temp.PREFIX,'func_proc_p1',sprintf('subpipe_%03u',pipe_idx_todel.P1{i}) )) ); % _seg gets deleted too
            end
        end
        % then delete all specified segs [SPECIAL CASE - pick through P1/Warp subfolders]
        if ~isempty(cell_idx_todel.Seg)
            for i=1:numel(pipe_idx_todel.Seg)
                % for all extant p1 dirs --> delete seg subpipe if found
                e = dir(sprintf('%s/%s/func_proc_p1/subpipe_*',outpath,InputStruct_temp.PREFIX));
                for j=1:numel(e)
                    unix(sprintf('rm -rf %s/%s/seg_subpipe_%03u',e(j).folder,e(j).name,pipe_idx_todel.Seg{i}));
                end
                e = dir(sprintf('%s/%s/anat_proc/subpipe_*',outpath,InputStruct_temp.PREFIX));
                for j=1:numel(e)
                    unix(sprintf('rm -rf %s/%s/seg_subpipe_%03u',e(j).folder,e(j).name,pipe_idx_todel.Seg{i}));
                end
            end
        end
        % then delete all specified warps
        if ~isempty(cell_idx_todel.Warp)
            for i=1:numel(pipe_idx_todel.Warp)
                unix(sprintf('rm -rf %s', fullfile( outpath,InputStruct_temp.PREFIX,'anat_proc',sprintf('subpipe_%03u',pipe_idx_todel.Warp{i}) )) ); % _seg gets deleted too
            end
        end
    end
end

% if it made it here successfully, delete the backup copy
delete(fullfile(outpath,'_pipe_manager','_opptmp_pipe_key.mat'));

disp('Step-0 (DELETION) Complete!');
