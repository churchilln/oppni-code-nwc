function [InputStruct_return] = interpret_seed_list( InputStruct, analysis_model, contrast_list_str)
% 
%  Takes raw seed data (in InputStruct.seed_unformat{run}) and analysis model / contrast list
%  and creates new task struct (in InputStruct.seed_info)
%  used to execute the analysis --> note that there is only one entry, applied to all runs, unlike task_info, which has one cell per run 
%
%  this includes:
%  .reorganizing the conditions to discard those of no-interest
%  .checking that all specified conditions are present
%  .reordering the remaining conditions alphabetically
%  .also including files ordered according to the contrasts
%
%  InputStruct_return.seed_info.contrast(n).label = (contrast string)
%  InputStruct_return.seed_info.contrast(n).location = (location of seedfile)
%  InputStruct_return.seed_info.contrast(n).type = (type of seedfile)
%

if analysis_model.uses_roifile>0

    if ~isfield(InputStruct,'seed_unformat')
        error(['For ',InputStruct.PREFIX,': unformatted seed_info file doesnt exist']);
    else
        % extract unformatted seedfile to do some modifications
        seed_info = InputStruct.seed_unformat;
    end

    if isempty(contrast_list_str) || strcmpi(contrast_list_str,'NONE')
        error('seed-based analysis, need to specify a contrast');
    elseif( ~isfield(seed_info,'seed') )
        error('seed-based analysis, seeds do not seem to be stored to file'); 
    end
    
    %% checking all conditions of interest exist - discard no-interest

    if contains(contrast_list_str,'-') || contains(contrast_list_str,'+')
       error('seed analysis does not support additive/subtractive operations (+,-) on seeds.'); 
    end

    % conditions being analyzed in pipelines
    conditions_of_interest = unique( regexp(contrast_list_str,',','split') );

    % (1) drop conditions not present in contrast(s) of interest
    clear no_interest;
    keepfields={};
    for i=1:numel(seed_info.seed)
        if( sum(strcmpi(seed_info.seed(i).label,conditions_of_interest))==0 )
            no_interest(i) = true; 
        else
            no_interest(i) = false;
            keepfields = [keepfields, {seed_info.seed(i).label}];
        end
    end        
    seed_info.seed( no_interest ) = []; %% delete fields of no-interest
    clear no_interest;

    % (2) check if all required conditions are present in file, and if "kept" conditions all have onsets/durations/etc. 
    for i=1:numel(conditions_of_interest)
        if sum(strcmpi(conditions_of_interest{i},keepfields))==0
            error('Condition "%s" in contrast not defined for file: %s',conditions_of_interest{i},InputStruct.SEED_filename);
        end
    end
    for i=1:numel(seed_info.seed)
        if  isempty( seed_info.seed(i).location ) || isempty( seed_info.seed(i).type )
            error('Condition "%s" in contrast has no location and/or type in file: %s',seed_info.seed(i).label,InputStruct.SEED_filename);                
        end
    end
    
    % (3) reorder conditions alphabetically, for ease of output checking
    [keepfields, ix] = sort( keepfields );
    seed_info_temp = seed_info;
    for i=1:numel(seed_info.seed)
        seed_info.seed(i) = seed_info_temp.seed(ix(i));
    end
    clear seed_info_temp ix;

    %% now interpret and reformat input conditions

    % number of "real" conditions
    num_kept_cond = numel(seed_info.seed);

    ca = regexp(contrast_list_str,',','split');
    for i=1:numel(ca)
        % port over entire seed struct, with indexing...
        ix = find( strcmpi(ca{i}, keepfields) );
        %-- 
        seed_info.contrast(i).label    = seed_info.seed(ix).label;
        seed_info.contrast(i).location = seed_info.seed(ix).location;
        seed_info.contrast(i).type     = seed_info.seed(ix).type;
        %--
        %seed_info.contrast(i) = seed_info.seed(ix);
        seed_info.contrast(i).c1 = ix;

        if numel(seed_info.contrast(i).c1) ~= 1
            error('something went wrong with seed contrast specification');
        end
        
        if strcmpi(seed_info.contrast(i).type,'multi')
            ismulti(i,1) = 1;
        elseif strcmpi(seed_info.contrast(i).type,'single')
            ismulti(i,1) = 0;
        else
            error('?!');
        end
        if sum(ismulti==1)>0
            if sum(ismulti==0)>0
                error('seed method cannot support mixed types (single+multi)!')
            elseif numel(seed_info.contrast)>1
                error('seed method cannot support more than one atlas (multi)!')
            end
        end
    end
else
    % case with no contrasts specified
    seed_info.contrast = [];  %% zero-field if no contrasts specified
end

% now storing parsed task_info into input structure
InputStruct.seed_info = seed_info;

% copy to output
InputStruct_return = InputStruct;
