function [InputStruct_return] = interpret_contrast_list( InputStruct, analysis_model, contrast_list_str)
%
%  Takes raw task data (in InputStruct.task_unformat{run}) and analysis model / contrast list
%  and creates new task struct (in InputStruct.task_contrast{run})
%  used to execute the analysis
%
%  this includes:
%  .reorganizing the conditions to discard those of no-interest
%  .checking that all specified conditions are present
%  .reordering the remaining conditions alphabetically
%  .standardizing stimulus units in msec
%  .adjusting onsets for dropped volumes
%  .creating the "simple" design matrix and convolving w HRF
%
%  InputStruct_return.task_info{run}.contrast(n).label = (contrast string)
%  InputStruct_return.task_info{run}.contrast(n).design_mat = (design matrix)
%  InputStruct_return.task_info{run}.contrast(n).c1(/c2) = (columns to contrast, regression analysis)
%  InputStruct_return.task_info{run}.contrast(n).o1(/o2) = (onset indices, event-related analysis)
%  InputStruct_return.task_info{run}.contrast(n).b1(/b2) = (binary task-on indices, 2class analysis)
%

for nr = 1:InputStruct.N_func
    
    if( analysis_model.uses_taskfile>0 )

        if ~isfield(InputStruct,'task_unformat')
            error(['For ',InputStruct.PREFIX,': unformatted task file of run ',num2str(nr),' doesnt exist']);
        else
            % extract unformatted taskfile to do some modifications
            task_info = InputStruct.task_unformat{nr};
        end

        if( isempty(contrast_list_str) || strcmpi(contrast_list_str,'NONE') )
            error('task-based analysis, need to specify a contrast');
        elseif( ~isfield(task_info,'cond') )
            error('task-based analysis, task conditions do not seem to be stored to file'); 
        end
            
        %% checking all conditions of interest exist - discard no-interest

        % conditions being analyzed in pipelines
        conditions_of_interest = unique( regexp(contrast_list_str,'[,+-]','split') );

        % (1) drop conditions not present in contrast(s) of interest
        clear no_interest;
        keepfields={};
        for i=1:numel(task_info.cond)
            if( sum(strcmpi(task_info.cond(i).label,conditions_of_interest))==0 )
                no_interest(i) = true; 
            else
                no_interest(i) = false;
                keepfields = [keepfields, {task_info.cond(i).label}];
            end
        end        
        task_info.cond( no_interest ) = []; %% delete fields of no-interest
        clear no_interest;
        if numel(keepfields) ~= numel(task_info.cond)
           error('number of kept fields does not match number of conditions??');
        end

        % (2) check if all required conditions are present in file, and if "kept" conditions all have onsets/durations/etc. 
        for i=1:numel(conditions_of_interest)
            if sum(strcmpi(conditions_of_interest{i},keepfields))==0
                error('Condition "%s" in contrast not defined for file: %s',conditions_of_interest{i},InputStruct.frun(nr).TASK_filename);
            end
        end
        for i=1:numel(task_info.cond)
            if  isempty( task_info.cond(i).onset ) || isempty( task_info.cond(i).duration )
                error('Condition "%s" in contrast has no onsets and/or durations in file: %s',task_info.cond(i).label,InputStruct.frun(nr).TASK_filename);                
            end
        end

        % (3) reorder conditions alphabetically, for ease of output checking
        [keepfields, ix] = sort( keepfields );
        task_info_temp = task_info;
        for i=1:numel(task_info.cond)
            task_info.cond(i) = task_info_temp.cond(ix(i));
        end
        clear task_info_temp ix;

        %% now interpret and reformat input conditions

        % number of "real" conditions
        num_kept_cond = numel(task_info.cond);

        % convert sec to msec / peel off "bad" time settings?
        for i = 1:length(task_info.cond)
            
            task_info.cond(i).onset = task_info.cond(i).onset*1000;
            task_info.cond(i).duration = task_info.cond(i).duration*1000;
            
            %-- removal of *extremely* offset condition evens
            
            % remove events that terminate earlier than 8 sec. before adjusted start time (as HRF amp dropped back to 50% amp after 8 sec.)
            ind_temp = (task_info.cond(i).onset+task_info.cond(i).duration) < (InputStruct.frun(nr).DROP_first*InputStruct.TR_MSEC - 8000);
            task_info.cond(i).onset(ind_temp) = [];
            task_info.cond(i).duration(ind_temp) = [];
            % remove events that start later than 3 sec. before adjusted end time (as HRF amp reached 50% amp after 3 sec.)
            ind_temp = task_info.cond(i).onset > ( InputStruct.frun(nr).Nt_adj*InputStruct.TR_MSEC - InputStruct.frun(nr).DROP_last*InputStruct.TR_MSEC - 3000 );
            task_info.cond(i).onset(ind_temp) = [];
            task_info.cond(i).duration(ind_temp) = [];
            % catch cases where this this step strips off all valid onsets
            if( isempty(task_info.cond(i).onset) ) 
                InputStruct.frun(nr).FUNC_filename,
                InputStruct.frun(nr).TASK_filename,
                error(['condition "',task_info.cond(i).label,'" has exclusively blocks that start before the start/after the end of the fMRI run! check your onsets!!']); 
            end
        end

        %--------------- start creating design mat ----------------%

        % .we create design matrix, with all ORIGINAL conditions included 
        % .start off with sub-sampled matrix and interpolate to correct TR intervals

        Nsubs  = max([1 round(InputStruct.TR_MSEC/100)]); % subsampled design matrix (#samples per TR) w catch in case TR<100ms
        design_raw = zeros( InputStruct.frun(nr).Nt_raw*Nsubs, num_kept_cond);     % initialize design matrix
        %fill-in stim blocks
        for cond_counter = 1:num_kept_cond %% for each real (non-combined) condition
            for onset_counter = 1:length(task_info.cond(cond_counter).onset)
                st = round(task_info.cond(cond_counter).onset(onset_counter)./(InputStruct.TR_MSEC/Nsubs)) + 1;
                ed = st + round(task_info.cond(cond_counter).duration(onset_counter)./(InputStruct.TR_MSEC/Nsubs));
                if ed>size(design_raw,1)
                    ed = size(design_raw,1);
                end
                design_raw(st:ed,cond_counter) = 1;
            end
        end
        % convolve with HRF - design_to_hrf function requires that we convert TR to sec.
        design_mat = design_to_hrf( design_raw, (InputStruct.TR_MSEC/Nsubs)/1000, [5.0 15.0] );
        % now, adjusting for "dropped" initial/later scans:
        design_mat = design_mat( (InputStruct.frun(nr).DROP_first * Nsubs + 1):(end - InputStruct.frun(nr).DROP_last * Nsubs), : );
        design_raw = design_raw( (InputStruct.frun(nr).DROP_first * Nsubs + 1):(end - InputStruct.frun(nr).DROP_last * Nsubs), : );
        if size(design_mat,1) ~= InputStruct.frun(nr).Nt_adj*Nsubs || size(design_raw,1) ~= InputStruct.frun(nr).Nt_adj*Nsubs
            error('excess overhang in design mat -- pls inspect')
            % trimming back in case of conv-overhang
            % design_mat = design_mat(1:InputStruct.frun(nr).Nt_adj*Nsubs,:);
        end
        % now, subsample to mid-TR to get HRF at "real" fMRI sampling rate
        task_info.design_mat    = design_mat( round(Nsubs/2): Nsubs : end, : );
        task_info.design_noconv = design_raw( round(Nsubs/2): Nsubs : end, : );
        % renorm s.t. inner product = average on-condition activity
        task_info.design_mat = bsxfun(@rdivide, task_info.design_mat, sum(task_info.design_mat,1)+eps);

        %--------------- done creating design mat ----------------%
        
        % now for onsets, adjust start time based on DROP_first --- this will give us negative values for v. early onsets!
        for i = 1:length(task_info.cond)
            task_info.cond(i).onset = task_info.cond(i).onset - InputStruct.frun(nr).DROP_first*InputStruct.TR_MSEC;
        end

        %% now creating structures for specified analysis model

        ca = regexp(contrast_list_str,',','split');
        for i=1:numel(ca)
            %--
            task_info.contrast(i).c1 = [];
            task_info.contrast(i).o1 = [];
            task_info.contrast(i).b1 = [];
            task_info.contrast(i).c2 = [];
            task_info.contrast(i).o2 = [];
            task_info.contrast(i).b2 = [];
            %--
            task_info.contrast(i).label = ca{i}; % store label
            cb=regexp(ca{i},'-','split');
            %--
            cc=regexp(cb{1},'+','split');
            task_info.contrast(i).c1=[];
            for j=1:numel(cc)
                ic = find( strcmpi(cc{j}, keepfields) );
                tmpo = ceil( task_info.cond(ic).onset ./ InputStruct.TR_MSEC );
                tmpb = task_info.design_mat(:,ic);
                tmpb = (tmpb - min(tmpb))./(max(tmpb)-min(tmpb));
                task_info.contrast(i).c1 = [task_info.contrast(i).c1; ic(:)];
                task_info.contrast(i).o1 = [task_info.contrast(i).o1 tmpo(:)];
                task_info.contrast(i).b1 = [task_info.contrast(i).b1; find(tmpb(:)>0.90) ];
            end
            %--
            if numel(cb)==1
                task_info.contrast(i).c2 = [];
                task_info.contrast(i).o2 = [];
                task_info.contrast(i).b2 = [];
            elseif numel(cb)==2
                %--
                cc=regexp(cb{2},'+','split');
                task_info.contrast(i).c2=[];
                for j=1:numel(cc)
                    ic = find( strcmpi(cc{j}, keepfields) );
                    tmpo = ceil( task_info.cond(ic).onset ./ InputStruct.TR_MSEC );
                    tmpb = task_info.design_mat(:,ic);
                    tmpb = (tmpb - min(tmpb))./(max(tmpb)-min(tmpb));
                    task_info.contrast(i).c2 = [task_info.contrast(i).c2; ic(:)];
                    task_info.contrast(i).o2 = [task_info.contrast(i).o2 tmpo(:)];
                    task_info.contrast(i).b2 = [task_info.contrast(i).b2; find(tmpb(:)>0.90) ];
                end
                %--
            else
                error('too many subtractive differences?')
            end
        end
    else
        % case with no contrasts specified
        task_info.contrast = [];
    end
    
    % now storing parsed task_info into input structure
    InputStruct.task_info{nr} = task_info;
end

% copy to output
InputStruct_return = InputStruct;
