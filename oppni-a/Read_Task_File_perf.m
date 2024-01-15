function TaskStruct = Read_Task_File_perf( taskfile ) % ,'ANALYSIS','CONTRAST','VOXRES'
% *PERF (nochange)
TaskStruct=[];

if( ~isempty(taskfile) )

    fid   = fopen(taskfile);
    if fid==-1
        error('cannot open task file:\n\t%s\n',taskfile);
    end
    tline = fgetl(fid);
    if ~ischar(tline)
        error('task file is empty:\n\t%s\n',taskfile);
    end

    n_name     = 0;
    n_onsets   = 0;
    n_duration = 0;
    namelist   = {}; %% to check for duplicates
    
    %% read input file
    while ischar(tline) 

        % [] bracketing fields
        istart = strfind(tline,'[')+1;
        iend   = strfind(tline,']')-1;   

        if( ~isempty(istart) && ~isempty(iend) )
            
            %% contrast based analysis fields...

            % check for condition-specific options
            if contains(upper(tline),'STIM_LABEL')
                n_name = n_name+1;
                if( numel(regexp( tline(istart:iend), '[^a-zA-Z0-9_]' ))>0 )
                    error('split file info: %s condition names can only include alphanumeric values and underscores.',taskfile);
                end
                TaskStruct.cond(n_name).label = tline(istart:iend);
                namelist = [namelist, {TaskStruct.cond(n_name).label}];
            end
            % 
            if contains(upper(tline),'STIM_ONSET')% || ~isempty(strfind(upper(tline),'ONSETS')))

                n_onsets = n_onsets+1;
                % stores numeric onsets -- if non, this is an empty cell
                e = regexp( tline(istart:iend), ',','split' );
                onsetsArray{n_onsets} = [];
                for e_counter_temp = 1:length(e)
                    onsetsArray{n_onsets} = [onsetsArray{n_onsets} str2num(e{e_counter_temp})];
                end
            end      
            % 
            if contains(upper(tline),'STIM_DURATION')% || ~isempty(strfind(upper(tline),'DURATIONS')))

                n_duration = n_duration+1;
                % stores numeric durations -- if non, this is an empty cell
                e = regexp( tline(istart:iend), ',','split' );                
                durationArray{n_duration} = [];
                for e_counter_temp = 1:length(e)
                     durationArray{n_duration}= [ durationArray{n_duration} str2num(e{e_counter_temp})];
                end
            end
        end

        tline = fgetl(fid);
        if isempty(tline)
            tline = fgetl(fid);
        end
    end
    fclose(fid); 

    if( length(unique(namelist)) < length(namelist) )
        error('duplicate conditions found in task design of %s',taskfile);
    end

    % contrast-based quality checks
    if( (n_name == n_onsets) && (n_name == n_duration) && (n_duration == n_onsets) )

        for(n=1:n_name)

            if( length(onsetsArray{n}) == length(durationArray{n}) )
                %
                TaskStruct.cond(n).onset    = onsetsArray{n};                
                TaskStruct.cond(n).duration = durationArray{n};
            else
                disp(taskfile);
                error('number of STIM_ONSET and STIM_DURATION values must be the same for every condition');
            end
        end
    else
        error('every task condition in split-info. requires a STIM_LABEL, list of STIM_ONSET and STIM_DURATION');
    end    
end
