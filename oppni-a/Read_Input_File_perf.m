function InputStruct = Read_Input_File_perf(inputfile)
% *PERF
% InputStruct(subj).(field) = string or numeric value
% InputStruct(subj).arun(func.run).(field) = cell array of strings or numeric values 
% InputStruct(subj).prun(func.run).(field) = cell array of strings or numeric values 
%

InputStruct=[];
if ~contains(inputfile,'=')
    fid   = fopen(inputfile);
    if fid==-1
        error('cannot open input file:\n\t%s\n',inputfile);
    end
    tline = fgetl(fid);
    if ~ischar(tline)
        error('input file is empty:\n\t%s\n',inputfile);
    end
else
    tline = inputfile; %% just directly assign the string
end

%% read each line of input file, parse into relevant structure
ns=0;
while ischar(tline) 

    ns = ns+1;
    tline     = regexprep(tline, '\t', ' ');             %% convert tabs to spaces
    ispaces   = [strfind( tline, ' ' )-1 length(tline)]; %% index all spaces (-1), and eol

    %% prefix...
    fieldname = 'PREFIX';
    if contains( upper(tline), strcat(fieldname,'=') )
        isrt = strfind( upper(tline),[fieldname,'='] ) + (numel(fieldname)+1);
        iend = ispaces(ispaces>=isrt);
        InputStruct(ns).PREFIX = tline(isrt:iend(1)); % take first ending
    else
        error('input file line does not contain mandatory %s field\n\t%s\n',fieldname,tline);
    end

    fieldname = 'LAB_METHOD'; % [PASL/CASL/PCASL]-[2D/3D]
    if contains( upper(tline), strcat(fieldname,'=') )
        isrt = strfind( upper(tline),[fieldname,'='] ) + (numel(fieldname)+1);
        iend = ispaces(ispaces>=isrt);
        InputStruct(ns).LAB_METHOD = tline(isrt:iend(1)); % take first ending
    else
        error('input file line does not contain mandatory %s field\n\t%s\n',fieldname,tline);
    end

    fieldname = 'LAB_ORDER'; % TC or CT
    if contains( upper(tline), strcat(fieldname,'=') )
        isrt = strfind( upper(tline),[fieldname,'='] ) + (numel(fieldname)+1);
        iend = ispaces(ispaces>=isrt);
        InputStruct(ns).LAB_ORDER = tline(isrt:iend(1)); % take first ending
    else
        error('input file line does not contain mandatory %s field\n\t%s\n',fieldname,tline);
    end

    fieldname = 'TR_MSEC';
    if contains( upper(tline), strcat(fieldname,'=') )
        isrt = strfind( upper(tline),[fieldname,'='] ) + (numel(fieldname)+1);
        iend = ispaces(ispaces>=isrt);
        InputStruct(ns).TR_MSEC = str2num(tline(isrt:iend(1))); % take first ending
    else
        error('input file line does not contain mandatory %s field\n\t%s\n',fieldname,tline);
    end

    fieldname = 'TE_MSEC';
    if contains( upper(tline), strcat(fieldname,'=') )
        isrt = strfind( upper(tline),[fieldname,'='] ) + (numel(fieldname)+1);
        iend = ispaces(ispaces>=isrt);
        InputStruct(ns).TE_MSEC = str2num(tline(isrt:iend(1))); % take first ending
    else
        error('input file line does not contain mandatory %s field\n\t%s\n',fieldname,tline);
    end

    fieldname = 'LDS_MSEC'; % labl / delay / slice times ... in millisec
    % Labl=TI1 for PASL / label-time for CASL; Delay=TI/TI2 for PASL / PLD for CASL
    % Slicetimes=0 for 3d acq / -1 for auto-estim
    if contains( upper(tline), strcat(fieldname,'=') )
        isrt = strfind( upper(tline),[fieldname,'='] ) + (numel(fieldname)+1);
        iend = ispaces(ispaces>=isrt);
        InputStruct(ns).LDS_MSEC = str2num(tline(isrt:iend(1))); % take first ending
    else
        error('input file line does not contain mandatory %s field\n\t%s\n',fieldname,tline);
    end

    %% perfusion data and derivatives...

    fieldname = 'PERF';
    if contains( upper(tline), strcat(fieldname,'=') )
        isrt = strfind( upper(tline),[fieldname,'='] ) + (numel(fieldname)+1);
        iend = ispaces(ispaces>=isrt);
        filestring_temp = tline(isrt:iend(1)); % take first ending
        if  ~contains(filestring_temp,',')
             filestring_temp = {filestring_temp};
        else filestring_temp = regexp(filestring_temp,',','split'); 
        end
        InputStruct(ns).N_perf = numel( filestring_temp ); % number of functional runs
        for nr=1:InputStruct(ns).N_perf
            InputStruct(ns).prun(nr).PERF_filename = filestring_temp{nr};
        end
    else
        error('input file line does not contain mandatory %s field\n\t%s\n',fieldname,tline);
    end

    fieldname = 'TASK';
    if contains( upper(tline), strcat(fieldname,'=') )
        isrt = strfind( upper(tline),[fieldname,'='] ) + (numel(fieldname)+1);
        iend = ispaces(ispaces>=isrt);
        filestring_temp = tline(isrt:iend(1)); % take first ending
        if  ~contains(filestring_temp,',')
             filestring_temp = {filestring_temp};
        else filestring_temp = regexp(filestring_temp,',','split'); 
        end
        if numel(filestring_temp) ~= InputStruct(ns).N_perf
            error('number of TASK args does not match number of perf runs. Check line:\n\t%s\n',tline);
        end
        for nr=1:InputStruct(ns).N_perf
            InputStruct(ns).prun(nr).TASK_filename = filestring_temp{nr};
        end
    else
        warning('input file line does not contain optional %s field - no task-based analysis\n\t%s\n',fieldname,tline);
        InputStruct(ns).prun(nr).TASK_filename = [];
    end

    fieldname = 'M0REF'; % comma-separated list; either string list of filenames / or / numeric indices
    if contains( upper(tline), strcat(fieldname,'=') )
        isrt = strfind( upper(tline),[fieldname,'='] ) + (numel(fieldname)+1);
        iend = ispaces(ispaces>=isrt);
        filestring_temp = tline(isrt:iend(1)); % take first ending
        if  ~contains(filestring_temp,',')
             filestring_temp = {filestring_temp};
        else filestring_temp = regexp(filestring_temp,',','split'); 
        end

        if contains(filestring_temp,'NONE')
            warning('input file line does not contain optional %s field - any absolute CBF estimates will be based on a "control" scan\n\t%s\n',fieldname,tline);
            InputStruct(ns).N_m0ref=InputStruct(ns).N_perf; % we will be pulling one "control" scan from each perf run
            InputStruct(ns).m0loc = 'none';
        else
            InputStruct(ns).N_m0ref = numel( filestring_temp ); % number of functional runs
            numre=0;
            for nr=1:InputStruct(ns).N_m0ref
                if contains(filestring_temp{nr},'.nii')
                    InputStruct(ns).mrun(nr).M0REF_filename = filestring_temp{nr};
                elseif numel(filestring_temp{nr})<=3 % format as numeric
                    InputStruct(ns).mrun(nr).M0REF_filename = str2num(filestring_temp{nr});
                    numre=numre+1;
                else
                    error('m0ref format unclear??')
                end
                if numre>0 
                    if InputStruct(ns).N_m0ref ~= InputStruct(ns).N_perf
                        error('if calibration scan(s) part of the perf-weighted run(s), need an index number for each!');
                    end
                    InputStruct(ns).m0loc='infile';
                else
                    InputStruct(ns).m0loc='separate';
                end
            end
        end
    else
        warning('input file line does not contain optional %s field - any absolute CBF estimates will be based on a "control" scan\n\t%s\n',fieldname,tline);
        InputStruct(ns).N_m0ref=InputStruct(ns).N_perf; % we will be pulling one "control" scan from each perf run
        InputStruct(ns).m0loc = 'none';
    end

    fieldname = 'PWDROP'; % assumed to be a ginle global value
    if contains( upper(tline), strcat(fieldname,'=') )
        isrt = strfind( upper(tline),[fieldname,'='] ) + (numel(fieldname)+1);
        iend = ispaces(ispaces>=isrt);
        InputStruct(ns).PWDROP = tline(isrt:iend(1)); % take first ending
    else
        InputStruct(ns).PWDROP = 'NONE';
    end

    %% anatomical data and derivatives...

    fieldname = 'ANAT';
    if contains( upper(tline), strcat(fieldname,'=') )
        isrt = strfind( upper(tline),[fieldname,'='] ) + (numel(fieldname)+1);
        iend = ispaces(ispaces>=isrt);
        filestring_temp = tline(isrt:iend(1)); % take first ending
        if  ~contains(filestring_temp,',')
             filestring_temp = {filestring_temp};
        else filestring_temp = regexp(filestring_temp,',','split'); 
        end
        InputStruct(ns).N_anat = numel( filestring_temp ); % number of perfusion runs
        for nr=1:InputStruct(ns).N_anat
            InputStruct(ns).arun(nr).ANAT_filename   = filestring_temp{nr};
        end
    else
        error('input file line does not contain mandatory %s field\n\t%s\n',fieldname,tline);
    end

    fieldname = 'ZCLIP';
    if contains( upper(tline), strcat(fieldname,'=') )
        isrt = strfind( upper(tline),[fieldname,'='] ) + (numel(fieldname)+1);
        iend = ispaces(ispaces>=isrt);
        filestring_temp = tline(isrt:iend(1)); % take first ending
        if  ~contains(filestring_temp,',')
             filestring_temp = {filestring_temp};
        else filestring_temp = regexp(filestring_temp,',','split'); 
        end
        for nr=1:InputStruct(ns).N_anat
            if strcmpi(filestring_temp{nr},'AUTO')
                InputStruct(ns).arun(nr).ZCLIP_thr  = 'AUTO';
            else
                InputStruct(ns).arun(nr).ZCLIP_thr  = str2num(filestring_temp{nr});
            end
        end
    else
        warning('input file line does not contain optional %s field - no clipping of T1 scan\n\t%s\n',fieldname,tline);
        for nr=1:InputStruct(ns).N_anat
            InputStruct(ns).arun(nr).ZCLIP_thr  = 'NONE';
        end
    end

    if( isempty(strfind(inputfile,'=')) )    
        tline = fgetl(fid);
        if isempty(tline)
            tline = fgetl(fid);
        end
    else
        tline=[]; %% empty the string to terminate
    end
end

if( isempty(strfind(inputfile,'=')) )    
    fclose(fid); 
end
