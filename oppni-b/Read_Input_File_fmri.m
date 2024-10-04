function InputStruct = Read_Input_File_fmri(inputfile)
%
% InputStruct(subj).(field) = string or numeric value
% InputStruct(subj).arun(func.run).(field) = cell array of strings or numeric values 
% InputStruct(subj).frun(func.run).(field) = cell array of strings or numeric values 
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

    fieldname = 'TPATTERN';
    if contains( upper(tline), strcat(fieldname,'=') )
        isrt = strfind( upper(tline),[fieldname,'='] ) + (numel(fieldname)+1);
        iend = ispaces(ispaces>=isrt);
        InputStruct(ns).TPATTERN = tline(isrt:iend(1)); % take first ending
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

    fieldname = 'SEED';
    if contains( upper(tline), strcat(fieldname,'=') )
        isrt = strfind( upper(tline),[fieldname,'='] ) + (numel(fieldname)+1);
        iend = ispaces(ispaces>=isrt);
        InputStruct(ns).SEED_filename = tline(isrt:iend(1)); % take first ending
    else
        warning('input file line does not contain optional %s field - no seed-based analysis\n\t%s\n',fieldname,tline);
        InputStruct(ns).SEED_filename = [];
    end

    %% functional data and derivatives...

    fieldname = 'FUNC';
    if contains( upper(tline), strcat(fieldname,'=') )
        isrt = strfind( upper(tline),[fieldname,'='] ) + (numel(fieldname)+1);
        iend = ispaces(ispaces>=isrt);
        filestring_temp = tline(isrt:iend(1)); % take first ending
        if  ~contains(filestring_temp,',')
             filestring_temp = {filestring_temp};
        else filestring_temp = regexp(filestring_temp,',','split'); 
        end
        InputStruct(ns).N_func = numel( filestring_temp ); % number of functional runs
        for nr=1:InputStruct(ns).N_func
            InputStruct(ns).frun(nr).FUNC_filename = filestring_temp{nr};
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
        if numel(filestring_temp) ~= InputStruct(ns).N_func
            error('number of TASK args does not match number of functional runs. Check line:\n\t%s\n',tline);
        end
        for nr=1:InputStruct(ns).N_func
            InputStruct(ns).frun(nr).TASK_filename = filestring_temp{nr};
        end
    else
        warning('input file line does not contain optional %s field - no task-based analysis\n\t%s\n',fieldname,tline);
        InputStruct(ns).frun(nr).TASK_filename = [];
    end

    fieldname = 'PHYSIO';
    if contains( upper(tline), strcat(fieldname,'=') )
        isrt = strfind( upper(tline),[fieldname,'='] ) + (numel(fieldname)+1);
        iend = ispaces(ispaces>=isrt);
        filestring_temp = tline(isrt:iend(1)); % take first ending
        if  ~contains(filestring_temp,',')
             filestring_temp = {filestring_temp};
        else filestring_temp = regexp(filestring_temp,',','split'); 
        end
        if numel(filestring_temp) ~= InputStruct(ns).N_func
            error('number of PHYSIO args does not match number of functional runs. Check line:\n\t%s\n',tline);
        end
        for nr=1:InputStruct(ns).N_func
            InputStruct(ns).frun(nr).PHYSIO_filename = filestring_temp{nr};
        end
    else
        warning('input file line does not contain optional %s field - no model-based physio correction\n\t%s\n',fieldname,tline);
        for nr=1:InputStruct(ns).N_func
            InputStruct(ns).frun(nr).PHYSIO_filename = [];
        end
    end

    fieldname = 'PHYSAMP_MSEC';
    if contains( upper(tline), strcat(fieldname,'=') )
        isrt = strfind( upper(tline),[fieldname,'='] ) + (numel(fieldname)+1);
        iend = ispaces(ispaces>=isrt);
        InputStruct(ns).PHYSAMP_PR_MSEC = str2num(tline(isrt:iend(1))); % take first ending
    else
        warning('input file line does not contain optional %s field - some physio correction may not work\n\t%s\n',fieldname,tline);
        InputStruct(ns).PHYSAMP_PR_MSEC = [];
    end

    fieldname = 'DROP';
    if contains( upper(tline), strcat(fieldname,'=') )
        isrt = strfind( upper(tline),[fieldname,'='] ) + (numel(fieldname)+1);
        iend = ispaces(ispaces>=isrt);
        filestring_temp = tline(isrt:iend(1)); % take first ending
        ibrak1 = strfind(filestring_temp,'[');
        ibrak2 = strfind(filestring_temp,']');
        if( numel(ibrak1)~=numel(ibrak2) )
            error('mistake in DROP call. Check line:\n\t%s\n',tline);
        elseif numel(ibrak1) ~= InputStruct(ns).N_func
            error('number of DROP args does not match number of functional runs. Check line:\n\t%s\n',tline);
        else
            for(j=1:numel(ibrak1))
                filestring_temp2{j} = filestring_temp( ibrak1(j)+1 : ibrak2(j)-1 );
                icomm = strfind(filestring_temp2{j},',');
                if( isempty(icomm) || numel(icomm)>1 ) error('mistake in DROP call. Check line:\n\t%s\n',tline); end
            end
            filestring_temp = filestring_temp2; clear filestring_temp2 icomm;
        end                 
        for nr=1:InputStruct(ns).N_func
            dropcell = regexp(filestring_temp{nr},',','split');
            InputStruct(ns).frun(nr).DROP_first      = str2num(dropcell{1});
            InputStruct(ns).frun(nr).DROP_last       = str2num(dropcell{2});
        end
    else
        warning('input file line does not contain mandatory %s field - no volumes will be dropped\n\t%s\n',fieldname,tline);
        for nr=1:InputStruct(ns).N_func
            InputStruct(ns).frun(nr).DROP_first      = 0;
            InputStruct(ns).frun(nr).DROP_last       = 0;
        end
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
        InputStruct(ns).N_anat = numel( filestring_temp ); % number of functional runs
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
