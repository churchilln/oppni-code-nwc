function InputStruct = Read_Input_File_diff(inputfile)
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

    fieldname = 'REV_MODE'; % format of reverse phase encoding (if present)
    if contains( upper(tline), strcat(fieldname,'=') )
        isrt = strfind( upper(tline),[fieldname,'='] ) + (numel(fieldname)+1);
        iend = ispaces(ispaces>=isrt);
        InputStruct(ns).REV_MODE = tline(isrt:iend(1)); % take first ending
    else
        error('input file line does not contain mandatory %s field\n\t%s\n',fieldname,tline);
    end

    fieldname = 'PE_FWD'; % direction of forward phase encoding
    if contains( upper(tline), strcat(fieldname,'=') )
        isrt = strfind( upper(tline),[fieldname,'='] ) + (numel(fieldname)+1);
        iend = ispaces(ispaces>=isrt);
        InputStruct(ns).PE_FWD = tline(isrt:iend(1)); % take first ending
    else
        error('input file line does not contain mandatory %s field\n\t%s\n',fieldname,tline);
    end

    fieldname = 'PE_REV'; % direction of reverse phase encoding *
    if contains( upper(tline), strcat(fieldname,'=') )
        isrt = strfind( upper(tline),[fieldname,'='] ) + (numel(fieldname)+1);
        iend = ispaces(ispaces>=isrt);
        InputStruct(ns).PE_REV = tline(isrt:iend(1)); % take first ending
    else
        warning('input file line does not contain optional %s field - no reverse phase encoding supported\n\t%s\n',fieldname,tline);
        InputStruct(ns).PE_REV = [];
    end

    fieldname = 'TRO_MSEC'; % total readout time in milliseconds
    if contains( upper(tline), strcat(fieldname,'=') )
        isrt = strfind( upper(tline),[fieldname,'='] ) + (numel(fieldname)+1);
        iend = ispaces(ispaces>=isrt);
        InputStruct(ns).TRO_MSEC = str2num(tline(isrt:iend(1))); % take first ending
    else
        error('input file line does not contain mandatory %s field\n\t%s\n',fieldname,tline);
    end

% %     fieldname = 'B0REF_IX'; % index of [forward] or [forward,reverse] pe volume used as reference
% %     if contains( upper(tline), strcat(fieldname,'=') )
% %         isrt = strfind( upper(tline),[fieldname,'='] ) + (numel(fieldname)+1);
% %         iend = ispaces(ispaces>=isrt);
% %         InputStruct(ns).PE_FWD = tline(isrt:iend(1)); % take first ending
% %     else
% %         error('input file line does not contain mandatory %s field\n\t%s\n',fieldname,tline);
% %     end
% %     % "F" or "L" means take first or last of b=0 scans in catted set
% %     % "F+(N)" or "L-(N)" means take the first/last +- (N) scans

    %% diffusion data and derivatives...

    fieldname = 'DIFF_FWD';
    if contains( upper(tline), strcat(fieldname,'=') )
        isrt = strfind( upper(tline),[fieldname,'='] ) + (numel(fieldname)+1);
        iend = ispaces(ispaces>=isrt);
        filestring_temp = tline(isrt:iend(1)); % take first ending
        if  ~contains(filestring_temp,',')
             filestring_temp = {filestring_temp};
        else filestring_temp = regexp(filestring_temp,',','split'); 
        end
        InputStruct(ns).N_diff_fwd = numel( filestring_temp ); % number of functional runs
        for nr=1:InputStruct(ns).N_diff_fwd
            InputStruct(ns).dfrun(nr).DIFF_FWD_filename = filestring_temp{nr};
        end
    else
        error('input file line does not contain mandatory %s field\n\t%s\n',fieldname,tline);
    end

    fieldname = 'DIFF_REV';
    if contains( upper(tline), strcat(fieldname,'=') )
        
        if isempty( InputStruct(ns).PE_REV )
            error('you specified reverse phase encoding data but not the direction! Need to supply a valid PE_REV field!')
        end

        isrt = strfind( upper(tline),[fieldname,'='] ) + (numel(fieldname)+1);
        iend = ispaces(ispaces>=isrt);
        filestring_temp = tline(isrt:iend(1)); % take first ending
        if  ~contains(filestring_temp,',')
             filestring_temp = {filestring_temp};
        else filestring_temp = regexp(filestring_temp,',','split'); 
        end
        InputStruct(ns).N_diff_rev = numel( filestring_temp ); % number of functional runs
        for nr=1:InputStruct(ns).N_diff_rev
            InputStruct(ns).drrun(nr).DIFF_REV_filename = filestring_temp{nr};
        end
    else
        warning('input file line does not contain optional %s field - no epi distortion correction possible\n\t%s\n',fieldname,tline);
        InputStruct(ns).N_diff_rev = 0;
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
            InputStruct(ns).arun(nr).ZCLIP_thr  = NaN;
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
