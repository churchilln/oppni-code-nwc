function ParamStruct = Read_Params_File_fmri( paramsfile )
%
% ParamStruct.(field) = string

ParamStruct=[];
if ~contains(paramsfile,'=')
    fid   = fopen(paramsfile);
    if fid==-1
        error('cannot open params file:\n\t%s\n',paramsfile);
    end
    newline = fgetl(fid);
    if ~ischar(newline)
        error('params file is empty:\n\t%s\n',paramsfile);
    end
    newline(newline==' ') = [];
    paramstring=[];
    while ischar(newline)
        %
        paramstring = [paramstring ' ' newline];
        newline      = fgetl(fid);
        newline(newline==' ') = [];
    end
    fclose(fid);
else
    paramstring = paramsfile; %% just directly assign the string
end

steplist  = {'ANALYSIS','CONTRAST','VOXRES','TEMPLATE','INIMOT','ROIMASK','TRUNC_ANL','GMMASK_ANL'};
mandatory = [1,         1,          1,       1,        0,       0,        0,          0           ];
procstyle = [0,         0,          0,       0,        1,       1,        0,          0           ]; % -if multiple positional args, in same style as proc

ileft  = strfind( paramstring, '[' );
iright = strfind( paramstring, ']' );
for(s=1:length(steplist))
    iStep = strfind( upper(paramstring), strcat(steplist{s},'=') );
    if isempty(iStep)
        if mandatory(s)==1
            error('Param field %s has not been defined, please check the spelling!\n This field must be specified to continue\n',steplist{s});
        elseif mandatory(s)==0
            warning('Param field %s has not been defined but is optional, please check the spelling!\n Will use defaults\n',steplist{s});
        else
            error('?');
        end
    else
        bleft       = ileft (ileft >iStep);   bleft= bleft(1);
        bright      = iright(iright>iStep);  bright=bright(1);

        if procstyle(s)==0
            ParamStruct.(steplist{s}) = paramstring((bleft+1):(bright-1));
        else
            pipeargs    = regexp(  paramstring((bleft+1):(bright-1))  ,',','split');
            if sum(cellfun(@isempty,pipeargs))>0
                error('pipe step %s has empty positional argument?',steplist{s});
            end
            ParamStruct.(steplist{s}) = pipeargs;
        end
    end
end

%--> append empty field in position 2 if undeclared / for general modularity
for(s=1:length(steplist))
    if procstyle(s)==1 && numel(ParamStruct.(steplist{s}))==1
        ParamStruct.(steplist{s}){2} = [];
    end
end