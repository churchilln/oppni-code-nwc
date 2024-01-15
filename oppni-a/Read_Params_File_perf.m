function ParamStruct = Read_Params_File_perf( paramsfile )
% *PERF
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

% PERF_MODEL -> code module used to run perfusion
% CALIB_MATCH -> matching calib scans by run, if multi-run data?
% ANALYSIS/CONTRAST/VOXRES/TEMPLATE --> the usual
% KM.. --> kinetic modelling parameters
steplist = {'PERF_MODEL','SUBTRACT_TYPE','ANALYSIS','CONTRAST','VOXRES','TEMPLATE','INIMOT','ROIMASK','KM_BLOODT1','KM_BLOODT2S','KM_PARCOEF','KM_PDRAT','KM_LABEFF','KM_QTI'};
mandatory= [ 1,          1,              1,         1,         1,       1,         1,       1,        0,           0,            0,           0,         0,          0,      ];

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
        ParamStruct.(steplist{s}) = paramstring((bleft+1):(bright-1));
    end
end
