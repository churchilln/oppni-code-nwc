function analysis_model_out = check_fmri_analysis_model( analysis_model )
%

% declaring path
CODE_PATH = fileparts(which('P0_fmri_populateDirectories.m'));
if CODE_PATH(end)~='/'
    CODE_PATH = [CODE_PATH '/'];
end

if( strcmpi(analysis_model,'NONE') )
    
    %-- no analysis being performed
    analysis_model_out.model_name   = 'NONE';
    analysis_model_out.design_type  = 'nocontrast';
    analysis_model_out.metric_def   = 'no';
    analysis_model_out.uses_roifile  = 0;
    analysis_model_out.uses_taskfile = 0;
else
    
    module_dir = [CODE_PATH 'analysis_modules'];
    
    e = dir( module_dir );
    kq=0;
    for(i=1:length(e))
       if( ~isempty(strfind(e(i).name,'.m')) )
          kq=kq+1;
          [path module_list{kq} ext] = fileparts( e(i).name );
       end
    end

    ix = find( strcmpi( analysis_model, module_list ) );

    if( isempty(ix) )
        module_csep=[];
        for(i=1:numel(module_list)) module_csep=[module_csep,', ',module_list{i}]; end
        error('The analysis model %s not found among module list:\n\n%s\n',analysis_model,module_csep(2:end));
    else
        analysis_model = module_list{ix};
    end

    currPath=pwd;                   % get current path
    cd(module_dir);                 % jump to module directory
    p   = str2func(analysis_model); % get function handle
    cd(currPath);                   % jump back to current path
    tmp = p();                      % storing attributes
    analysis_model_out          = tmp.attributes;
    analysis_model_out.filepath = module_dir; % store absolute file path
end
