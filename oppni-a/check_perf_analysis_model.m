function analysis_model_out = check_perf_analysis_model( analysis_model )
%

% declaring path
CODE_PATH = fileparts(which('P0_perf_populateDirectories.m'));
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
    
    module_dir = [CODE_PATH 'analysis_modules/'];
    if ~exist([module_dir,analysis_model,'.m'],'file')
        e = dir( [module_dir, '*.m'] );
        modlist=[];
        for i=1:numel(e)
            modlist = [modlist, ', ',e(i).name(1:end-2)];
        end
        error('The analysis model %s not found among module list:\n\n%s\n\n',analysis_model,modlist(3:end));
    else
        currPath=pwd;                   % get current path
        cd(module_dir);                 % jump to module directory
        p   = str2func(analysis_model); % get function handle
        cd(currPath);                   % jump back to current path
        tmp = p();                      % storing attributes
        analysis_model_out          = tmp.attributes;
        analysis_model_out.filepath = module_dir; % store absolute file path
    end
end
