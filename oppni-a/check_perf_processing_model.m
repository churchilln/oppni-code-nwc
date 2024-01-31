function pipeline_model_out = check_perf_processing_model( pipestruct )
%

% declaring path
CODE_PATH = fileparts(which('P0_perf_populateDirectories.m'));
if CODE_PATH(end)~='/'
    CODE_PATH = [CODE_PATH '/'];
end

% pipe steps
steplist = {'PNAME','AMASK','AWARP','ASEG','TCFILT','PWALIGN','PRESMO','PWWARP','POSTSMO'}; % ==> perf_est after presmo

% append unique pipeline name

for i=2:numel(steplist) % for each processing step

    if strcmpi(pipestruct.(steplist{i}){1},'OFF')
        pipeline_model_out.filepath = [];
        pipeline_model_out.model_name = 'OFF';
    else
        module_dir = [CODE_PATH 'processing_modules/' steplist{i},'/'];
        pa_model   = [lower(steplist{i}),'_',pipestruct.(steplist{i}){1}];
    
        if ~exist([module_dir,pa_model,'.m'],'file')
            e = dir( [module_dir, '*.m'] );
            modlist=[];
            for i=1:numel(e)
                modlist = [modlist, ', ',e(i).name(1:end-2)];
            end
            error('The analysis model %s not found among module list:\n\n%s\n\n',pa_model,modlist(3:end));
        else
    %         currPath=pwd;                   % get current path
    %         cd(module_dir);                 % jump to module directory
    %         p   = str2func(pa_model); % get function handle
    %         cd(currPath);                   % jump back to current path
    %         tmp = p();                      % storing attributes
    %         analysis_model_out          = tmp.attributes;
    %         analysis_model_out.filepath = module_dir; % store absolute file path
            pipeline_model_out.(steplist{i}).filepath   = module_dir;
            pipeline_model_out.(steplist{i}).model_name = pa_model;

%             currPath=pwd;                   % get current path
%             cd(module_dir);                 % jump to module directory
%             p   = str2func([lower(steplist{i}),'_HEAD']); % get function handle
%             cd(currPath);                   % jump back to current path
%             tmp = p();                      % storing attributes
% 
%             if numel(pipestruct.(steplist{i})) ~= (tmp.attributes.extra_fields_number + 1)
%                 error('pipeline step %s is missing extra fields %s\n',steplist{i},tmp.attributes.extra_fields_descript);
%             end
        end
    end
end
