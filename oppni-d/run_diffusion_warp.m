function run_diffusion_warp(subject_dir, python_script, template_fa)
%
% . run diffusion scalar warping with TBSS-cleaned FA and ANTs
%

if ~exist(subject_dir,'dir')
    error('subject directory not found:\n\t%s\n',subject_dir);
end
if ~exist(python_script,'file')
    error('python diffusion warp script not found:\n\t%s\n',python_script);
end
if ~exist(template_fa,'file')
    error('FA template not found:\n\t%s\n',template_fa);
end

python_exec = getenv('OPPNI_PYTHON');
if isempty(python_exec)
    python_exec = 'python3';
end

% construct python command
threads = getenv('OPPNI_DIFF_WARP_THREADS');
cmd = sprintf('%s "%s" --subject-dir "%s" --template "%s"', ...
    python_exec,python_script,subject_dir,template_fa);
if ~isempty(threads)
    cmd = sprintf('%s --threads "%s"',cmd,threads);
end

fprintf('\nrunning diffusion warp:\n%s\n\n',cmd);

% execute python command and report failure logs
[status,result] = system(cmd,'-echo');
if status ~= 0
    error('diffusion warp failed for:\n\t%s\n\n%s',subject_dir,result);
end

fprintf('\ndiffusion warp complete.\n');
