function undist_AF1( Funcfile, prefix, odir, acqpar, refvol_cell, ParamCell )
%
% .undist_AF1:
% .distortion correction using AFNI utilities
% .uses afni_proc.py blip correction

% odir => sprintf('%s/prewarp',opath2f)
% prefix => sprintf('func%u',nr)
% Funcfile => sprintf('%s/prewarp/func%u_tshift.nii.gz',opath2f,nr)

outfile = sprintf('%s/%s_undist.nii.gz',odir,prefix);

if ~exist(outfile,'file')

    if numel(refvol_cell)<1 || isempty(refvol_cell{1})
        error('UNDIST AF1 needs one reverse phase-encode volume in DIST');
    end
    pa_file = refvol_cell{1};

    CODE_PATH = fileparts(which('P0_fmri_populateDirectories.m'));
    python_script = fullfile(CODE_PATH,'scriptheap','distortion','run_afni_blip.py');
    subj_id = sprintf('%s_afni_blip',regexprep(prefix,'[^A-Za-z0-9_]','_'));
    pref = fullfile(odir,[prefix '_undist_AF1']);

    cmd = sprintf(['python3 "%s" ' ...
        '--ap "%s" ' ...
        '--pa "%s" ' ...
        '--subj-id "%s" ' ...
        '--outdir "%s"'], ...
        python_script,Funcfile,pa_file,subj_id,pref);

    [status,result] = system(cmd);

    if status ~= 0
        error('UNDIST AF1 failed:\n%s',result);
    end

    corrected_file = fullfile(pref,[subj_id '.results'],sprintf('pb02.%s.r01.blip+orig',subj_id));
    if ~exist([corrected_file '.HEAD'],'file')
        error('expected UNDIST AF1 output not found:\n\t%s\n',corrected_file);
    end

    cmd = sprintf('3dAFNItoNIFTI -prefix "%s" "%s"',outfile,corrected_file);
    [status,result] = system(cmd);
    if status ~= 0
        error('UNDIST AF1 conversion failed:\n%s',result);
    end

    if ~exist(outfile,'file')
        error('failed to create UNDIST output:\n\t%s\n',outfile);
    end
else
    disp('afni-blip undist already exists!')
end
