function undist_FL1( Funcfile, prefix, odir, acqpar, refvol_cell, ParamCell )
%
% .undist_FL1:
% .distortion correction using FSL utilities
% .uses topup and applytopup modules

% odir => sprintf('%s/prewarp',opath2f)
% prefix => sprintf('func%u',nr)
% Funcfile => sprintf('%s/prewarp/func%u_tshift.nii.gz',opath2f,nr)

outfile = sprintf('%s/%s_undist.nii.gz',odir,prefix);

if ~exist(outfile,'file')

    if numel(refvol_cell)<1 || isempty(refvol_cell{1})
        error('UNDIST FL1 needs one reverse phase-encode volume in DIST');
    end
    pa_file = refvol_cell{1};

    ap_pe  = require_dist_field(acqpar,'AP_PE_DIR');
    pa_pe  = require_dist_field(acqpar,'PA_PE_DIR');
    ap_tro = require_dist_field(acqpar,'AP_TRO');
    pa_tro = require_dist_field(acqpar,'PA_TRO');

    CODE_PATH = fileparts(which('P0_fmri_populateDirectories.m'));
    python_script = fullfile(CODE_PATH,'scriptheap','distortion','run_topup_ap_pa.py');
    pref = fullfile(odir,[prefix '_undist_FL1']);

    cmd = sprintf(['python3 "%s" ' ...
        '--ap "%s" ' ...
        '--pa "%s" ' ...
        '--ap-pe-dir "%s" ' ...
        '--pa-pe-dir "%s" ' ...
        '--ap-total-readout-time "%s" ' ...
        '--pa-total-readout-time "%s" ' ...
        '--outdir "%s"'], ...
        python_script,Funcfile,pa_file,ap_pe,pa_pe,ap_tro,pa_tro,pref);

    [status,result] = system(cmd);

    if status ~= 0
        error('UNDIST FL1 failed:\n%s',result);
    end

    corrected_file = fullfile(pref,'AP_topup_corrected.nii.gz');
    if ~exist(corrected_file,'file')
        error('expected UNDIST FL1 output not found:\n\t%s\n',corrected_file);
    end

    copyfile(corrected_file,outfile);
    if ~exist(outfile,'file')
        error('failed to create UNDIST output:\n\t%s\n',outfile);
    end
else
    disp('fsl-topup undist already exists!')
end
end

function fieldval = require_dist_field(acqpar,fieldname)
if ~isfield(acqpar,fieldname) || isempty(acqpar.(fieldname))
    error('UNDIST FL1 needs %s in DIST parameter file',fieldname);
end
fieldval = acqpar.(fieldname);
end
