function undist_FL2( Funcfile, prefix, odir, acqpar, refvol_cell, ParamCell )
%
% .undist_FL2:
% .distortion correction using FSL utilities
% .uses fieldmap preparation and fugue modules

% odir => sprintf('%s/prewarp',opath2f)
% prefix => sprintf('func%u',nr)
% Funcfile => sprintf('%s/prewarp/func%u_tshift.nii.gz',opath2f,nr)

outfile = sprintf('%s/%s_undist.nii.gz',odir,prefix);

if ~exist(outfile,'file')

    if numel(refvol_cell)<2 || isempty(refvol_cell{1}) || isempty(refvol_cell{2})
        error('UNDIST FL2 needs magnitude and phasediff volumes in DIST');
    end
    mag_file = refvol_cell{1};
    phasediff_file = refvol_cell{2};

    ap_pe  = require_dist_field(acqpar,'AP_PE_DIR');
    echosp = require_dist_field(acqpar,'ECHOSP');
    te1    = require_dist_field(acqpar,'TE1');
    te2    = require_dist_field(acqpar,'TE2');

    CODE_PATH = fileparts(which('P0_fmri_populateDirectories.m'));
    python_script = fullfile(CODE_PATH,'scriptheap','distortion','run_fieldmap_fsl.py');
    pref = fullfile(odir,[prefix '_undist_FL2']);

    cmd = sprintf(['python3 "%s" ' ...
        '--ap "%s" ' ...
        '--mag "%s" ' ...
        '--phasediff "%s" ' ...
        '--phase-encoding-direction "%s" ' ...
        '--effective-echo-spacing "%s" ' ...
        '--echo-time-1 "%s" ' ...
        '--echo-time-2 "%s" ' ...
        '--outdir "%s"'], ...
        python_script,Funcfile,mag_file,phasediff_file,ap_pe,echosp,te1,te2,pref);

    [status,result] = system(cmd);

    if status ~= 0
        error('UNDIST FL2 failed:\n%s',result);
    end

    corrected_file = fullfile(pref,'AP_fieldmap_corrected.nii.gz');
    if ~exist(corrected_file,'file')
        error('expected UNDIST FL2 output not found:\n\t%s\n',corrected_file);
    end

    copyfile(corrected_file,outfile);
    if ~exist(outfile,'file')
        error('failed to create UNDIST output:\n\t%s\n',outfile);
    end
else
    disp('fsl-fieldmap undist already exists!')
end
end

function fieldval = require_dist_field(acqpar,fieldname)
if ~isfield(acqpar,fieldname) || isempty(acqpar.(fieldname))
    error('UNDIST FL2 needs %s in DIST parameter file',fieldname);
end
fieldval = acqpar.(fieldname);
end
