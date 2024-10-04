function ricor_AF1( Funcfile, prefix, odir, pulsfile, respfile, acqpar, ParamCell)

pref = [odir,'/__opptmp_p2func_ricor'];

spec_case = {'alt+z','alt+z2','alt-z','alt-z2','seq+z','seq-z'};

if ~exist(sprintf('%s/%s_ricor.nii.gz',odir,prefix),'file')

    % build directory struct recursively
    unix(sprintf('mkdir -p %s',pref));

    % embed the slice timing information!
    if sum(strcmpi(acqpar.tpatt,spec_case))>0
        disp('special case arg');
        unix(sprintf('3dTcat -prefix %s/func_tinhdr.nii.gz -tpattern %s %s',pref,acqpar.tpatt,Funcfile));
    else
        disp('custom arg');
        unix(sprintf('3dTcat -prefix %s/func_tinhdr.nii.gz -tpattern @%s %s',pref,acqpar.tpatt,Funcfile));
    end
    % now ricor it
    unix(sprintf('3dretroicor -order 2 -prefix %s/%s_ricor.nii.gz -card %s -resp %s %s/func_tinhdr.nii.gz',odir,prefix,pulsfile,respfile,pref));

    unix(sprintf('rm -rf %s',pref));
else
    disp('afni-ricor already exists!')
end
