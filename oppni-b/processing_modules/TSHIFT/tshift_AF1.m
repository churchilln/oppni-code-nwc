function tshift_AF1( Funcfile, prefix, odir, tpatt )

% odir => sprintf('%s/prewarp',opath2f)
% prefix => sprintf('func%u',nr)
% Funcfile => sprintf('%s/func%u%s.nii',opath0,nr,drop_tag)

if ~exist(sprintf('%s/%s_tshift.nii.gz',odir,prefix),'file')

    unix(sprintf('3dTshift -prefix %s/%s_tshift.nii.gz -tpattern %s %s',odir,prefix,tpatt,Funcfile));
else
    disp('afni-tshift already exists!')
end
