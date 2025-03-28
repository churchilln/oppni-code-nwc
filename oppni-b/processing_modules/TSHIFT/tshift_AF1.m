function tshift_AF1( Funcfile, prefix, odir, tpatt, ParamCell )
%
% .tshift_AF1:
% .slice-timing correction using AFNI utilities
% .uses 3dTshift module

% odir => sprintf('%s/prewarp',opath2f)
% prefix => sprintf('func%u',nr)
% Funcfile => sprintf('%s/func%u%s.nii',opath0,nr,drop_tag)

spec_case = {'alt+z','alt+z2','alt-z','alt-z2','seq+z','seq-z'};

if ~exist(sprintf('%s/%s_tshift.nii.gz',odir,prefix),'file')
    if sum(strcmpi(tpatt,spec_case))>0
        disp('special case arg');
        unix(sprintf('3dTshift -prefix %s/%s_tshift.nii.gz -tpattern %s %s',odir,prefix,tpatt,Funcfile));
    else
        disp('custom arg');
        % slice offsets in msec, from bottom to top
        unix(sprintf('3dTshift -prefix %s/%s_tshift.nii.gz -tpattern @%s %s',odir,prefix,tpatt,Funcfile));
    end
else
    disp('afni-tshift already exists!')
end
