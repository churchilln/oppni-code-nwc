function postsmo_AF1( Funcfile, prefix, odir, ParamCell )

fwhm = ParamCell{1};

% odir => opath4f
% prefix => sprintf('func%u',nr)
% Funcfile => sprintf('%s/postwarp/func%u_warped.nii',opath3f,nr)
% fwhm => PipeStruct_aug.SMOOTH{2}

if ~exist(sprintf('%s/%s_postsmo.nii.gz',odir,prefix),'file')
    
    unix(sprintf('3dmerge -prefix %s/%s_postsmo.nii.gz -doall -1blur_fwhm %s %s',odir,prefix,fwhm,Funcfile));
else
    disp('afni-pre-smoothe already exists!')
end
