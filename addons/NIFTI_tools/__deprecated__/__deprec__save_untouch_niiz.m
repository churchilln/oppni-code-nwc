function save_untouch_niiz(nii, filename, overwrite)
%
% NWC wrapper --> checks if gzip specified, if yes, zips results
if nargin<3
    overwrite = 1;
end

if exist(filename,'file')
    if overwrite==0
        warning('save_untouh_niiz skipped - will not overwrite existing file %s\n',filename);
    else
        unix(sprintf('rm %s',filename))
    end
end

if strcmp(filename(end-6:end),'.nii.gz')
    save_untouch_nii(nii,filename(1:end-3));
    unix(sprintf('gzip %s',filename(1:end-3)));
    if ~exist(filename,'file')
        error('save_untouch_niiz failed in file %s: failed creating zipped nii',filename);
    end
elseif strcmp(filename(end-3:end),'.nii')
    save_untouch_nii(nii,filename);
else
    error('savenii - unrecognized suffix!')
end
