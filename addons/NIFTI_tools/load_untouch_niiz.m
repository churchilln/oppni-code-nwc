function nii = load_untouch_niiz(filename)
%
% NWC wrapper --> checks if gzip specified, if yes, unzips for reading then rezips

if ~exist(filename,'file')
    error('file as requested does not exist: %s',filename);
end

if strcmp(filename(end-6:end),'.nii.gz')
    unix(sprintf('gunzip %s',filename));
    if ~exist(filename(1:end-3),'file')
        error('load_untouch_niiz failed in file %s: failed creating unzipped nii',filename);
    end
    nii=load_untouch_nii(filename(1:end-3));
    unix(sprintf('gzip %s',filename(1:end-3)));
    if ~exist(filename,'file')
        error('load_untouch_niiz failed in file %s: failed creating rezipped nii',filename);
    end
elseif strcmp(filename(end-3:end),'.nii')
    nii=load_untouch_nii(filename);
else
    error('loadnii - unrecognized suffix!')
end
