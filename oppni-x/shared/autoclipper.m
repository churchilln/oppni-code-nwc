function zval = autoclipper( anatfile_in, HtMeanSd, alpha )
% *
% . quick script to iteratively trim neck and thereby identify
% . minimum z-axis plane for image clipping
%

if ~isempty(dir('__opptmp_autoclip_anat*'))
    error('cannot run autoclip - tempfiles suggest already in progress and/or early termination');
end
if nargin<2
    % based on sample of ~100 participants, manually masked n ac-pc aligned, mean~136.5 mm / sd~5.4
    HtMeanSd = [136.5 5.4];    
end
if nargin<3
    alpha    = 0.005;
end

% takes zipped or unzipped -- defaults to gzip
if contains( anatfile_in, '.nii.gz' )
    unix(sprintf('cp %s __opptmp_autoclip_anat.nii.gz',anatfile_in));
elseif contains( anatfile_in, '.nii' )
    unix(sprintf('cp %s __opptmp_autoclip_anat.nii',anatfile_in));
    unix('gzip __opptmp_autoclip_anat.nii')
else
    error('Unrecognized datatype of file: %s for autoclip\n',anatfile_in)
end

% iteration1
unix('bet __opptmp_autoclip_anat.nii.gz __opptmp_autoclip_anat_seg.nii.gz -m');
M=load_untouch_niiz('__opptmp_autoclip_anat_seg_mask.nii.gz');
f = find( squeeze(sum(sum(M.img,1),2))>1);
cax = M.hdr.hist.qoffset_z + f(1)*M.hdr.dime.pixdim(4);
hsz = (f(end)-f(1))*M.hdr.dime.pixdim(4);
zsz = (hsz - HtMeanSd(1))/HtMeanSd(2);
if normcdf( -abs(zsz) )*2 < alpha
    warning('autoclip firstpass -- found abnormal size of h=%d mm, outside of normal range at z=%f\n',hsz,zsz);
end
unix(sprintf('@clip_volume -input __opptmp_autoclip_anat.nii.gz -below %.02f -prefix __opptmp_autoclip_anat_clp1.nii.gz', cax(1) ));
unix('rm __opptmp_autoclip_anat_seg.nii* __opptmp_autoclip_anat_seg_mask.nii*');

% iteration2
unix(sprintf('bet __opptmp_autoclip_anat_clp1.nii.gz __opptmp_autoclip_anat_seg.nii.gz -m'));
M=load_untouch_niiz('__opptmp_autoclip_anat_seg_mask.nii.gz');
f = find( squeeze(sum(sum(M.img,1),2))>1);
cax = M.hdr.hist.qoffset_z + f(1)*M.hdr.dime.pixdim(4);
hsz = (f(end)-f(1))*M.hdr.dime.pixdim(4);
zsz = (hsz - HtMeanSd(1))/HtMeanSd(2);
if normcdf( -abs(zsz) )*2 < alpha
    warning('autoclip secondpass -- found abnormal size of h=%d mm, outside of normal range at z=%f\n',hsz,zsz);
end
unix(sprintf('@clip_volume -input __opptmp_autoclip_anat.nii.gz -below %.02f -prefix __opptmp_autoclip_anat_clp2.nii.gz', cax(1) ));
unix('rm __opptmp_autoclip_anat_seg.nii* __opptmp_autoclip_anat_seg_mask.nii*');

% iteration3
unix(sprintf('bet __opptmp_autoclip_anat_clp2.nii.gz __opptmp_autoclip_anat_seg.nii.gz -m'));
M=load_untouch_niiz('__opptmp_autoclip_anat_seg_mask.nii.gz');
f = find( squeeze(sum(sum(M.img,1),2))>1);
cax = M.hdr.hist.qoffset_z + f(1)*M.hdr.dime.pixdim(4);
hsz = (f(end)-f(1))*M.hdr.dime.pixdim(4);
zsz = (hsz - HtMeanSd(1))/HtMeanSd(2);
if normcdf( -abs(zsz) )*2 < alpha
    error('autoclip lastpass -- found abnormal size of h=%d mm, outside of normal range at z=%f\n',hsz,zsz);
end
unix(sprintf('@clip_volume -input __opptmp_autoclip_anat.nii.gz -below %.02f -prefix __opptmp_autoclip_anat_clp3.nii.gz', cax(1) ));
unix('rm __opptmp_autoclip_anat_seg.nii* __opptmp_autoclip_anat_seg_mask.nii*');

zval = cax(1) - 2; % drop 2mm below threshold ... seems to work consistently with anat data

% tempfile cleanup
unix('rm -r __opptmp_autoclip_anat*');
