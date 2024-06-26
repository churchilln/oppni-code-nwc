function CreateROI(dwifile, maskfile, outputfile)
%
% function CreateROI(dwifile, maskfile)
%
% This function converts 4-D DWI volume into the data format suitable
% for subsequent NODDI fitting.
%
% Inputs:
%
% dwifile: the 4-D DWI volume in NIfTI or Analyze format
%
% maskfile: the brain mask volume in NIfTI or Analyze format which
% specifies the voxels to include for fitting
%
% outputfile: the mat file to store the resulting data for subsequent
% fitting
%
% author: Gary Hui Zhang (gary.zhang@ucl.ac.uk)
%

%%% [NWC edits]: commented below where utilities in SPM toolbox (niftimatlib-1.2)
%%%               are replaced with NIFTI-matlab utilities

% % % first check if niftimatlib is available
% % if (exist('nifti') ~= 2)
% %     error('niftimatlib does not appear to be installed or included in the search path');
% %     return;
% % end

% load the DWI volume
fprintf('loading the DWI volume : %s\n', dwifile);
dwi = load_untouch_nii(dwifile); % = nifti(dwifile);               [NWC ed]
xsize = dwi.hdr.dime.dim(2); % = dwi.dat.dim(1);                   [NWC ed]
ysize = dwi.hdr.dime.dim(3); % = dwi.dat.dim(2);                   [NWC ed]
zsize = dwi.hdr.dime.dim(4); % = dwi.dat.dim(3);                   [NWC ed]
ndirs = dwi.hdr.dime.dim(5); % = dwi.dat.dim(4);                   [NWC ed]

% convert the data from scanner order to voxel order
dwi = dwi.img; % = dwi.dat(:,:,:,:); [NWC ed]
dwi = permute(dwi,[4 1 2 3]);

% load the brain mask volume
fprintf('loading the brain mask : %s\n', maskfile);
mask = load_untouch_nii(maskfile); % = nifti(maskfile);            [NWC ed]
mask = mask.img(:,:,:); % = mask.dat(:,:,:);                       [NWC ed]

% create an ROI that is in voxel order and contains just the voxels in the
% brain mask
fprintf('creating the output ROI ...\n');

% first get the number of voxels first
% to more efficiently allocate the memory
count=0;
for i=1:xsize
    for j=1:ysize
        for k=1:zsize
            if mask(i,j,k) > 0
                count = count + 1;
                mask(i,j,k) = count;
            end
        end
    end
end
roi = zeros(count,ndirs);
idx = zeros(count,3);

% next construct the ROI
count=0;
for i=1:xsize
    for j=1:ysize
        for k=1:zsize
            if mask(i,j,k) > 0
                count = count + 1;
                roi(count,:) = dwi(:,i,j,k);
                idx(count,:) = [i j k];
            end
        end
    end
end

% save the ROI
fprintf('saving the output ROI as %s\n', outputfile);
save(outputfile, 'roi', 'mask', 'idx');

disp('done');

