function SaveAsNIfTI(data, target, output)

% function SaveAsNIfTI(data, nifti)
%
% Input:
%
% data: the data array to be saved to disk
%
% target: the NIfTI object specifying the target volume specification
%
% output: the filename for the output NIfTI file
%

% following the example in
% http://niftilib.sourceforge.net/mat_api_html/README.txt
%
% author: Gary Hui Zhang (gary.zhang@ucl.ac.uk)
%

nii=target;                                                       %[NWC ed]
nii.img = data;                                                   %[NWC ed]
nii.hdr.dime.datatype = 16;                                       %[NWC ed]
nii.hdr.hist = target.hdr.hist;                                   %[NWC ed]
nii.hdr.dime.dim(5) = 1;                                          %[NWC ed]
save_untouch_nii(nii,output);                                     %[NWC ed]

% dat = file_array;                                                [NWC ed]
% dat.fname = output;                                              [NWC ed] 
% dat.dim = target.dim;                                            [NWC ed]
% dat.dtype = 'FLOAT64-LE';                                        [NWC ed]
% dat.offset = ceil(348/8)*8;                                      [NWC ed]
% 
% N = nifti;                                                       [NWC ed]
% N.dat = dat;                                                     [NWC ed]
% N.mat = target.mat;                                              [NWC ed]
% N.mat_intent = target.mat_intent;                                [NWC ed]
% N.mat0 = target.mat0;                                            [NWC ed]
% N.mat0_intent = target.mat0_intent;                              [NWC ed]
% 
% create(N);                                                       [NWC ed]
% 
% N.dat(:,:,:,:) = data;                                           [NWC ed]
