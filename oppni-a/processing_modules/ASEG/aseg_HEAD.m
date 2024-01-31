%
% -------------------------------------------------------------------------
%          module: ASEG
% -------------------------------------------------------------------------
%  pipeline stage: Seg
%      applied to: Anatomical data
%         summary: segments anatomical image into GM/WM/CSF
%   function call: aseg_<name>( Anatmasked, Mask, odir, ParamCell )
%                     Anatmasked = path+full name of masked anatomical image 
%                           Mask = path+full name of binary brain mask
%                           odir = directory to store corrected file (output)
%                     (ParamCell = supplemental parameters)
%
