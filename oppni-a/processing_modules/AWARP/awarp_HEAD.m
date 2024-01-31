%
% -------------------------------------------------------------------------
%          module: AWARP
% -------------------------------------------------------------------------
%  pipeline stage: Warp
%      applied to: Anatomical data
%         summary: warping anatomical image to a template
%   function call: awarp_<name>(Adataset, Maskdset, Basedset, odir, ParamCell)
%                       Adataset = path+full name of input anatomical file
%                       Maskdset = path+full name of brain mask
%                       Basedset = path+full name of reference template
%                           odir = directory to store corrected file (output)
%                     (ParamCell = supplemental parameters)
%
