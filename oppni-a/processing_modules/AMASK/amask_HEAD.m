%
% -------------------------------------------------------------------------
%          module: AMASK
% -------------------------------------------------------------------------
%  pipeline stage: Warp
%      applied to: Anatomical data
%         summary: masking / skull-stripping of anatomical image 
%   function call: amask_<name>(Adataset, Basedset, odir, ParamCell)
%                       Adataset = path+full name of input anatomical file
%                       Basedset = path+full name of reference template
%                           odir = directory to store corrected file (output)
%                     (ParamCell = supplemental parameters)
%