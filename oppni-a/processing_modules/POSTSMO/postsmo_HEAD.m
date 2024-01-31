%
% -------------------------------------------------------------------------
%          module: SMOOTH
% -------------------------------------------------------------------------
%  pipeline stage: P1
%      applied to: fMRI data
%         summary: spatially smooths data to improve signal-to-noise ratio
%   function call: smooth_<name>( Funcfile, prefix, odir, ParamCell )
%                       Funcfile = path+full name of uncorrected fMRI file
%                         prefix = prefix given to the corrected file (output)
%                           odir = directory to store corrected file (output)
%                     (ParamCell = supplemental parameters)
%
