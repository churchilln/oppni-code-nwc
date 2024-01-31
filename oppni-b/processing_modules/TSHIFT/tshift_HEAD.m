%
% -------------------------------------------------------------------------
%          module: TSHIFT
% -------------------------------------------------------------------------
%  pipeline stage: P1
%      applied to: fMRI data
%         summary: adjusts for delay in timing offsets between axial slices by "shifting" (interpolating)
%   function call: tshift_<name>( Funcfile, prefix, odir, tpatt, ParamCell )
%                       Funcfile = path+full name of uncorrected fMRI file
%                         prefix = prefix given to the corrected file (output)
%                           odir = directory to store corrected file (output)
%                          tpatt = string specifying  the slice-timing pattern
%                     (ParamCell = supplemental parameters)
%
