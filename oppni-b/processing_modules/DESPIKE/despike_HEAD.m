%
% -------------------------------------------------------------------------
%          module: DESPIKE
% -------------------------------------------------------------------------
%  pipeline stage: P1
%      applied to: fMRI data
%         summary: removal of outlier brain volumes, correcting for abrupt head movements 
%   function call: despike_<name>( Funcfile, prefix, odir, base, ParamCell )
%                       Funcfile = path+full name of uncorrected fMRI file
%                         prefix = prefix given to the corrected file (output)
%                           odir = directory to store corrected file (output)
%                           base = integer specifying reference volume for motion estimation 
%                     (ParamCell = supplemental parameters)
%