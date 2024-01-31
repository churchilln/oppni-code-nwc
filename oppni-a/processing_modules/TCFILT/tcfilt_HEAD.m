%
% -------------------------------------------------------------------------
%          module: TCFILT
% -------------------------------------------------------------------------
%  pipeline stage: P1
%      applied to: perf data
%         summary: removal of outlier perf data, correcting for abrupt head movements 
%   function call: tcfilt_<name>( Funcfile, prefix, odir, base, ParamCell )
%                       Funcfile = path+full name of uncorrected fMRI file
%                         prefix = prefix given to the corrected file (output)
%                           odir = directory to store corrected file (output)
%                           base = integer specifying reference volume for motion estimation 
%                     (ParamCell = supplemental parameters)
%
