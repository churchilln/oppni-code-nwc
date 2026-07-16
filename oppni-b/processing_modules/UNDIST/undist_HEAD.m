%
% -------------------------------------------------------------------------
%          module: UNDIST
% -------------------------------------------------------------------------
%  pipeline stage: P1
%      applied to: fMRI data
%         summary: corrects susceptibility-related distortions in functional data
%   function call: undist_<name>( Funcfile, prefix, odir, acqpar, refvol_cell, ParamCell )
%                       Funcfile = path+full name of uncorrected fMRI file
%                         prefix = prefix given to the corrected file (output)
%                           odir = directory to store corrected file (output)
%                         acqpar = structure containing distortion parameters
%                    refvol_cell = cell array containing reverse PE or fieldmap files
%                     (ParamCell = supplemental parameters)
%
%
