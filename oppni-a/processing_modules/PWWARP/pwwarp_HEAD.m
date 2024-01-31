%
% -------------------------------------------------------------------------
%          module: PWWARP
% -------------------------------------------------------------------------
%  pipeline stage: P1
%      applied to: fMRI data
%         summary: warping of functional data into anatomical template space
%   function call: pwwarp_<name>( Funcfile_set, prefix_set, odir1, odir2, Basefile_masked, Anatloc, ParamCell )
%                   Funcfile_set = cell array containing path+full name of each uncorrected fMRI file run
%                     prefix_set = cell array containing prefix given to each of the corrected files (output)
%                          odir1 = directory to store first set of processing files (output)
%                          odir2 = directory to store final processed datas (output)
%                     (ParamCell = supplemental parameters)
