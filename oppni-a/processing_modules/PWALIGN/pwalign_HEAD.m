%
% -------------------------------------------------------------------------
%          module: PWWARP
% -------------------------------------------------------------------------
%  pipeline stage: P1
%      applied to: fMRI data
%         summary: rigid alignment of perf data to ref-volume
%   function call: pwalign_<name>( Perffile, Basefile, prefix, odir1, ParamCell )
%                   Funcfile_set = cell array containing path+full name of each uncorrected fMRI file run
%                     prefix_set = cell array containing prefix given to each of the corrected files (output)
%                          odir1 = directory to store first set of processing files (output)
%                          odir2 = directory to store final processed datas (output)
%                     (ParamCell = supplemental parameters)
