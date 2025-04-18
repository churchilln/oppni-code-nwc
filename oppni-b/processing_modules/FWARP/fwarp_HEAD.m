%
% -------------------------------------------------------------------------
%          module: FWARP
% -------------------------------------------------------------------------
%  pipeline stage: P1
%      applied to: fMRI data
%         summary: warping of functional data into anatomical template space
%   function call: fwarp_<name>( Funcfile_set, prefix_set, odir1, odir2, base_set, Anatloc, ParamCell )
%                   Funcfile_set = cell array containing path+full name of each uncorrected fMRI file run
%                     prefix_set = cell array containing prefix given to each of the corrected files (output)
%                          odir1 = directory to store first set of processing files (output)
%                          odir2 = directory to store final processed datas (output)
%                     (ParamCell = supplemental parameters)


function output = fwarp_HEAD()

%------------------ Model Attributes (mandatory fields) ------------------%
output.attributes.model_name    = 'FWARP';
output.attributes.pipe_stage    = 'P1';
output.attributes.description   = 'functional image warping to template';
output.attributes.extra_fields_number    = 0;
output.attributes.extra_fields_descript  = [];

disp('no inputs - returning attributes');
%-------------------------------------------------------------------------%