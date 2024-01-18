function output = smooth_HEAD()

%------------------ Model Attributes (mandatory fields) ------------------%
output.attributes.model_name    = 'SMOOTH';
output.attributes.pipe_stage    = 'P1';
output.attributes.description   = 'spatial smoothing of functional data';
output.attributes.extra_fields_number    = 1;
output.attributes.extra_fields_descript  = 'scale';

disp('no inputs - returning attributes');
%-------------------------------------------------------------------------%