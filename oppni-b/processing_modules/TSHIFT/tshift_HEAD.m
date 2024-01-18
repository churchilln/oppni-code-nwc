function output = tshift_HEAD()

%------------------ Model Attributes (mandatory fields) ------------------%
output.attributes.model_name    = 'TSHIFT';
output.attributes.pipe_stage    = 'P1';
output.attributes.description   = 'slice timing correction by "shifting" slices in time';
output.attributes.extra_fields_number    = 0;
output.attributes.extra_fields_descript  = [];

disp('no inputs - returning attributes');
%-------------------------------------------------------------------------%