function output = gsreg_HEAD()

%------------------ Model Attributes (mandatory fields) ------------------%
output.attributes.model_name    = 'GSREG';
output.attributes.pipe_stage    = 'P2';
output.attributes.description   = 'global signal regression - removal of global activation effects';
output.attributes.extra_fields_number    = 1;
output.attributes.extra_fields_descript  = 'type';

disp('no inputs - returning attributes');
%-------------------------------------------------------------------------%