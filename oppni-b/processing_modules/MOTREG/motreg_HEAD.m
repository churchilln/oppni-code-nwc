function output = motreg_HEAD()

%------------------ Model Attributes (mandatory fields) ------------------%
output.attributes.model_name    = 'MOTREG';
output.attributes.pipe_stage    = 'P2';
output.attributes.description   = 'motion parameter regression - removal of residual motion effects';
output.attributes.extra_fields_number    = 2;
output.attributes.extra_fields_descript  = 'type,order';

disp('no inputs - returning attributes');
%-------------------------------------------------------------------------%