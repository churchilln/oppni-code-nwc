function output = fwarp_HEAD()

%------------------ Model Attributes (mandatory fields) ------------------%
output.attributes.model_name    = 'FWARP';
output.attributes.pipe_stage    = 'P1';
output.attributes.description   = 'functional image warping to template';
output.attributes.extra_fields_number    = 0;
output.attributes.extra_fields_descript  = [];

disp('no inputs - returning attributes');
%-------------------------------------------------------------------------%