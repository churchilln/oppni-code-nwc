function output = awarp_HEAD()

%------------------ Model Attributes (mandatory fields) ------------------%
output.attributes.model_name    = 'AWARP';
output.attributes.pipe_stage    = 'Warp';
output.attributes.description   = 'anatomical image warping to template';
output.attributes.extra_fields_number    = 0;
output.attributes.extra_fields_descript  = [];

disp('no inputs - returning attributes');
%-------------------------------------------------------------------------%