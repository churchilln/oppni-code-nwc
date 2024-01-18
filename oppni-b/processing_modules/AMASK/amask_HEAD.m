function output = amask_HEAD()

%------------------ Model Attributes (mandatory fields) ------------------%
output.attributes.model_name    = 'AMASK';
output.attributes.pipe_stage    = 'Warp';
output.attributes.description   = 'anatomical image masking/skull stripping';
output.attributes.extra_fields_number    = 0;
output.attributes.extra_fields_descript  = [];

disp('no inputs - returning attributes');
%-------------------------------------------------------------------------%
