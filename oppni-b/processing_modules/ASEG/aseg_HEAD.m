function output = aseg_HEAD()

%------------------ Model Attributes (mandatory fields) ------------------%
output.attributes.model_name    = 'ASEG';
output.attributes.pipe_stage    = 'Seg';
output.attributes.description   = 'anatomical image segmentation into GM/WM/CSF';
output.attributes.extra_fields_number    = 0;
output.attributes.extra_fields_descript  = [];

disp('no inputs - returning attributes');
%-------------------------------------------------------------------------%