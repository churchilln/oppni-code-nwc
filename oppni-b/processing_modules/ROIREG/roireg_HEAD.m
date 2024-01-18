function output = roireg_HEAD()

%------------------ Model Attributes (mandatory fields) ------------------%
output.attributes.model_name    = 'ROIREG';
output.attributes.pipe_stage    = 'P2';
output.attributes.description   = 'WM/CSF ROI-based regression - removal of residual physiologic effects';
output.attributes.extra_fields_number    = 2;
output.attributes.extra_fields_descript  = 'type,order';

disp('no inputs - returning attributes');
%-------------------------------------------------------------------------%