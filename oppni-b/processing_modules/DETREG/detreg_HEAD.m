function output = detreg_HEAD()

%------------------ Model Attributes (mandatory fields) ------------------%
output.attributes.model_name    = 'DETREG';
output.attributes.pipe_stage    = 'P2';
output.attributes.description   = 'detrending regression - removal of low-frequency signal fluctuations';
output.attributes.extra_fields_number    = 1;
output.attributes.extra_fields_descript  = 'order';

disp('no inputs - returning attributes');
%-------------------------------------------------------------------------%