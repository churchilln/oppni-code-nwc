function output = lopass_HEAD()

%------------------ Model Attributes (mandatory fields) ------------------%
output.attributes.model_name    = 'LOPASS';
output.attributes.pipe_stage    = 'P2';
output.attributes.description   = 'low-pass filtering of functional data - removing high-freuqency fluctuations';
output.attributes.extra_fields_number    = 0;
output.attributes.extra_fields_descript  = [];

disp('no inputs - returning attributes');
%-------------------------------------------------------------------------%