function output = despike_HEAD()

%------------------ Model Attributes (mandatory fields) ------------------%
output.attributes.model_name    = 'DESPIKE';
output.attributes.pipe_stage    = 'P1';
output.attributes.description   = 'removal of outlier volumes ("spikes") in functional data';
output.attributes.extra_fields_number    = 1;
output.attributes.extra_fields_descript  = 'type';

disp('no inputs - returning attributes');
%-------------------------------------------------------------------------%