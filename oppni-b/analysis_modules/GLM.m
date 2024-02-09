function output = GLM( datacell, task_struc, params )
%
% =========================================================================
% MODULE_GLM: module that performs General Linear Model analyses
% =========================================================================
%
%   Syntax:
%           output = module_GLM( datamat, split_info )
%
% ------------------------------------------------------------------------%

%------------------ Model Attributes (mandatory fields) ------------------%
output.attributes.model_name    = 'GLM';
output.attributes.design_type   = 'regression';
output.attributes.model_type    = 'univariate';
output.attributes.image_types   = 'BetaCtr';
output.attributes.mat2d_types   = [];
output.attributes.stats_types   = [];
output.attributes.uses_roifile  = 0;
output.attributes.uses_taskfile = 1;

if(nargin==0)
    disp('no inputs - returning attributes');
    return;
end
%----------------------- Default Parameter Checks ------------------------%
%-------------------------------------------------------------------------%

% cell-ify if data matrix is input
if ~iscell(datacell)
    datacell = {datacell};
end
if ~iscell(task_struc)
    task_struc = {task_struc};
end
if numel(datacell) ~= numel(task_struc)
    error('number of datamat/taskfile runs do not match');
end

%-- option1: cat everything, do the analysis

datacat = [];
xdescat = [];
for nr=1:numel(datacell)
   dtmp = 100 * datacell{nr}./mean(datacell{nr},2); % scaling as percent of baseline
   dtmp = dtmp - mean(dtmp,2); % subtract out mean response
   datacat = [datacat, dtmp];
   xdescat = [xdescat; task_struc{nr}.design_mat];
   %
   output.temp.design_mat{nr} = task_struc{nr}.design_mat;
end

% raw betas - no intercept, since mean-subtracted
betamaps = (datacat * xdescat) /( xdescat'*xdescat );

for i=1:numel(task_struc{1}.contrast)

    % 
    bsum1 = mean( betamaps(:,task_struc{1}.contrast(i).c1) ,2);
    if isempty(task_struc{1}.contrast(i).c2)
        bsum2 = 0;
    else
        bsum2 = mean( betamaps(:,task_struc{1}.contrast(i).c2) ,2);
    end
    output.image.BetaCtr(:,i) = bsum1-bsum2;
end
