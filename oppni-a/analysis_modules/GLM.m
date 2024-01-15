function output = GLM( datamat, task_info, params )
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
output.attributes.num_comp      = 'multi_component';
output.attributes.spms          = 'rSPM,Beta';
output.attributes.metrics       = 'R';
output.attributes.uses_roifile  = 0;
output.attributes.uses_taskfile = 1;

if(nargin==0)
    disp('no inputs - returning attributes');
    return;
end
%----------------------- Default Parameter Checks ------------------------%
%-------------------------------------------------------------------------%

% unit normalized bold timeseries
ximag = (datamat - mean(datamat,2))./(std(datmat,0,2)+eps);
% full design matrix - no mean, since matrix is centered
xdes  = task_info.design_mat;
xdes  = bsxfun(@rdivide, xdes, sum(xdes,1)+eps);


% raw betas
output.image.beta = (dataVol * xdes) /( xdes'*xdes );

for i=1:numel(task_info.contrast)

    % 
    bsum1 = mean( output.image.beta(:,task_info.contrast(i).c1) ,2);
    if isempty(task_info.contrast(i).c2)
        bsum2 = 0;
    else
        bsum2 = mean( output.image.beta(:,task_info.contrast(i).c2) ,2);
    end
    output.image.beta_contrast(:,i) = bsum1-bsum2;
end
