function output = GCONN( datacell, params )
%
% =========================================================================
% MODULE_GCONN: module that performs global connectivity analysis
% =========================================================================
%
%   Syntax:
%           output = module_GCONN( datamat, split_info )
%
% ------------------------------------------------------------------------%

%------------------ Model Attributes (mandatory fields) ------------------%
output.attributes.model_name    = 'GCONN';
output.attributes.design_type   = 'nocontrast';
output.attributes.model_type    = 'univariate';
output.attributes.image_types   = 'Gconn,zGconn';
output.attributes.mat2d_types   = [];
output.attributes.stats_types   = [];
output.attributes.uses_roifile  = 0;
output.attributes.uses_taskfile = 0;

if(nargin==0)
    disp('no inputs - returning attributes');
    return;
end
%----------------------- Default Parameter Checks ------------------------%
%-------------------------------------------------------------------------%

Nr = numel(datacell);

for ir=1:Nr

    datamat = datacell{ir};

    % unit normalized on timeseries
    ximag = datamat;
    ximag = bsxfun(@minus,  ximag, mean(ximag,2));
    ximag = bsxfun(@rdivide,ximag, sqrt(sum(ximag.^2,2))+eps);
    [u,l,v] = svd(ximag,'econ');

    % global connectivity
    output.image.Gconn(:,ir) = sqrt( sum(bsxfun(@times,diag(l.^2)',u).^2,2)./size(datamat,1) );
    output.image.zGconn(:,ir) = zscore(output.image.Gconn);
end

% averaging across runs
output.image.Gconn(:,ir)  = mean(output.image.Gconn(:,ir),3);
output.image.zGconn(:,ir)  = mean(output.image.zGconn(:,ir),3);

