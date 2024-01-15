function output = CONN( datacell, seed_struc, params )
%
% =========================================================================
% MODULE_CONN: module that performs seed-based connectivity analysis
% =========================================================================
%
%   Syntax:
%           output = module_CONN( datamat, split_info )
%
% ------------------------------------------------------------------------%

%------------------ Model Attributes (mandatory fields) ------------------%
output.attributes.model_name    = 'CONN';
output.attributes.design_type   = 'connectivity';
output.attributes.model_type    = 'univariate';
output.attributes.image_types   = 'SeedCorr';
output.attributes.mat2d_types   = 'ConnMat';
output.attributes.stats_types   = [];
output.attributes.uses_roifile  = 1;
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
    ximag = (datamat - mean(datamat,2))./(std(datamat,0,2)+eps);
    xseed = seed_struc.seedmat;
    
    % expand case of multiparc data
    if strcmpi( seed_struc.contrast(1).type, 'multi') && size(xseed,2)==1

        labvals = sort(unique(seedmat(seedmat>0)),1);
        ns2 = numel(labvals);
        xseed_new = zeros(size(xseed,1),ns2);
        for i=1:ns2
            xseed_new(:,i) = double(seedmat==labvals(i));
        end
        xseed = xseed_new; clear xseed_new;
    end
    
    % seed timecourse
    tt  = zscore( ximag'*xseed );
    % correlation seed map
    output.image.SeedCorr(:,:,ir) = (ximag * tt)./(size(ximag,2)-1);
    % network connectivity matrix
    output.mat2d.ConnMat(:,:,ir) = corr(tt);

end

output.image.SeedCorr = mean(output.image.SeedCorr,3);
output.mat2d.ConnMat = mean(output.mat2d.ConnMat,3);