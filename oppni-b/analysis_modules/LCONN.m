function output = LCONN( datacell, params )
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
output.attributes.model_name    = 'LCONN';
output.attributes.design_type   = 'nocontrast';
output.attributes.model_type    = 'univariate';
output.attributes.image_types   = 'Lconn,zLconn';
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

mask = params.mask;


Nr = numel(datacell);

for ir=1:Nr

    datamat = datacell{ir};

    datamat = bsxfun(@minus,datamat,mean(datamat,2));
    datamat = bsxfun(@rdivide,datamat,sqrt(sum(datamat.^2,2)));
    
    TMPVOL = zeros( [size(mask) size(datamat,2)] );
    for t= 1:size(TMPVOL,4)
        tmp = mask;
        tmp(tmp>0) = datamat(:,t);
        TMPVOL(:,:,:,t) = tmp;
    end
    [n1,n2,n3]=size(mask);
    
    corsum=0;
    msksum=0;
    
    % d3
    mask_shf = cat(3, zeros(n1,n2,1), mask(:,:,1:end-1));
    TVOL_shf = cat(3, zeros(n1,n2,1,size(TMPVOL,4)), TMPVOL(:,:,1:end-1,:));
        corsum = corsum + sum(TMPVOL.*TVOL_shf,4);
        msksum = msksum + mask_shf;
    
    mask_shf = cat(3, mask(:,:,2:end), zeros(n1,n2,1));
    TVOL_shf = cat(3, TMPVOL(:,:,2:end,:), zeros(n1,n2,1,size(TMPVOL,4)));
        corsum = corsum + sum(TMPVOL.*TVOL_shf,4);
        msksum = msksum + mask_shf;
    
    % d2
    mask_shf = cat(2, zeros(n1,1,n3), mask(:,1:end-1,:));
    TVOL_shf = cat(2, zeros(n1,1,n3,size(TMPVOL,4)), TMPVOL(:,1:end-1,:,:));
        corsum = corsum + sum(TMPVOL.*TVOL_shf,4);
        msksum = msksum + mask_shf;
    
    mask_shf = cat(2, mask(:,2:end,:), zeros(n1,1,n3));
    TVOL_shf = cat(2, TMPVOL(:,2:end,:,:), zeros(n1,1,n3,size(TMPVOL,4)));
        corsum = corsum + sum(TMPVOL.*TVOL_shf,4);
        msksum = msksum + mask_shf;
    
    % d1
    mask_shf = cat(1, zeros(1,n2,n3), mask(1:end-1,:,:));
    TVOL_shf = cat(1, zeros(1,n2,n3,size(TMPVOL,4)), TMPVOL(1:end-1,:,:,:));
        corsum = corsum + sum(TMPVOL.*TVOL_shf,4);
        msksum = msksum + mask_shf;
    
    mask_shf = cat(1, mask(2:end,:,:), zeros(1,n2,n3));
    TVOL_shf = cat(1, TMPVOL(2:end,:,:,:), zeros(1,n2,n3,size(TMPVOL,4)));
        corsum = corsum + sum(TMPVOL.*TVOL_shf,4);
        msksum = msksum + mask_shf;
    
    % average local connectivity
    spm = corsum./msksum;
    % local connectivity
    output.image.Lconn(:,ir) = spm(mask>0);
    output.image.zLconn(:,ir) = zscore(  output.image.Lconn(:,ir) );

end

% averaging across runs
output.image.Lconn(:,ir)  = mean(output.image.Lconn(:,ir),3);
output.image.zLconn(:,ir)  = mean(output.image.zLconn(:,ir),3);

