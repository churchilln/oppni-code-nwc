function [Xreg,stat] = gsreg_OP1( datamat, ParamCell )
%
% .gsreg_OP1:
% .global signal estimation for denoising
% .implements both global average and PCA approaches

modelstr = ParamCell{1};
stat = [];

if strcmpi(modelstr,'GAV')
    Xreg = mean(zscore(datamat'),2);
elseif strcmpi(modelstr,'PC1')
    [~,~,v] = svd( zscore(datamat')','econ');
    Xreg = v(:,1);
end

stat = [1];