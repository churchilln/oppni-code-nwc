function [Xreg,stat] = gsreg_OP1( datamat, ParamCell )

modelstr = ParamCell{1};
stat = [];

if strcmpi(modelstr,'GAV')
    Xreg = mean(zscore(datamat'),2);
elseif strcmpi(modelstr,'PC1')
    [~,~,v] = svd( zscore(datamat')','econ');
    Xreg = v(:,1);
end

stat = [1];