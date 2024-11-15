function [Xreg,stat] = motreg_OP1( mpefile, ParamCell )

modelstr = ParamCell{1}; 
decompstr = ParamCell{2};

stat = [];

modelstr  = regexp(modelstr,'+','split');
decompstr = regexp(decompstr,'-','split');

% load the file
X = load(mpefile,'-ascii');
% mean-centering for quadratic case (doesnt matter elsewhere really)
X = bsxfun(@minus,X,mean(X,1));

Xreg=[];

for i=1:numel(modelstr)
    
    if strcmpi(modelstr{i},'L')
        Xreg = [zscore(X)];
    elseif strcmpi(modelstr{i},'LQ')
        Xreg = [zscore(X) zscore(X.^2)];
    elseif strcmpi(modelstr{i},'LD')
        Xreg = [zscore(X) zscore( [zeros(1,size(X,2)); X(2:end,:)-X(1:end-1,:)] )];
    elseif strcmpi(modelstr{i},'LQD') % linear and quad., +derivative of both
        xtmp = X.^2;
        Xreg = [zscore(X) zscore(xtmp), zscore([zeros(1,size(X,2)); X(2:end,:)-X(1:end-1,:)]) zscore( [zeros(1,size(xtmp,2)); xtmp(2:end,:)-xtmp(1:end-1,:)] )];
    elseif strcmpi(modelstr{i},'LDQ') % linear and deriv, +quadratic of both
        xtmp = [zeros(1,size(X,2)); X(2:end,:)-X(1:end-1,:)]; % deriv.
        Xreg = [zscore(X) zscore(xtmp), zscore(X.^2) zscore(xtmp.^2)];
    else
        error('unrecognized motreg arg %s',modelstr{i})
    end
end

rr = rank(Xreg);
[u,l,~]=svd(Xreg);

stat = [rr];

if strcmpi(decompstr{1},'FULL')
    disp('motreg, undecomped!')
    Xreg = u(:,1:rr);
elseif strcmpi(decompstr{1},'NUMPC')
    np = str2num(decompstr{2});
    Xreg = u(:,1:np);
elseif strcmpi(decompstr{1},'PCTPC')
    nv = str2num(decompstr{2});
    np = sum( cumsum( diag(l.^2)./trace(l.^2) ) < nv ) +1;
    Xreg = u(:,1:np);
else
    error('unrecognized decomp args %s and %s\n',decompstr{1},decompstr{2});
end
