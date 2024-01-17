function xreg = mpe_regressor_builder( mpefile, typestr )

if strcmpi(typestr,'0') || strcmpi(typestr,'OFF')
    xreg=[];
else
    if strcmpi(typestr,'1') || strcmpi(typestr,'ON')
        typestr = {'L','PC2'};
    else
        typestr = regexp(typestr,'+','split');
    end
    
    % load the file
    X = load(mpefile,'-ascii');
    % mean-centering for quadratic case (doesnt matter elsewhere really)
    X = bsxfun(@minus,X,mean(X,1));

    xreg=[];
    
    do_pca=0;
    for i=1:numel(typestr)
        
        if strcmpi(typestr{i},'L')
            xreg = [xreg zscore(X)];
        elseif strcmpi(typestr{i},'Q')
            xreg = [xreg zscore(X.^2)];
        elseif strcmpi(typestr{i},'LD') || strcmpi(typestr{i},'DL')
            xreg = [xreg zscore( [zeros(1,size(X,2)); X(2:end,:)-X(1:end-1,:)] )];
        elseif strcmpi(typestr{i},'QD')
           xtmp = X.^2;
           xreg = [xreg zscore( [zeros(1,size(xtmp,2)); xtmp(2:end,:)-xtmp(1:end-1,:)] )];
        elseif strcmpi(typestr{i},'DQ')
            xtmp = [zeros(1,size(X,2)); X(2:end,:)-X(1:end-1,:)];
            xreg = [xreg zscore(xtmp.^2)];
        elseif contains(typestr{i},'PC')
            do_pca=1;
            ixa=strfind(typestr{i},'PC');
            numpc=str2num(typestr{i}(ixa+2:end));
        else
            error('unrecognized motreg arg %s',typestr{i})
        end
    end
    
    if do_pca>0
        [u,l,v]=svd(xreg);
        xreg = zscore(u(:,1:numpc));
    end
end