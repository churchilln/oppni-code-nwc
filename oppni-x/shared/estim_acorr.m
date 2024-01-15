function out = estim_acorr( vol, mask, voxsz, maxwid )

% -- recommended maxwid=30

% predeclared
vol( mask==0 ) = NaN;
% get permissible range
xtmp = squeeze( sum(sum(mask,2),3) ); maxbnd(1) = voxsz(1)*( sum(xtmp/max(xtmp)>0.5) );
xtmp = squeeze( sum(sum(mask,1),3) ); maxbnd(2) = voxsz(2)*( sum(xtmp/max(xtmp)>0.5) );
xtmp = squeeze( sum(sum(mask,1),2) ); maxbnd(3) = voxsz(3)*( sum(xtmp/max(xtmp)>0.5) );
maxbnd = min([maxwid maxbnd]);

%% dim-x
    maxlag = floor( maxbnd / voxsz(1) );
    c = 1;
    w = NaN;
    for(lag=1:maxlag)
        a = reshape( vol(1:(end-lag),:,:),[],1);
        b = reshape( vol((1+lag):end,:,:),[],1);
        ix = isfinite(a) & isfinite(b);
        if(sum(ix)>3)
            c(lag+1,1) = corr(a(ix),b(ix));
            w(lag+1,1) = sum(ix);
        else
            c(lag+1,1) = NaN; 
            w(lag+1,1) = sum(ix);
        end
    end
    % autocorrelation decay
    out.ACD(1) = (c(3)-2*c(2)+c(1))./voxsz(1);
    % curve fitting, weighted monoexponential decay
    x = voxsz(1) * (1:(maxlag-1))';
    y = c(2:maxlag);
    q = w(2:maxlag);
    %
    des = x.*q;
    tar = log(y).*q;
    out.TAU(1) = sum(des.*tar) / sum(des.*des);

%% dim-y
    maxlag = floor( maxbnd / voxsz(2) );
    c = 1;
    w = NaN;
    for(lag=1:maxlag)
        a = reshape( vol(:,1:(end-lag),:),[],1);
        b = reshape( vol(:,(1+lag):end,:),[],1);
        ix = isfinite(a) & isfinite(b);
        if(sum(ix)>3)
            c(lag+1,1) = corr(a(ix),b(ix));
            w(lag+1,1) = sum(ix);
        else
            c(lag+1,1) = NaN; 
            w(lag+1,1) = sum(ix);
        end
    end
    % autocorrelation decay
    out.ACD(2) = (c(3)-2*c(2)+c(1))./voxsz(2);
    % curve fitting, weighted monoexponential decay
    x = voxsz(2) * (1:(maxlag-1))';
    y = c(2:maxlag);
    q = w(2:maxlag);
    %
    des = x.*q;
    tar = log(y).*q;
    out.TAU(2) = sum(des.*tar) / sum(des.*des);

%% dim-z
    maxlag = floor( maxbnd / voxsz(3) );
    c = 1;
    w = NaN;
    for(lag=1:maxlag)
        a = reshape( vol(:,:,1:(end-lag)),[],1);
        b = reshape( vol(:,:,(1+lag):end),[],1);
        ix = isfinite(a) & isfinite(b);
        if(sum(ix)>3)
            c(lag+1,1) = corr(a(ix),b(ix));
            w(lag+1,1) = sum(ix);
        else
            c(lag+1,1) = NaN; 
            w(lag+1,1) = sum(ix);
        end
    end
    % autocorrelation decay
    out.ACD(3) = (c(3)-2*c(2)+c(1))./voxsz(3);
    % curve fitting, weighted monoexponential decay
    x = voxsz(3) * (1:(maxlag-1))';
    y = c(2:maxlag);
    q = w(2:maxlag);
    %
    des = x.*q;
    tar = log(y).*q;
    out.TAU(3) = sum(des.*tar) / sum(des.*des);  
    
    % averaged
    out.ACD(4) = mean(out.ACD(1:3));
    out.TAU(4) = mean(out.TAU(1:3));

    out.ACD = -out.ACD;
    out.TAU = -out.TAU;

    %% currently unused parametric models...
%         x = voxsz(1) * (1:(maxlag-1))';
%     y = c(2:maxlag);
%     q = w(2:maxlag);
%     % curve fitting, weighted monoexponential decay
%     des = x.*q;
%     tar = log(y).*q;
%     out.TAU(1) = (des'*des) \ (des'*tar);
%     % curve fitting, quadratic curvature
%     des = [x x.^2] .* [q q];
%     tar = (y-1) .* q;
%     Beto = (des'*des) \ (des'*tar);
%     xt = 1:maxbnd; k = 2*Beto(2) ./ (1 + (2.*Beto(2).*xt + Beto(1)).^2).^(3/2);
%     out.QMC(1) = max(k);
