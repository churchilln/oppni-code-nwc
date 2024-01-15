function out = batch_outlier_testing( X, testtype, alpha, adjust )

[N,K] = size(X);

for k=1:K

    if strcmpi(testtype{k},'norm')
        pvalset(:,k) = 2*normcdf( -abs(zscore(X(:,k))) ); % 2-tailed

    elseif strcmpi(testtype{k},'norm-')
        pvalset(:,k) =   normcdf( zscore(X(:,k)) ); % 1-tailed (lower)

    elseif strcmpi(testtype{k},'norm+')
        pvalset(:,k) =   normcdf( -zscore(X(:,k)) ); % 1-tailed (upper)

    elseif strcmpi(testtype{k},'gam')
        xtmp    = X(:,k)./max(X(:,k));
        xtmp(xtmp<eps)=eps; % cannot be zeros
        par_ab  = gamfit( xtmp );
        pvalset(:,k) = 1-gamcdf( xtmp, par_ab(1), par_ab(2) ); % 1-tailed (upper)

    elseif strcmpi(testtype{k},'1-gam')
        xtmp    = 1 - (  X(:,k)./max(X(:,k))  );
        xtmp(xtmp<eps)=eps; % cannot be zeros
        par_ab  = gamfit( xtmp );
        pvalset(:,k) = 1-gamcdf( xtmp, par_ab(1), par_ab(2) ); % 1-tailed (upper)

    elseif strcmpi(testtype{k},'corr')
        xtmp = 0.5*( log(1+X(:,k)) - log(1-X(:,k)) );
        pvalset(:,k) = 2*normcdf( -abs(zscore(xtmp)) ); % 2-tailed

    else
        error('cannot identify appropriate stat test!')
    end
        
    % thresholding
    %
    if    ( strcmpi(adjust,'Uncorrected') )  
        tvalset(:,k) = pvalset(:,k) < alpha;
    elseif( strcmpi(adjust,'FDR') )       
        [~,tvalset(:,k)] = fdr( pvalset(:,k),'p',alpha,0);
    elseif( strcmpi(adjust,'Bonferroni') )   
        tvalset(:,k) = pvalset(:,k) < alpha/numel(pvalset(:,k));
    end

end

out.pvl = pvalset;
out.thr = tvalset;
