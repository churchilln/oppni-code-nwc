function out = Permute_gl( datamat, design, estim )
%
% . Bootstrapped analysis, compares 1-group rel. 0, or between 2 groups
%

NBOOT = 2000;

if(nargin<3) estim='mean'; end

if    ( isempty(design) || numel(unique(design))==1 )
    
    disp('Permutation, 1-sample...');
    
    % parameters
    n    = size(datamat,2);
    disp('running resampling...');
    parest = mean(datamat,2); % observed parameter
    prmmat = zeros( size(datamat,1), NBOOT );
    for(bsr=1:NBOOT)
        [bsr NBOOT],
        signr = sign( rand(n,1)-0.5 );
        switch estim % select bootstrapping statistic
            case 'mean'
                prmmat(:,bsr) = mean(datamat * diag(signr),2);
            case 'stdev'
                prmmat(:,bsr) = std(datamat * diag(signr),0,2);
            case 'robmean'
                prmmat(:,bsr) = location_mest( datamat * diag(signr), 1.5, 50 );
            case 'robstdev'
                prmmat(:,bsr) = 1.49*median( abs(bsxfun(@minus,datamat * diag(signr),median(datamat * diag(signr),2))),2);
            otherwise
                error(['undefined parameter input:', estim]);
        end        
    end
    out.perm   = parest./std(prmmat,0,2);
    prmmat     = bsxfun(@minus,prmmat,parest);
    out.perm_p = 2*min([mean(prmmat>0,2), mean(prmmat<0,2)],[],2);
    
    out.testname = 'permutation_1samp';

elseif( numel(unique(design))>1 && size(design,1)>1 && size(design,2)==1 )
    
    ix = unique(design);
    datamat1 = datamat(:,design==ix(1));
    datamat2 = datamat(:,design==ix(2));
    disp('Permutation, 2-sample (unpaired)...');

    % parameters
    n1 = size(datamat1,2);
    n2 = size(datamat2,2);
    disp('running resampling...');
    parest = mean(datamat2,2)-mean(datamat1,2); % observed parameter
    prmmat = zeros( size(datamat,1), NBOOT );
    for(bsr=1:NBOOT)
        [bsr NBOOT],
        list = randperm(n1+n2);
        prmmat(:,bsr) = mean( datamat(:,list(1:n2)),2 ) - mean( datamat(:,list(n2+1:n2+n1)),2 );
        
        switch estim % select bootstrapping statistic
            case 'mean'
                prmmat(:,bsr) = mean( datamat(:,list(1:n2)),2 ) - mean( datamat(:,list(n2+1:n2+n1)),2 );
            case 'stdev'
                prmmat(:,bsr) = std( datamat(:,list(1:n2)),0,2 ) - std( datamat(:,list(n2+1:n2+n1)),0,2 );
            case 'robmean'
                prmmat(:,bsr) = location_mest( datamat(:,list(1:n2)), 1.5, 50 ) - location_mest( datamat(:,list(n2+1:n2+n1)), 1.5, 50 );
            case 'robstdev'
                a = 1.49* ( median( abs(bsxfun(@minus,datamat(:,list(1:n2)),median(datamat(:,list(1:n2)),2))),2) );
                b = 1.49* ( median( abs(bsxfun(@minus,datamat(:,list(n2+1:n2+n1)),median(datamat(:,list(n2+1:n2+n1)),2))),2) );
                prmmat(:,bsr) = a-b;
            otherwise
                error(['undefined parameter input:', estim]);
        end        
        
    end
    out.perm   = parest./std(prmmat,0,2);
    prmmat     = bsxfun(@minus,prmmat,parest);
    out.perm_p = 2*min([mean(prmmat>0,2), mean(prmmat<0,2)],[],2);
    
    out.testname = 'permutation_2samp_unpair';
    
elseif( numel(unique(design))>1 && size(design,1)>1 && size(design,2)==2 )
    
    ix = unique(design(:,1));
    datamat1 = datamat(:,design(:,1)==ix(1)); 
    datamat2 = datamat(:,design(:,1)==ix(2)); 
    
    clear datamat; %% wipe out old datamat
   
    des1     = design(design(:,1)==ix(1),2);
    des2     = design(design(:,1)==ix(2),2);
    
    if( size(datamat1,2) ~= size(datamat2,2) ) 
        error('unequal splits, cannot match');
    else
        datamat = zeros( size(datamat1) );
    end
    
    tmp=zeros(size(datamat2));
    for(i=1:length(des1))
        tmp(:,i) = datamat2(:,des2==des1(i)); 
    end
    datamat2 = tmp; clear tmp; % reorder to match datamat1
 
    disp('Permutation, 2-sample (paired)...');

    % parameters
    n    = size(datamat,2);
    disp('running resampling...');
    parest = mean(datamat,2); % observed parameter
    prmmat = zeros( size(datamat,1), NBOOT );
    for(bsr=1:NBOOT)
        [bsr NBOOT],
        signr = sign( rand(n,1)-0.5 );
        prmmat(:,bsr) = mean( datamat * diag(signr) ,2);
        
        signr = sign( rand(n,1)-0.5 );
        dataprm1 = datamat1 * diag(signr) + datamat2 * diag((1-signr));
        dataprm2 = datamat1 * diag((1-signr)) + datamat2 * diag(signr);
        
        switch estim % select bootstrapping statistic
            case 'mean'
                prmmat(:,bsr) = mean(dataprm2,2) - mean(dataprm1,2);
            case 'stdev'
                prmmat(:,bsr) = std(dataprm2,0,2) - std(dataprm1,0,2);
            case 'robmean'
                prmmat(:,bsr) = location_mest( dataprm2, 1.5, 50 ) - location_mest( dataprm1, 1.5, 50 );
            case 'robstdev'
                prmmat(:,bsr) = 1.49* ( median( abs(bsxfun(@minus,dataprm2,median(dataprm2,2))),2)  -  median( abs(bsxfun(@minus,dataprm1,median(dataprm1,2))),2) );
            otherwise
                error(['undefined parameter input:', estim]);
        end        
        
    end
    out.perm   = parest./std(prmmat,0,2);
    prmmat     = bsxfun(@minus,prmmat,parest);
    out.perm_p = 2*min([mean(prmmat>0,2), mean(prmmat<0,2)],[],2);
    
    out.testname = 'permutation_2samp_pair';
    
end

