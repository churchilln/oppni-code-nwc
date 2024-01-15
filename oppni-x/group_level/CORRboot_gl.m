function out = CORRboot_gl( datamat, design )

econd = eps;
NBOOT = 1000;

    disp('CORR, bootstrapped...');
 
    % parameters
    n    = size(datamat,2);
    k    = size(design, 2);

    % whole-data statistic
    list = 1:n;
    D = datamat(:,list);
    D = bsxfun(@minus,  D,mean(D,2));
    D = bsxfun(@rdivide,D+econd,sqrt(sum(D.^2,2))+econd);
    y = design(list,:);
    y = bsxfun(@minus,  y,mean(y,1));        
    y = bsxfun(@rdivide,y+econd,sqrt(sum(y.^2))+econd);        
    out.corr = D*y;
    
    disp('running resampling...');
    bsrmat = zeros( size(datamat,1), k, NBOOT );
    for(bsr=1:NBOOT)
        [bsr NBOOT],
        list = ceil(n*rand(n,1));
        D = datamat(:,list);
        D = bsxfun(@minus,  D,mean(D,2));
        D = bsxfun(@rdivide,D+econd,sqrt(sum(D.^2,2))+econd);
        y = design(list,:);
        y = bsxfun(@minus,  y,mean(y,1));        
        y = bsxfun(@rdivide,y+econd,sqrt(sum(y.^2))+econd);        
        bsrmat(:,:,bsr) = D*y;
    end
    out.bsr     = mean(bsrmat,3)./std(bsrmat,0,3);
    out.bsr_p   = 2*min(cat(3,sum(bsrmat>0,3), sum(bsrmat<0,3)),[],3)./bsr;
   
    out.testname = 'corr_bootstrap';
