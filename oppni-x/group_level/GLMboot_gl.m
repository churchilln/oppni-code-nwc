function out = GLMboot_gl( datamat, design )

econd = eps;
NBOOT = 1000;

    disp('GLM, bootstrapped...');
 
    % parameters
    n    = size(datamat,2);
    k    = size(design, 2);
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
        op   = GLM_model_fmri( D, 0, [],y, [] );
        bsrmat(:,:,bsr) = op.Beta_signl;
    end
    out.bsr     = mean(bsrmat,3)./std(bsrmat,0,3);
    out.bsr_p   = 2*min(cat(3,sum(bsrmat>0,3), sum(bsrmat<0,3)),[],3)./bsr;
   
    out.testname = 'glm_bootstrap';
