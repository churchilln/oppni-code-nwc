function out = GLM_ph( datamat, design )

econd = eps;

%%
disp('GLM, t-statistic...');

% augment design mat plus scale vec with intercept
design = [ones(size(design,1),1), design];

% parameters
n    = size(datamat,2);
k    = size(design, 2);

    % sample
    D = datamat;
    y = design;
    % run ols regression (sans intercept)
    Beta = D * (y / (y'*y));

    % now, estimate t-statistics on signal
    residvar    = var(D - D_estim,0,2);
    % catch instances of zero-variance
    residvar(var(D,0,2)==0) = econd;
    % --
    Xinvdiag    = diag(inv(y'*y));
    Tmap = Beta ./ sqrt( residvar * Xinvdiag(:)' );


    % parameters
    out.tstat    = Tmap;
    out.tstat_p  = 2.*tcdf( -abs(out.tstat), n-k-1 ); %2-tailed likelihood 
    
    out.testname = 'glm_tstat';