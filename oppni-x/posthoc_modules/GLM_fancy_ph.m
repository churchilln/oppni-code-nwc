function out = GLM_fancy_ph( datamat, design, yscal, xscal, contrmat )

% 1=center, 2=unitnorm, 3=both
if nargin<3 || isempty(yscal)
    yscal = 3;
end
if nargin<4 || isempty(xscal)
    for k=1:size(design,2)
        if numel(unique(design(:,k)))>2
            xscal(k,1) = 3;
        else
            xscal(k,1) = 0; % binary is unmodified
        end
    end
elseif numel(xscal)==1 && size(design,2)>1
    xscal = repmat( xscal, size(design,2), 1 );
end

if nargin<5
    contrmat = [];
end

econd = eps;

%%
disp('GLM, t-statistic...');

% augment design mat plus scale vec with intercept
design = [ones(size(design,1),1), design];
xscal = [0; xscal];

% parameters
n    = size(datamat,2);
k    = size(design, 2);

    % sample
    D = datamat;
    y = design;
    % standardized - zero-mean, unit var
    if yscal==1 || yscal==3
        D = bsxfun(@minus,  D,mean(D,2));
    end
    if yscal==2 || yscal==3
        D = bsxfun(@rdivide,D+econd,sqrt(sum(D.^2,2))+econd);
    end
    if sum(xscal==1 | xscal==3)>0
        yav = mean(y);
        yav( xscal~=1 & xscal~=3 ) = 0;
        y = bsxfun(@minus,  y,yav); 
    end
    if sum(xscal==2 | xscal==3)>0
        ynm = sqrt(sum(y.^2));
        ynm( xscal~=2 & xscal~=3 ) = 1;
        y = bsxfun(@rdivide,y+econd,ynm+econd);
    end
    % run ols regression--
    Beta = D * (y / (y'*y));
    D_estim = Beta * y';

    % now, estimate t-statistics on signal
    residvar    = var(D - D_estim,0,2);
    % catch instances of zero-variance
    residvar(var(D,0,2)==0) = econd;
    % --
    Xinvdiag    = diag(inv(y'*y));
    Tmap = Beta ./ sqrt( residvar * Xinvdiag(:)' );

    %---

    
    if isempty(contrmat)
        out.tstat    = Tmap(:,2:end); % discard intercept
        out.tstat_p  = 2.*tcdf( -abs(out.tstat), n-k-1 ); %2-tailed likelihood
        out.contr   = [];
    elseif iscell(contrmat) % cell array, each row is contrast; for cth, subtract contrmat{c,2}-contrmat{c,1} 
        for c = 1:size(contrmat,1)
            co = zeros(size(y,2),1);
            co(contrmat{c,1}+1)=-1; % adjust nottion to include intercept
            co(contrmat{c,2}+1)= 1;
            bcontr = Beta * co;
            Tmap(:,c) = (Beta * co) ./ sqrt( residvar * co'*((y'*y)\co) );
        end
        out.tstat    = Tmap; % discard intercept
        out.tstat_p  = 2.*tcdf( -abs(out.tstat), n-k-1 ); %2-tailed likelihood
        out.contr = contrmat;
    else % matrix, each row givdes the contrasts
        for c = 1:size(contrmat,1)
            co = [0 contrmat(c,:)]'; % pad to include intercept
            bcontr = Beta * co;
            Tmap(:,c) = (Beta * co) ./ sqrt( residvar * co'*((y'*y)\co) );
        end
        out.tstat    = Tmap; % discard intercept
        out.tstat_p  = 2.*tcdf( -abs(out.tstat), n-k-1 ); %2-tailed likelihood
        out.contr = contrmat;
    end


    % parameters - discard intercept!

    out.tstat    = Tmap(:,2:end);
    out.tstat_p  = 2.*tcdf( -abs(out.tstat), n-k-1 ); %2-tailed likelihood 
    
    out.testname = 'glm_tstat';
