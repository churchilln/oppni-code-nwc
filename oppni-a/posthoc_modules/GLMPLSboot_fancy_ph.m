function out = GLMPLSboot_fancy_ph( datamat, design, yscal, xscal, contrmat )
% allows you to take in multiple "datamat" blocks, regress and combine
NBlok = size(datamat,3);

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
NBOOT = 2000;

    disp('GLM, bootstrapped...');

    % augment design mat plus scale vec with intercept
    design = [ones(size(design,1),1), design];
    xscal = [0; xscal];
 
    % parameters
    n    = size(datamat,2);
    disp('running resampling...');
    bsrmat = zeros( size(datamat,1), size(design,2), NBOOT );

    for(bsr=1:NBOOT)
        list = ceil(n*rand(n,1));
        % boot sample
        D = datamat(:,list);
        y = design(list,:);
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
        % run ols regression (sans intercept)
        bsrmat(:,:,bsr) = D * (y / (y'*y));
    end
    bsrmat=bsrmat(:,2:end,:); % discarding intercept
 
    if isempty(contrmat)
        out.bsr     = mean(bsrmat,3)./std(bsrmat,0,3);
        out.bsr_p   = 2*min(cat(3,sum(bsrmat>0,3), sum(bsrmat<0,3)),[],3)./bsr;
        out.contr   = [];
    else
        for c = 1:size(contrmat,1)
            b = bsrmat(:,contrmat{c,2},:);
            if contrmat{c,1}==0
                a = 0;
            else
                a = bsrmat(:,contrmat{c,1},:);
            end
            bcontr = mean(b,2)-mean(a,2);
            out.bsr(:,c)   = mean(bcontr,3)./std(bcontr,0,3);
            out.bsr_p(:,c) = 2*min(cat(3,sum(bcontr>0,3), sum(bcontr<0,3)),[],3)./bsr;
        end
        out.contr = contrmat;
    end
    
    out.testname = 'glm_fancy_bootstrap';


% NOTES:
%
%
% - variance normalization of the outcome (y) seems to stabilize bootstrap
%   estimates and enhance effect sizes; also allows comparison across
%   voxels in terms of relative effect. Only needs deactivation if you
%   need absolute effects
% - centering on y only matters if you don't include an intercept!
%   not essential if it is a ddefault part of the model!
% 
% - variance normalization + centering of predictors (x) ...
%   for [continuous], centering and scaling tend to improv identifiability
%   and matrix condition. Also gives standardized coefficient values.
%   discards information about absolute effects though
%   centering > only affects how intercept is modelled(!); uncentered tries
%   to find average of "zeros" in mode / centered = average of all
%
% - for [categorical], centering similarly has a minimal effect, as long as
% there is an intercept in the model. Scaling also doesnt matter when only
% 1 vector BUT !!!!!....
%
% scaling can be problematic when multiple conditions to contrast against, with varying prevalence;
% lower prevalence conditions get smaller weights due to smaller vector norms
% may underestimate coeffs of effect, leading to trivial differences IF you are doing second-order contrasts 
% but also injects some variability in coeff estimates(!)
%
% In sum: models should by default (a) center & rescale Y
%                (b) include intercept term, and (c) center & rescale continuous
%                variables, and (d) center & rescale binary variables ONLY
%                MAYBE if a single predictor +interaction model, otherwise leave unchanged
