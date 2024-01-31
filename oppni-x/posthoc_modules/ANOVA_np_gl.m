function out = ANOVA_np_gl( X0 )
%
% non-parametric anova ... order (V x S x T)
%

iter_perm = 2000;
iter_boot = 2000;
[nv,ns,nt]=size(X0);

RAT = modelfit_TR( X0 );

% -----------------longitudinal model, perm-----------------

pset=0; 
Xprm=X0; RATprm_avg = 0;
for(iter=1:iter_perm)
    iter,
    for(s=1:ns)
        Xprm(:,s,:) = X0(:,s,randperm(4));
    end
    RATprm = modelfit_TR( Xprm );
    pset = pset + double( RATprm>RAT )./iter_perm;
    RATprm_avg = RATprm_avg + RATprm./iter_perm;
end

% -----------------longitudinal model, boot-----------------

pset2=0;
Xboot=X0; RATboot_avg = 0;
for(iter=1:iter_boot)
    iter,
    list = ceil( ns*rand(ns,1) );
    Xboot = X0(:,list,:);
    RATboot = modelfit_TR( Xboot );
    pset2 = pset2 + double( RATboot<RATprm_avg )./iter_boot;
    RATboot_avg = RATboot_avg + RATboot./iter_boot;
end

out.p_prm = pset;
out.p_boot = pset2;
out.Ratio = RATboot_avg;

%%
function RAT = modelfit_TR( X )

[nv,ns,nt]=size(X);

mut = mean(X,2); % mean per time
mus = mean(X,3); % mean per subject
mug = sum(sum(X,2),3)./(ns*nt); % global mean
ssT = ns*sum( bsxfun(@minus, mut,mug).^2 ,3); % ss-err (time)
ssR = sum( sum( bsxfun(@plus,bsxfun(@minus,bsxfun(@minus,X,mut),mus),mug).^2, 2),3); % ss-err (resid)
RAT = ssT./ssR; % ss-err ratio
