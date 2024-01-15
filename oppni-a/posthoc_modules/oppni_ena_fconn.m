function out = oppni_ena_fconn( Datacell, Mask, Template, Designmat, odir )

[NS,NR] = size(Datacell);

M=load_untouch_nii(Mask);
mask=double(M.img);

unix(sprintf('3dresample -master %s -input %s -prefix %s_template.nii',Mask,Template,odir));
T=load_untouch_nii(sprintf('%s_template.nii',odir))
templ = nifti_to_mat(T,M);

for s=1:NS
    for r=1:NR
        if r>1
            error('multirun unsupported for now');
        end

        V=load_untouch_nii(Datacell{s,r});
        volmat = nifti_to_mat( V,M ); nt=size(volmat,2);
        volmat = zscore(volmat')';
        tser   = zscore( volmat'*templ );
        fconn(:,s) = (volmat * tser)./(nt-1);
    end
end

zconn = 0.5*( log(1+fconn) - log(1-fconn));

out = GLM_ph( zconn, Designmat );

VOLMAT = zeros( [size(mask) size(out.tstat,2)]);
for p=1:size(out.tstat,2)
    tmp=mask;tmp(tmp>0) = out.tstat(:,p);
    VOLMAT(:,:,:,p) = tmp;
end
nii = make_nii(VOLMAT,M.hdr.dime.pixdim(2:4));
nii.hdr.hist = M.hdr.hist;
save_nii(nii,sprintf('%s_ttest.nii',odir));



%%
function out = GLM_ph( datamat, design )

    disp('GLM, t-statistic...');

    n    = size(datamat,2);
    k    = size(design, 2);
    op   = GLM_model_fmri( datamat, 0, [], design, [], [] );
    
    % parameters
    out.tstat    = op.Tmap_signl;
    out.tstat_p  = 2.*tcdf( -abs(out.tstat), n-k-1 ); %2-tailed likelihood
    
    out.testname = 'glm_tstat';
