function [Xreg,stat] = roireg_OP2( functional_run, roi_paths, modelstr, decompstr )
%
% updated version --> does either IND or GRP pca of regions

maskpath = roi_paths{1}; % expect : <maskpath>/func_<type>_mask_grp.nii
parcpath = roi_paths{2}; % expect : <parcpath>/Ugrp_<type>.nii

stat = [];

modelstr  = regexp(modelstr,'+','split');
decompstr = regexp(decompstr,'+','split');

Xreg=[];

str={'WM','CSF','tSD'};
for i=1:numel(modelstr)
    match=0;
    for j=1:numel(str)
        if contains(modelstr{i},[str{j}]) %% we have a match on type!
            match=match+1;
            numpc=str2num( modelstr{i}((numel(str{j})+1):end) );
            tissmask = [maskpath,'/func_',str{j},'_mask_grp.nii'];
            weitmask = [parcpath,'/Ugrp_',str{j},'.nii'];

            unix('rm __opptmp_roireg_volblur.nii');

            unix(sprintf('3dBlurInMask -prefix __opptmp_roireg_volblur.nii -fwhm 6 -input %s -mask %s',functional_run,tissmask));
            V=load_untouch_niiz('__opptmp_roireg_volblur.nii');
            M=load_untouch_niiz(tissmask);
            volmat = nifti_to_mat(V,M);
            if strcmpi(decompstr{1},'IND') && strcmpi(decompstr{2},'PCA')
                [~,~,v] = svd( zscore(volmat')','econ');
            elseif strcmpi(decompstr{1},'GRP') && strcmpi(decompstr{2},'PCA')
                W = load_untouch_nii(weitmask);
                wtmat = nifti_to_mat(W,M);
                v = zscore(volmat')*wtmat;
                v = bsxfun(@times,v,sqrt(sum(v.^2)));
            else
                error('unrecognized roireg decomposition style')
            end
            Xreg = [Xreg, zscore(v(:,1:numpc))];
            
            unix('rm __opptmp_roireg_volblur.nii');
        end
    end
    if match==0
        error('ROIREG failed - cannot identify tissue type');
    elseif match>=2
        error('ROIREG failed - somehow matched to multiple tissue types?')
    end
end

rr = rank(Xreg);
[u,l,~]=svd(Xreg);
Xreg = u(:,1:rr);

stat = [rr];