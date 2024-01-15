function [Xreg,stat] = roireg_OP1( functional_run, wcv_masks, wcv_weights, modelstr, decompstr )
%
% original version --> does either IND or GRP pca of regions


stat = [];

modelstr  = regexp(modelstr,'+','split');
decompstr = regexp(decompstr,'-','split');

Xreg=[];

str={'W','C','V'};
for i=1:numel(modelstr)

    for j=1:numel(str)
        if strcmpi(modelstr{i}(1),[str{j}])

            numpc=str2num(modelstr{i}(2:end));

            unix('rm __opptmp_roireg_volblur.nii');

            unix(sprintf('3dBlurInMask -prefix __opptmp_roireg_volblur.nii -fwhm 6 -input %s -mask %s',functional_run,wcv_masks{j}));
            V=load_untouch_niiz('__opptmp_roireg_volblur.nii');
            M=load_untouch_niiz(wcv_masks{j});
            volmat = nifti_to_mat(V,M);
            if strcmpi(decompstr{1},'IND') && strcmpi(decompstr{2},'PCA')
                [~,~,v] = svd( zscore(volmat')','econ');
            elseif strcmpi(decompstr{1},'GRP') && strcmpi(decompstr{2},'PCA')
                W = load_untouch_nii(wcv_weights{j});
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
end

rr = rank(Xreg);
[u,l,~]=svd(Xreg);
Xreg = u(:,1:rr);

stat = [rr];
