function [Xreg,stat] = roireg_OP1( functional_run, roi_paths, ParamCell )
%
% ParamCell = {'WM'/'WM+CSF+tSD'}
%
% now revised -- simplest model averaging over tissue voxedls / no subregions or anything...

modelstr = ParamCell{1}; 

maskpath = roi_paths{1}; % expect : <maskpath>/func_<type>_mask_grp.nii
parcpath = roi_paths{2}; % expect : <parcpath>/Ugrp_<type>.nii

stat = [];

modelstr  = regexp(modelstr,'+','split');

Xreg=[];

str={'WM','CSF','tSD'};
for i=1:numel(modelstr)
    match=0;
    for j=1:numel(str)
        if contains(modelstr{i},[str{j}]) %% we have a match on type!
            match=match+1;
            tissmask = [maskpath,'/func_',str{j},'_mask_grp.nii'];

            unix('rm __opptmp_roireg_volblur.nii');

            unix(sprintf('3dBlurInMask -prefix __opptmp_roireg_volblur.nii -fwhm 6 -input %s -mask %s',functional_run,tissmask));
            V=load_untouch_niiz('__opptmp_roireg_volblur.nii');
            M=load_untouch_niiz(tissmask);
            volmat = nifti_to_mat(V,M);
            v = mean(zscore(volmat'),2);
            Xreg = [Xreg, zscore(v)];
            
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
