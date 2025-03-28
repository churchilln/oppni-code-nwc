function [Xreg,stat] = roireg_OP3( functional_run, roi_paths, ParamCell )
%
% .roireg_OP3:
% .regression of physio signals from WM/CSF/high-variance areas
% .uses previously specified masks
%
% *** deprecated --> new version places intermediate files in correct
% folder


% ParamCell = {'CSF_LV_bilat'/'CSF_LV_left+CSF_LV_right+WM_CnS','PCA'/'AVG'}
%
% updated --> does either average or PCA-1 on set of tissue subregions

modelstr  = ParamCell{1}; % list of ROIs
decompstr = ParamCell{2}; % decomposition method

maskpath = roi_paths{1}; % expect : <maskpath>/func_<type>_mask_grp.nii
parcpath = roi_paths{2}; % expect : <parcpath>/roimask_resam_<type>_<subseg>.nii

stat = [];

modelstr  = regexp(modelstr,'+','split');
decompstr = regexp(decompstr,'-','split');

Xreg=[];

str={'WM','CSF','tSD'};
for i=1:numel(modelstr)
    match=0;
    for j=1:numel(str)
        if contains(modelstr{i},[str{j}]) %% we have a match on type!
            match=match+1;
            tissmask = [maskpath,'/func_',str{j},'_mask_grp.nii'];
            parcmask = [parcpath,'/roimask_resam_',modelstr{i},'.nii'];
            intrmask = '__opptmp_roireg_intersectmsk.nii';
            
            unix('rm __opptmp_roireg_intersectmsk.nii __opptmp_roireg_volblur.nii');
            V1=load_untouch_niiz(tissmask);
            V2=load_untouch_niiz(parcmask);
            V1.img = double(V1.img).*double(V2.img);
            save_untouch_niiz(V1,intrmask); clear V1 V2;

            unix(sprintf('3dBlurInMask -prefix __opptmp_roireg_volblur.nii -fwhm 6 -input %s -mask %s',functional_run,intrmask));
            V=load_untouch_niiz('__opptmp_roireg_volblur.nii');
            M=load_untouch_niiz(intrmask);
            volmat = nifti_to_mat(V,M);
            if strcmpi(decompstr{1},'AVG')
                v = mean(zscore(volmat'),2);
            elseif strcmpi(decompstr{1},'PCA')
                [~,~,v] = svd( zscore(volmat')','econ');
                v=v(:,1);
            elseif strcmpi(decompstr{1},'NUMPC')
                np = str2num(decompstr{2});
                [~,~,v] = svd( zscore(volmat')','econ');
                v=v(:,1:np);
            else
                error('unrecognized roireg decomposition style')
            end
            Xreg = [Xreg, zscore(v)];
            
            unix('rm __opptmp_roireg_intersectmsk.nii __opptmp_roireg_volblur.nii');
        end
    end
    if match==0
        error('ROIREG failed - cannot identify tissue type');
    elseif match>=2
        error('ROIREG failed - somehow matched to multiple tissue types?')
    end
end

rr = rank(Xreg);

if ~strcmpi(decompstr{1},'NUMPC')
    [u,l,~]=svd(Xreg);
    Xreg = u(:,1:rr);
else
    % for NUMPC, we retain unrotated components, so we can reconstruct ROI-specific contributions
    % .no modification to xreg
end

stat = [rr];
