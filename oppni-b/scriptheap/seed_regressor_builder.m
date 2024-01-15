function xreg = seed_regressor_builder( functional_run, wcv_masks, wcv_weights, typestr )

if strcmpi(typestr,'0') || strcmpi(typestr,'OFF')
    xreg=[];
else
    if strcmpi(typestr,'1') || strcmpi(typestr,'ON')
        typestr = {'WPC2','CPC2','GRP'};
    else
        typestr = regexp(typestr,'+','split');
    end
    
    xreg=[];
    
    iflag=0;
    for i=1:numel(typestr)
        if strcmpi(typestr{i},'IND')
            iflag=1;
        elseif strcmpi(typestr{i},'GRP')
            iflag=-1;
        end
    end
    if iflag==0 
        error('roireg needs either IND or GRP argument');
    end
    
    str={'W','C','V'};
    for i=1:numel(typestr)
    
        for j=1:numel(str)
            if contains(typestr{i},[str{j},'PC'])
    
                ixa=strfind(typestr{i},[str{j},'PC']);
                numpc=str2num(typestr{i}(ixa+3:end));
                %
                unix('rm tmp_volblur.nii');
                unix(sprintf('3dBlurInMask -prefix tmp_volblur.nii -fwhm 6 -input %s -mask %s',functional_run,wcv_masks{j}));
                V=load_untouch_nii('tmp_volblur.nii');
                M=load_untouch_nii(wcv_masks{j});
                volmat = nifti_to_mat(V,M);
                if iflag==1
                    [~,~,v] = svd( zscore(volmat')','econ');
                elseif iflag==-1
                    W = load_untouch_nii(wcv_weights{j});
                    wtmat = nifti_to_mat(W,M);
                    v = zscore(volmat')*wtmat;
                    v = bsxfun(@times,v,sqrt(sum(v.^2)));
                end
                xreg = [xreg, v(:,1:numpc)];
                unix('rm tmp_volblur.nii');
            end
        end
    end

end
