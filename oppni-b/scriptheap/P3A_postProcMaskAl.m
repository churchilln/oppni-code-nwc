function P3_postProcessingQC( inputfile, pipelinefile, paramlist, outpath )
%
% . this script:
% . runs da qc
%

% initializing structure and running checkes
[subject_list, InputStruct_aug, ParamStruct_aug] = P0_fmri_populateDirectories( inputfile, pipelinefile, paramlist, outpath );
% now augmenting outpath... do this after P0!
outpath = fullfile(outpath,'fmri_proc'); % subdir should be fmri_proc
% check for missing files from previous step too
% File_Existence_Checker(InputStruct_aug,outpath,2);  %%-----> NEEDS ANOTHER STEPPPPPPPPPPPP

% warp-specific path
catpath = [ParamStruct_aug.WarpPrefix{1},'.',ParamStruct_aug.VOXRES,'mm'];

% func
Mref_func = load_untouch_nii([outpath,'/_group_level/masks/',catpath,'/func_brain_mask_grp.nii']);
V    = load_untouch_nii([outpath,'/_group_level/brain_maps/',catpath,'/func_tAV_grp.nii']);
avec = nifti_to_mat(V,Mref_func);
V    = load_untouch_nii([outpath,'/_group_level/brain_maps/',catpath,'/func_tSD_grp.nii']);
svec = nifti_to_mat(V,Mref_func);

% anat
Mref_anat = load_untouch_nii([outpath,'/_group_level/masks/',catpath,'/anat_brain_mask_grp.nii']);
V    = load_untouch_nii([outpath,'/_group_level/brain_maps/',catpath,'/anat_brain_grp.nii']);
bvec = nifti_to_mat(V,Mref_anat);
V    = load_untouch_nii([outpath,'/_group_level/brain_maps/',catpath,'/anat_CSF_grp.nii']);
tmat(:,1) = nifti_to_mat(V,Mref_anat);
V    = load_untouch_nii([outpath,'/_group_level/brain_maps/',catpath,'/anat_GM_grp.nii']);
tmat(:,2) = nifti_to_mat(V,Mref_anat);
V    = load_untouch_nii([outpath,'/_group_level/brain_maps/',catpath,'/anat_WM_grp.nii']);
tmat(:,3) = nifti_to_mat(V,Mref_anat);

for ns=1:7%numel(subject_list)

    % check existence of subject specific struct file
    if exist(fullfile( outpath,subject_list{ns},'InputStruct_ssa.mat'),'file')
        load( fullfile( outpath,subject_list{ns},'InputStruct_ssa.mat') )
    else
        error('cannot find Input struct file for subject: %s \n');
    end

    fprintf('\n===> qc collection. now on subj %u/%u: %s...\n',ns,numel(subject_list),subject_list{ns}),

    opath1 = fullfile(outpath,InputStruct_ssa.PREFIX,'func_proc_p1');
    opath3a= fullfile(outpath,InputStruct_ssa.PREFIX,'anat_proc',ParamStruct_aug.WarpPrefix{1});
    opath3f= fullfile(outpath,InputStruct_ssa.PREFIX,'func_proc_p1',[ParamStruct_aug.WarpPrefix{1},'.',ParamStruct_aug.VOXRES,'mm']);
    opath4f= fullfile(outpath,InputStruct_ssa.PREFIX,'func_proc_p2',[ParamStruct_aug.WarpPrefix{1},'.',ParamStruct_aug.VOXRES,'mm']);

    catpath = [ParamStruct_aug.WarpPrefix{1},'.',ParamStruct_aug.VOXRES,'mm'];

    QCStruct_quant(ns).PREFIX = InputStruct_ssa.PREFIX;
    QCStruct_quant(ns).N_anat = InputStruct_ssa.N_anat;
    QCStruct_quant(ns).N_func = InputStruct_ssa.N_func;

    %% ANATOMICAL

    nr=1; % anatomical metrics

    M=load_untouch_nii(sprintf('%s/seg/anat_mask.nii',opath3a));
    qcmat_a(ns,1) = sum(M.img(:));
    M=load_untouch_nii(sprintf('%s/anat_warped_mask.nii',opath3a));
    qcmat_a(ns,2) = sum(M.img(:));
    % jaccard overlap viz groupmask
    a=Mref_anat.img(:);
    b=M.img(:);
    qcmat_a(ns,3) = sum( a.*b ) ./ sum( double(a | b) );

    V = load_untouch_nii(sprintf('%s/anat_warped.nii',opath3a));
    qcmat_a(ns,4) = corr( bvec, nifti_to_mat(V,Mref_anat));
    % signal scaling
    qcmat_a(ns,5) = mean( nifti_to_mat(V,Mref_anat) );
    qcmat_a(ns,6) =  std( nifti_to_mat(V,Mref_anat) );
    
    V = load_untouch_nii(sprintf('%s/seg/anat_seg_CSF_warped.nii',opath3a));
    qcmat_a(ns,7) = corr( tmat(:,1), nifti_to_mat(V,Mref_anat));
    qcmat_a(ns,8) = mean( nifti_to_mat(V,Mref_anat) );
    V = load_untouch_nii(sprintf('%s/seg/anat_seg_GM_warped.nii',opath3a));
    qcmat_a(ns,9) = corr(  tmat(:,2), nifti_to_mat(V,Mref_anat));
    qcmat_a(ns,10) = mean( nifti_to_mat(V,Mref_anat) );
    V = load_untouch_nii(sprintf('%s/seg/anat_seg_WM_warped.nii',opath3a));
    qcmat_a(ns,11) = corr(  tmat(:,3), nifti_to_mat(V,Mref_anat));
    qcmat_a(ns,12) = mean( nifti_to_mat(V,Mref_anat) );


	V0= load_untouch_nii([opath3a,'/anatUAC.x.nii']);
    imag0= double(V0.img);
    %
    V = load_untouch_nii([opath3a,'/anatSS.x.nii.nii']);
	mask = double(V.img>eps);
	imag = double(V.img);
	imvc = imag(mask>0);
    
%     %--- quick plotting
%     ixnz = find( squeeze( sum(sum(imag>0,1),2))>0);
%     ixnz = ixnz([1 end]);
%     slc = round(linspace(1,ixnz(2)-ixnz(1),20));
%     bnd = prctile( imag(:),[2.5 97.5]);
%     h0=figure;
%     set(gcf, 'Units', 'normalized');
%     set(gcf, 'Position', [0.05 0.15 0.90 0.75]);
%     %
%     for k=1:(numel(slc)-2)
%         A = flipud( imag0(:,:,slc(k+1)+ixnz(1))' );
%         B = flipud( imag( :,:,slc(k+1)+ixnz(1))' );
%         C = flipud( mask(:,:,slc(k+1)+ixnz(1))' );
%         
%         sfh = subplot(3,6,k   ); 
%         ha = imshow(A, [], 'Colormap', gray(256));
%         hold on;
%         cmap = jet(256);
%         rgbImage = ind2rgb(uint8(225 * C), cmap);
%         hb = imshow(rgbImage)
%         set(hb,'AlphaData',C.*0.3);
% 
%         sfh.Position = sfh.Position + [0 0 0.031 0.031];
%     end
%     text(10,-10,subject_list{ns})

    %--- trick plotting
    ixnz = find( squeeze( sum(sum(imag>0,1),2))>0);
    ixnz = ixnz([1 end]);
    slc = round(linspace(1,ixnz(2)-ixnz(1),20));
    bnd = prctile( imag(:),[2.5 97.5]);
    h0=figure;
    set(gcf, 'Units', 'normalized');
    set(gcf, 'Position', [0.05 0.15 0.90 0.75]);
    %
    for k=1:(numel(slc)-2)
        A = flipud( imag0(:,:,slc(k+1)+ixnz(1))' );
        B = flipud( imag( :,:,slc(k+1)+ixnz(1))' );
        
        sfh = subplot(3,6,k   ); 
        cmap = jet(256);
        rgbImage = ind2rgb(uint8(255 * mat2gray(B)), cmap);
        ha = imshow(A, [], 'Colormap', gray(256))
        hold on;
        hb = imshow(rgbImage)
        % Set opacity/transparency to something less than 1 (alpha).  
        % 1 is the default and it means the last image is opaque and the image below can't be seen.
        hb.AlphaData = 0.3;


        sfh.Position = sfh.Position + [0 0 0.031 0.031];
    end
    text(10,-10,subject_list{ns})


    
    return;
    %% FUNCTIONAL

    for nr=1:InputStruct_ssa.N_func % functional metrics
        
        % unwarped mask [mask quality / brain shape] - total brain volume
        M=load_untouch_nii(sprintf('%s/func%u_2std_mask.nii',opath1,nr));
        qcmat_f(ns,1) = sum(M.img(:));
        % warped mask [mask quality / brain shape] - total brain volume
        M=load_untouch_nii(sprintf('%s/func%u_warped_mask.nii',opath4f,nr));
        qcmat_f(ns,2) = sum(M.img(:));
        % jaccard overlap viz groupmask
        a=Mref_func.img(:);
        b=M.img(:);
        qcmat_f(ns,3) = sum( a.*b ) ./ sum( double(a | b) );

        V = load_untouch_nii(sprintf('%s/func%u_warped_tav.nii',opath4f,nr));
        qcmat_f(ns,4) = corr( avec, nifti_to_mat(V,Mref_func));
        V = load_untouch_nii(sprintf('%s/func%u_warped_tsd.nii',opath4f,nr));
        qcmat_f(ns,5) = corr( svec, nifti_to_mat(V,Mref_func));

    end

end

return;

% > metrics focused on end of processing chain and/or key steps to get there
%
% anatomical categories:
%
% *maskvol_unwarp, maskvol_warp, mask_ovl_grp --> quality/consistency of masking (can also catch poor snr or strange head sizes)
% *corr_brain --> quality of image warping
% corr_csf/gm/wm, frac_csf/gm/wm --> quality/consistency of tissue segmentation
% *brain_av,brain_sd -->  signal scaling/gain AFTER cleanup (eg site/scanner effects)
% 
% functional categories:
%
% maskvol_unwarp, maskvol_warp, mask_ovl_grp --> quality/consistency of masking (can also catch poor snr or strange head sizes)
% *corr_tav, corr_tsd --> quality of image warping
% mot_avg/max_tot/rel --> amount of head motion
% dvr_avg/max_tot/rel --> amount of remaining (potentially motion related) variance
% *tav_map_av,tav_map_sd,tsd_map_av,tsd_map_sd --> signal scaling/gain AFTER cleanup (eg site/scanner effects)
% tsn_av,tsn_sd --> signal-to-noise measure
% krt_peak, krt_avg, krt_frac --> density of outliers
% corr_gs_av, corr_gs_sd --> global signal effects
%
% > for mot and dvr, tried mean( x.^2 ) and mean( (x-mean(x)).^2 ) as
% indices of overall displacement; discarded in favour of "tot" which is v.
% similar (corr>0.95) for single run, and not influenced by mean offset
% effects in multirun data
%
% > for gs-corr, used combined metric ratio -- better distribution,
% relevant to how "uniformly high" global signal is
%
% >

disp('exporting quantitative qc files.')

% export the .mat structure to folder
save( fullfile(outpath,'_group_level','QC','qc.quant','QCStruct_quant.mat'), 'QCStruct_quant' );

% now go through array and export to arrays, print to file
ytmp = {'maskvol_unwarp','maskvol_warp','mask_ovl_grp','corr_tav','corr_tsd','mot_avg_tot','mot_avg_rel','mot_max_tot','mot_max_rel','dvr_avg_tot','dvr_avg_rel','dvr_max_tot','dvr_max_rel','tav_map_av','tav_map_sd','tsd_map_av','tsd_map_sd','tsn_av','tsn_sd','krt_peak','krt_avg','krt_frac','corr_gs_rat'};
ytst = {'norm',          'norm',        '1-gam',       '1-gam',   '1-gam',   'gam',        'gam',        'gam',        'gam',        'gam',        'gam',        'gam',        'gam',        'norm',      'norm',      'norm',      'norm',      'norm-', 'norm-', 'gam',     'gam',    'gam',     'gam'    };
yid  = {'ID'};
% functional first
kq=0; clear catid catmat;
for ns=1:numel(QCStruct_quant)
    for nr=1:QCStruct_quant(ns).N_func
        kq=kq+1;
        xtmp   = QCStruct_quant(ns).frun(nr);
        catid{kq} = strcat(QCStruct_quant(ns).PREFIX,'_run(',num2str(nr),')');
        cattmp = [];
        for iu=1:numel(ytmp)
            cattmp = [cattmp xtmp.(ytmp{iu})];
        end
        catmat(kq,:) = [cattmp];
    end
end
catid = pad(catid); yid{1}=pad('ID',numel(catid{1})); % pad out the IDs
filo = fopen( fullfile(outpath,'_group_level','QC','qc.quant','table_func_stats.txt'),'w');
fprintf(filo,'%s | %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s\n',yid{1},ytmp{:});
for i=1:size(catmat,1)
    fprintf(filo,'%s | %15.2f %15.2f %15.2f %15.2f %15.2f %15.2f %15.2f %15.2f %15.2f %15.2f %15.2f %15.2f %15.2f %15.2f %15.2f %15.2f %15.2f %15.2f %15.2f %15.2f %15.2f %15.2f %15.2f\n',catid{i},catmat(i,:));
end
fclose(filo);

out = batch_outlier_testing( catmat, ytst, 0.05, 'FDR' );

% % filo = fopen( fullfile(outpath,'_group_level','QC','qc.quant','outlier_func_stats.txt'),'w');
% % ko=0;
% % for i=1:size(catmat,2)
% %     xmed = median(catmat(:,i));
% %     ix = find( abs(catmat(:,i)-xmed) > 0.1*abs(xmed) );
% %     if ~isempty(ix)
% %         for j = 1:numel(ix)
% %             ko=ko+1;
% %             fprintf(filo,'%u. Line=%u/%u, ID=%s, metric=%s:  value is %.2f, median is %.2f\n',ko,ix(j),size(catmat,1),catid{ix(j)},ytmp{i}, catmat(ix(j),i),xmed);
% %         end
% %     end
% % end
% % fprintf(filo,'\n==> A total of %u functional outlier instances found!\n',ko);
% % fclose(filo);

% now go through array and export to arrays, print to file
ytmp = {'maskvol_unwarp','maskvol_warp','mask_ovl_grp','corr_brain','corr_csf','corr_gm','corr_wm','brain_av','brain_sd'};
yid  = {'ID'};
% anatomical nirst
kq=0; clear catid catmat;
for ns=1:numel(QCStruct_quant)
    for nr=1 %%%:QCStruct_quant(ns).N_anat
        kq=kq+1;
        xtmp   = QCStruct_quant(ns).arun(nr);
        catid{kq} = strcat(QCStruct_quant(ns).PREFIX,'_run(',num2str(nr),')');
        cattmp = [];
        for iu=1:numel(ytmp)
            cattmp = [cattmp xtmp.(ytmp{iu})];
        end
        catmat(kq,:) = [cattmp];
    end
end
catid = pad(catid); yid{1}=pad('ID',numel(catid{1})); % pad out the IDs
filo = fopen( fullfile(outpath,'_group_level','QC','qc.quant','table_anat_stats.txt'),'w');
fprintf(filo,'%s | %15s %15s %15s %15s %15s %15s %15s %15s %15s\n',yid{1},ytmp{:});
for i=1:size(catmat,1)
    fprintf(filo,'%s | %15.2f %15.2f %15.2f %15.2f %15.2f %15.2f %15.2f %15.2f %15.2f\n',catid{i},catmat(i,:));
end
fclose(filo);
%
% % filo = fopen( fullfile(outpath,'_group_level','QC','qc.quant','outlier_anat_stats.txt'),'w');
% % ko=0;
% % for i=1:size(catmat,2)
% %     xmed = median(catmat(:,i));
% %     ix = find( abs(catmat(:,i)-xmed) > 0.1*abs(xmed) );
% %     if ~isempty(ix)
% %         for j = 1:numel(ix)
% %             ko=ko+1;
% %             fprintf(filo,'%u. Line=%u/%u, ID=%s, metric=%s:  value is %.2f, median is %.2f\n',ko,ix(j),size(catmat,1),catid{ix(j)},ytmp{i}, catmat(ix(j),i),xmed);
% %         end
% %     end
% % end
% % fprintf(filo,'\n==> A total of %u anatomical outlier instances found!\n',ko);
% % fclose(filo);
