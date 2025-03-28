function P3_perf_motion_diagnostics( inputfile, pipelinefile, paramlist, outpath, qc_subj_idxes )
%

% declaring path
CODE_PATH = fileparts(which('P0_perf_populateDirectories.m'));
if CODE_PATH(end)~='/'
    CODE_PATH = [CODE_PATH '/'];
end

if nargin<5
    qc_subj_idxes=[];
end

% initializing structure and running checkes
[subject_list, InputStruct_aug, PipeStruct_aug, ParamStruct_aug] = P0_perf_populateDirectories( inputfile, pipelinefile, paramlist, outpath );
% now augmenting outpath... do this after P0!
outpath = fullfile(outpath,'perf_proc'); % subdir should be fmri_proc
if ~exist(outpath,'dir') error('perf proc directory does not exist!'); end
e=dir([outpath,'*']); % dir to get absolute
if ~strcmpi( e(1).name, 'perf_proc') error('first dir should be "perf_proc"'); end
outpath = fullfile(e(1).folder,e(1).name); % convert to absolute
% check for missing files from previous step too
File_Existence_Checker_perf(InputStruct_aug,outpath,1); 

% list of subjects for constructing group masks
if isempty(qc_subj_idxes)
    disp('using all subj for final qc!')
    qc_subj_idxes = 1:numel(subject_list);
elseif ischar(qc_subj_idxes) && exist(qc_subj_idxes,'file')
    disp('loading file list for QC testing!')
    x = load(qc_subj_idxes);
    qc_subj_idxes = x.qc_subj_idxes; clear x;
elseif ~isnumeric(qc_subj_idxes)
    error('unrecognized qc id format?')
else
    disp('using numeric list of subj values for final qc construction!')
end
% store information about masking sublist...
subject_list_forqc = subject_list(qc_subj_idxes);
if exist([outpath,'/_group_level/QC/qc.quant/pipe_',PipeStruct_aug.PNAME{1},'_qc_subj_idxes.mat'],'file')
    x=load([outpath,'/_group_level/QC/qc.quant/pipe_',PipeStruct_aug.PNAME{1},'_qc_subj_idxes.mat']);
    if     ~isempty( setdiff(subject_list_forqc,x.subject_list_forqc) ) 
        error('custom subject list for qcing :: subjects in new list not present in old! delete group level folders if you want to update!')
    elseif ~isempty( setdiff(x.subject_list_forqc,subject_list_forqc) )
        error('custom subject list for qcing :: subjects not in new list that are present in old! delete group level folders if you want to update!')
    else
        disp('custom subject list for qcing :: list is consistent with old one ... continuing without modification!')
    end
else
    save([outpath,'/_group_level/QC/qc.quant/pipe_',PipeStruct_aug.PNAME{1},'_qc_subj_idxes.mat'],'qc_subj_idxes','subject_list_forqc');
end

if ~exist(  fullfile(outpath,'_group_level','QC','qc.quant',['QCStruct_quant_pipe_',PipeStruct_aug.PNAME{1},'.mat'])  ,'file')
    
    % perf
    Mref_perf = load_untouch_niiz([outpath,'/_group_level/masks/pipe_',PipeStruct_aug.PNAME{1},'/perf_brain_mask_grp.nii']);
%     V    = load_untouch_niiz([outpath,'/_group_level/brain_maps/pipe_',PipeStruct_aug.PNAME{1},'/perf_tAV_grp.nii']);
%     avec = nifti_to_mat(V,Mref_perf);
%     V    = load_untouch_niiz([outpath,'/_group_level/brain_maps/pipe_',PipeStruct_aug.PNAME{1},'/perf_tSD_grp.nii']);
%     svec = nifti_to_mat(V,Mref_perf);
    
    % anat
    Mref_anat = load_untouch_niiz([outpath,'/_group_level/masks/pipe_',PipeStruct_aug.PNAME{1},'/anat_brain_mask_grp.nii']);
    V    = load_untouch_niiz([outpath,'/_group_level/brain_maps/pipe_',PipeStruct_aug.PNAME{1},'/anat_brain_grp.nii']);
    bvec = nifti_to_mat(V,Mref_anat);
    V    = load_untouch_niiz([outpath,'/_group_level/brain_maps/pipe_',PipeStruct_aug.PNAME{1},'/anat_pCSF_grp.nii']);
    tmat(:,1) = nifti_to_mat(V,Mref_anat);
    V    = load_untouch_niiz([outpath,'/_group_level/brain_maps/pipe_',PipeStruct_aug.PNAME{1},'/anat_pGM_grp.nii']);
    tmat(:,2) = nifti_to_mat(V,Mref_anat);
    V    = load_untouch_niiz([outpath,'/_group_level/brain_maps/pipe_',PipeStruct_aug.PNAME{1},'/anat_pWM_grp.nii']);
    tmat(:,3) = nifti_to_mat(V,Mref_anat);
    
    
    for ni=1:numel(qc_subj_idxes)
    
        % check existence of subject specific struct file
        if exist(fullfile( outpath,subject_list{qc_subj_idxes(ni)},'InputStruct_ssa.mat'),'file')
            load( fullfile( outpath,subject_list{qc_subj_idxes(ni)},'InputStruct_ssa.mat') )
        else
            error('cannot find Input struct file for subject: %s \n');
        end
    
        fprintf('\n===> qc-proc. now on subj %u/%u: %s...\n',ni,numel(qc_subj_idxes),subject_list{qc_subj_idxes(ni)}),
    
        % quick formatting stuff, again assuRImes that directory structure was already constructed in "P0" pipeline step 
        opath0 = fullfile(outpath,InputStruct_ssa.PREFIX,'rawdata');
        %
        opath1a = fullfile( outpath, InputStruct_ssa.PREFIX,'anat_proc');
        opath2a = fullfile( opath1a,sprintf('subpipe_%03u',    PipeStruct_aug.pipe_idx.Warp));
        opath3a = fullfile( opath2a,sprintf('seg_subpipe_%03u',PipeStruct_aug.pipe_idx.Seg ));
        %
        opath1p = fullfile(outpath,InputStruct_ssa.PREFIX,'phys_proc');
        %
        opath1f = fullfile(outpath,InputStruct_ssa.PREFIX,'perf_proc_p1');
        opath2f = fullfile(opath1f,sprintf('subpipe_%03u',PipeStruct_aug.pipe_idx.P1));
        opath3f = fullfile(opath2f,sprintf('seg_subpipe_%03u',PipeStruct_aug.pipe_idx.Seg ));
        opath4f = fullfile(outpath,InputStruct_ssa.PREFIX,'perf_proc_p2',['pipe_',PipeStruct_aug.PNAME{1}]);
    
        QCStruct_quant(ni).PREFIX = InputStruct_ssa.PREFIX;
        QCStruct_quant(ni).N_anat = InputStruct_ssa.N_anat;
        QCStruct_quant(ni).N_perf = InputStruct_ssa.N_perf;
    
        %% -- ANATOMICAL --
    
        nr=1; % anatomical metrics
    
        M=load_untouch_niiz(sprintf('%s/anatBrainMask.nii.gz',opath2a));
        QCStruct_quant(ni).arun(nr).maskvol_unwarp = sum(M.img(:));
        M=load_untouch_niiz(sprintf('%s/anatBrainMask_warped.nii.gz',opath2a));
        QCStruct_quant(ni).arun(nr).maskvol_warp = sum(M.img(:));
        % jaccard overlap viz groupmask
        a=Mref_anat.img(:);
        b=M.img(:);
        QCStruct_quant(ni).arun(nr).mask_ovl_grp = sum( a.*b ) ./ sum( double(a | b) );
    
        V = load_untouch_niiz(sprintf('%s/anat_warped.nii.gz',opath2a));
        QCStruct_quant(ni).arun(nr).corr_brain = corr( bvec, nifti_to_mat(V,Mref_anat));
        % signal scaling
        QCStruct_quant(ni).arun(nr).brain_av = mean( nifti_to_mat(V,Mref_anat) );
        QCStruct_quant(ni).arun(nr).brain_sd =  std( nifti_to_mat(V,Mref_anat) );
        
        V = load_untouch_niiz(sprintf('%s/anat_seg_CSF_warped.nii.gz',opath3a));
        QCStruct_quant(ni).arun(nr).corr_csf = corr( tmat(:,1), nifti_to_mat(V,Mref_anat));
        QCStruct_quant(ni).arun(nr).frac_csf = mean( nifti_to_mat(V,Mref_anat) );
        V = load_untouch_niiz(sprintf('%s/anat_seg_GM_warped.nii.gz',opath3a));
        QCStruct_quant(ni).arun(nr).corr_gm = corr(  tmat(:,2), nifti_to_mat(V,Mref_anat));
        QCStruct_quant(ni).arun(nr).frac_gm = mean( nifti_to_mat(V,Mref_anat) );
        V = load_untouch_niiz(sprintf('%s/anat_seg_WM_warped.nii.gz',opath3a));
        QCStruct_quant(ni).arun(nr).corr_wm = corr(  tmat(:,3), nifti_to_mat(V,Mref_anat));
        QCStruct_quant(ni).arun(nr).frac_wm = mean( nifti_to_mat(V,Mref_anat) );
    
        %% -- PERFUSIONAL --
    
        for nr=1:InputStruct_ssa.N_perf
    
            mpe = load(sprintf('%s/warp/perf%u_mpe',opath2f,nr),'-ascii');
    
            % absolute [motion] disp
            rmsd=[];
            for i=1:size(mpe,1)-1
                rmsd = [rmsd; mean(  (mpe(i+1:end,:)-mpe(i,:)).^2, 2  ) ];
            end
            QCStruct_quant(ni).prun(nr).mot_avg_tot = mean(rmsd);
            QCStruct_quant(ni).prun(nr).mot_max_tot = max(rmsd);
            % relative [motion] framewise disp
            rmsd = mean( (mpe(2:end,:)-mpe(1:end-1,:)).^2,2);
            QCStruct_quant(ni).prun(nr).mot_avg_rel = mean(rmsd);
            QCStruct_quant(ni).prun(nr).mot_max_rel = max(rmsd);
    
            % unwarped mask [mask quality / brain shape] - total brain volume
            M=load_untouch_niiz(sprintf('%s/prewarp/motref_smo_mask.nii.gz',opath2f)); Mqp = M; %--> hold this mask aside
            QCStruct_quant(ni).prun(nr).maskvol_unwarp = sum(M.img(:));
            % warped mask [mask quality / brain shape] - total brain volume
            M=load_untouch_niiz(sprintf('%s/postwarp/m0ref_warped_mask.nii.gz',opath2f));
            QCStruct_quant(ni).prun(nr).maskvol_warp = sum(M.img(:));
            % jaccard overlap viz groupmask
            a=Mref_perf.img(:);
            b=M.img(:);
            QCStruct_quant(ni).prun(nr).mask_ovl_grp = sum( a.*b ) ./ sum( double(a | b) );
    
            %skip for now --> we aren't aligning?
%             V = load_untouch_niiz(sprintf('%s/postwarp/perf%u_warped_tav.nii.gz',opath2f,nr));
%             QCStruct_quant(ni).prun(nr).corr_tav = corr( avec, nifti_to_mat(V,Mref_perf));
%             V = load_untouch_niiz(sprintf('%s/postwarp/perf%u_warped_tsd.nii.gz',opath2f,nr));
%             QCStruct_quant(ni).prun(nr).corr_tsd = corr( svec, nifti_to_mat(V,Mref_perf));
            QCStruct_quant(ni).prun(nr).corr_tav = 1;
            QCStruct_quant(ni).prun(nr).corr_tsd = 1;
    
            %V = load_untouch_niiz([opath4f,'/perf',num2str(nr),'_fullproc.nii.gz']);
            %volmat = nifti_to_mat(V,Mref_perf);
            V = load_untouch_niiz([opath2f,'/prewarp/perf',num2str(nr),'_presmo.nii.gz']);
            volmat = nifti_to_mat(V,Mqp);
    
            % dvars [motion] disp
            rmsd=[];
            for i=1:size(volmat,2)-1
                rmsd = [rmsd; mean(  (volmat(:,i+1:end)-volmat(:,i)).^2, 1  )' ];
            end
            QCStruct_quant(ni).prun(nr).dvr_avg_tot = mean(rmsd);
            QCStruct_quant(ni).prun(nr).dvr_max_tot = max(rmsd);
            % dvars [motion] framewise disp
            rmsd = mean( (volmat(:,2:end)-volmat(:,1:end-1)).^2,1);
            QCStruct_quant(ni).prun(nr).dvr_avg_rel = mean(rmsd);
            QCStruct_quant(ni).prun(nr).dvr_max_rel = max(rmsd);
            % signal scaling
            QCStruct_quant(ni).prun(nr).tav_map_av = mean( mean(volmat,2) );
            QCStruct_quant(ni).prun(nr).tav_map_sd =  std( mean(volmat,2) );
            QCStruct_quant(ni).prun(nr).tsd_map_av = mean( std(volmat,0,2) );
            QCStruct_quant(ni).prun(nr).tsd_map_sd =  std( std(volmat,0,2) );
            % tsnr [amt of noise] average and variability
            QCStruct_quant(ni).prun(nr).tsnr_av = mean(  mean(volmat,2)./(std(volmat,0,2)+eps) );
            QCStruct_quant(ni).prun(nr).tsnr_sd =  std(  mean(volmat,2)./(std(volmat,0,2)+eps) );
            % kurtosis [presence of outliers]
            krt = kurtosis( volmat,0,2 );
            thr_krt = prctile( kurtosis( randn(size(volmat,2),5000) ), 99.5 );
            QCStruct_quant(ni).prun(nr).krt_peak = prctile(krt,99.5);
            QCStruct_quant(ni).prun(nr).krt_avg  = mean(krt);
            QCStruct_quant(ni).prun(nr).krt_frac = mean(krt>thr_krt);
            % globalness of signale
            vnrm = bsxfun(@minus,volmat, mean(volmat,2));
            vnrm = bsxfun(@rdivide,volmat,sqrt(sum(volmat.^2,2)));
            [u,~,v] =svd( vnrm,'econ');
            rho = vnrm * (v(:,1) * sign(mean(u(:,1))) );
            QCStruct_quant(ni).prun(nr).corr_gs_rat = mean(rho)/std(rho);
    
        end
    end
    
    % for ni=1:numel(subject_list)
    %     for nr=1:InputStruct_ssa.N_perf
    %         qcmatm(ni,1) = QCStruct_quant(ni).prun(nr).mot_avg_tot;
    %         qcmatm(ni,2) = QCStruct_quant(ni).prun(nr).mot_max_tot;
    %         qcmatm(ni,3) = QCStruct_quant(ni).prun(nr).mot_avg_rel;
    %         qcmatm(ni,4) = QCStruct_quant(ni).prun(nr).mot_max_rel;
    %     end
    % end
    % out = batch_outlier_testing( qcmatm, {'gam','gam','gam','gam'}, 0.05, 'FDR' );
    % figure,imagesc( out.thr)
    % 
    % for ni=1:numel(subject_list)
    %     for nr=1:InputStruct_ssa.N_perf
    %         qcmatd(ni,1) = QCStruct_quant(ni).prun(nr).dvr_avg_tot;
    %         qcmatd(ni,2) = QCStruct_quant(ni).prun(nr).dvr_max_tot;
    %         qcmatd(ni,3) = QCStruct_quant(ni).prun(nr).dvr_avg_rel;
    %         qcmatd(ni,4) = QCStruct_quant(ni).prun(nr).dvr_max_rel;
    %     end
    % end
    % out = batch_outlier_testing( qcmatd, {'gam','gam','gam','gam'}, 0.05, 'FDR' );
    % figure,imagesc( out.thr)
    
    %% OMNIBUS SUMMARY...
    
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
    save( fullfile(outpath,'_group_level','QC','qc.quant',['QCStruct_quant_pipe_',PipeStruct_aug.PNAME{1},'.mat']), 'QCStruct_quant' );

else
    disp('qcfile found, reloading...');
    %
    load( fullfile(outpath,'_group_level','QC','qc.quant',['QCStruct_quant_pipe_',PipeStruct_aug.PNAME{1},'.mat']) );
end

%% Functional summary report...

% now go through array and export to arrays, print to file
ytmp = {'maskvol_unwarp','maskvol_warp','mask_ovl_grp','corr_tav','corr_tsd','mot_avg_tot','mot_avg_rel','mot_max_tot','mot_max_rel','dvr_avg_tot','dvr_avg_rel','dvr_max_tot','dvr_max_rel','tav_map_av','tav_map_sd','tsd_map_av','tsd_map_sd','tsnr_av','tsnr_sd','krt_peak','krt_avg','krt_frac','corr_gs_rat'};
ytst = {'norm',          'norm',        '1-gam',       '1-gam',   '1-gam',   'gam',        'gam',        'gam',        'gam',        'gam',        'gam',        'gam',        'gam',        'norm',      'norm',      'norm',      'norm',      'norm-', 'norm-', 'gam',     'gam',    'gam',     'gam'    };
yid  = {'ID'};
% functional first
kq=0; clear catid catmat;
for ni=1:numel(QCStruct_quant)
    for nr=1:QCStruct_quant(ni).N_perf
        kq=kq+1;
        xtmp   = QCStruct_quant(ni).prun(nr);
        catid{kq} = strcat(QCStruct_quant(ni).PREFIX,'_run(',num2str(nr),')');
        cattmp = [];
        for iu=1:numel(ytmp)
            cattmp = [cattmp xtmp.(ytmp{iu})];
        end
        catmat(kq,:) = [cattmp];
    end
end
catid = pad(catid); yid{1}=pad('ID',numel(catid{1})); % pad out the IDs
filo = fopen( fullfile(outpath,'_group_level','QC','qc.quant',['table_perf_stats_pipe_',PipeStruct_aug.PNAME{1},'.txt']),'w');
fprintf(filo,'%s | %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s\n',yid{1},ytmp{:});
for i=1:size(catmat,1)
    fprintf(filo,'%s | %15.2f %15.2f %15.2f %15.2f %15.2f %15.2f %15.2f %15.2f %15.2f %15.2f %15.2f %15.2f %15.2f %15.2f %15.2f %15.2f %15.2f %15.2f %15.2f %15.2f %15.2f %15.2f %15.2f\n',catid{i},catmat(i,:));
end
fclose(filo);

out = batch_outlier_testing( catmat, ytst, 0.05, 'FDR' );
out.thr(out.thr<0)=0;

clear abar; legcell = {'vol/shape','mot-mpe','mot-bold','bold-scal','tsnr','kurt','globsig'};
abar(:,1) = sum(out.thr(:,[ 1: 5]),2);
abar(:,2) = sum(out.thr(:,[ 6: 9]),2);
abar(:,3) = sum(out.thr(:,[10:13]),2);
abar(:,4) = sum(out.thr(:,[14:17]),2);
abar(:,5) = sum(out.thr(:,[18:19]),2);
abar(:,6) = sum(out.thr(:,[20:22]),2);
abar(:,7) = sum(out.thr(:,[   23]),2);
figure, bar( abar,'stacked' ); ylim([0 10]);
title('outlier counts - perfusional data');
legend(legcell);

ix = find( sum(abar,2)>0 );

fprintf('\n\nTotal of %u functional runs with outlier values (FDR=0.05):\n',numel(ix)),
if ~isempty(ix)
    for i=1:numel(ix)
        stro=[];
        ix2 = find( abar(ix(i),:)>0 );
        for j=1:numel(ix2)
            stro = [stro, sprintf(' %u (%s) /',abar(ix(i),ix2(j)), legcell{ix2(j)})];
        end
        fprintf('\t%u. %s with: %s\n',i,catid{ix(i)},stro(1:end-1));
    end
end

%% Anatomical summary report...

% now go through array and export to arrays, print to file
ytmp = {'maskvol_unwarp','maskvol_warp','mask_ovl_grp','corr_brain','corr_csf','corr_gm','corr_wm','brain_av','brain_sd'};
ytst = {'norm',          'norm',        '1-gam',       '1-gam',     '1-gam',   '1-gam',  '1-gam',  'norm',    'norm'    };
yid  = {'ID'};
% anatomical nirst
kq=0; clear catid catmat;
for ni=1:numel(QCStruct_quant)
    for nr=1 %%%:QCStruct_quant(ni).N_anat
        kq=kq+1;
        xtmp   = QCStruct_quant(ni).arun(nr);
        catid{kq} = strcat(QCStruct_quant(ni).PREFIX,'_run(',num2str(nr),')');
        cattmp = [];
        for iu=1:numel(ytmp)
            cattmp = [cattmp xtmp.(ytmp{iu})];
        end
        catmat(kq,:) = [cattmp];
    end
end
catid = pad(catid); yid{1}=pad('ID',numel(catid{1})); % pad out the IDs
filo = fopen( fullfile(outpath,'_group_level','QC','qc.quant',['table_anat_stats_pipe_',PipeStruct_aug.PNAME{1},'.txt']),'w');
fprintf(filo,'%s | %15s %15s %15s %15s %15s %15s %15s %15s %15s\n',yid{1},ytmp{:});
for i=1:size(catmat,1)
    fprintf(filo,'%s | %15.2f %15.2f %15.2f %15.2f %15.2f %15.2f %15.2f %15.2f %15.2f\n',catid{i},catmat(i,:));
end
fclose(filo);

out = batch_outlier_testing( catmat, ytst, 0.05, 'FDR' );
out.thr(out.thr<0)=0;

clear abar; legcell = {'vol/shape','segment','sig-scal'};
abar(:,1) = sum(out.thr(:,[ 1: 4]),2);
abar(:,2) = sum(out.thr(:,[ 5: 7]),2);
abar(:,3) = sum(out.thr(:,[ 8: 9]),2);
figure, bar( abar,'stacked' ); ylim([0 10]);
title('outlier counts - structural data');
legend(legcell);

ix = find( sum(abar,2)>0 );

fprintf('\n\nTotal of %u anatomical runs with outlier values (FDR=0.05):\n',numel(ix)),
if ~isempty(ix)
    for i=1:numel(ix)
        stro=[];
        ix2 = find( abar(ix(i),:)>0 );
        for j=1:numel(ix2)
            stro = [stro, sprintf(' %u (%s) /',abar(ix(i),ix2(j)), legcell{ix2(j)})];
        end
        fprintf('\t%u. %s with: %s\n',i,catid{ix(i)},stro(1:end-1));
    end
end
