close all;
clear;

% ==> local texture analysis running robustness checkses

list_s = [10:13 17:18, 49:54, 16];
list_c = [11101:11141 11143:11175, 12101:12141 12143:12175];
masklab = [list_s list_c];

subusblist=[72 14 120];


% hardcode absolute path
pathloc_fs = 'Neurocovid/NCOV_SUBS';
pathloc_af = 'Neurocovid/fmri_proc';
e = dir(sprintf('%s/sub-*-V01',pathloc_fs));
if isempty(e)
    error('no folders found!');
end

return;

% %running texture analysis, optimization (subcortical only)
% 
% for i=1:numel(e)
%     if ~strcmpi(e(i).name,'sub-NCOV1F006_ses-V01')
%         prfix = e(i).name;    
%         prfix,
% 
%         t1_file = sprintf('texturial_data/%s_t1.nii.gz',prfix);
%         mskfile = sprintf('texturial_data/%s_seg.nii.gz',prfix);
%         outfile = sprintf('texturial_matfiles_defrange/%s',prfix);
% 
%         texmap_wrapper3( t1_file, mskfile, [list_s], 0, outfile, [], [], [], [], 'KDE-G', 'optim' );
%     else
%         disp('not on list!');
%     end
% end

% %running texture analysis, optimization (subcortical only; different range settings)
% 
% for i=1:numel(e)
%     if ~strcmpi(e(i).name,'sub-NCOV1F006_ses-V01')
%         prfix = e(i).name;    
%         prfix,
% 
%         t1_file = sprintf('texturial_data/%s_t1.nii.gz',prfix);
%         mskfile = sprintf('texturial_data/%s_seg.nii.gz',prfix);
%         outfile = sprintf('texturial_matfiles_5pctrange/%s',prfix);
% 
%         texmap_wrapper3( t1_file, mskfile, [list_s], 0, outfile, [], [], ['PctRange-5'], [], 'KDE-G', 'optim' );
%     else
%         disp('not on list!');
%     end
% end

% %running texture analysis, fixed BW value (cortical only
% for i=1:numel(e)
%     if ~strcmpi(e(i).name,'sub-NCOV1F006_ses-V01')
%         prfix = e(i).name;    
%         prfix,
% 
%         t1_file = sprintf('texturial_data/%s_t1.nii.gz',prfix);
%         mskfile = sprintf('texturial_data/%s_seg.nii.gz',prfix);
%         outfile = sprintf('texturial_matfiles_defrange2/%s',prfix);
% 
%         texmap_wrapper3( t1_file, mskfile, [list_c], 0, outfile, [], [], [], [], 'KDE-G', [2.7*10^-4] );
%     else
%         disp('not on list!');
%     end
% end

%running texture analysis, optimization (subcortical only)

% for i=1:numel(e)
%     if ~strcmpi(e(i).name,'sub-NCOV1F006_ses-V01')
%         prfix = e(i).name;    
%         prfix,
% 
%         t1_file = sprintf('texturial_data/%s_t1.nii.gz',prfix);
%         mskfile = sprintf('texturial_data/%s_seg.nii.gz',prfix);
%         outfile = sprintf('texturial_matfiles_histo/%s',prfix);
% 
%         texmap_wrapper3( t1_file, mskfile, [masklab([1 40 80 120 160])], 0, outfile, [], [], [], [], 'HDE', 'optim' );
%     else
%         disp('not on list!');
%     end
% end


for i=1:numel(e)
    if ~strcmpi(e(i).name,'sub-NCOV1F006_ses-V01')
        prfix = e(i).name;    
        prfix,

        t1_file = sprintf('texturial_data/%s_t1.nii.gz',prfix);
        mskfile = sprintf('texturial_data/%s_seg.nii.gz',prfix);
        outfile = sprintf('texturial_matfiles_kerno/%s',prfix);

        texmap_wrapper3( t1_file, mskfile, [masklab([1 20 40 60 80 100 120 140 160])], 0, outfile, [], [], [], [], 'KDE-G', 'optim' ); %  
    else
        disp('not on list!');
    end
end

return;
%collecting datas1
kq=0;
for i=1:numel(e)
    if ~strcmpi(e(i).name,'sub-NCOV1F006_ses-V01')
        kq=kq+1;
        prfix = e(i).name; 
        prfix,
        outfile_full = sprintf('texturial_matfiles_histo/%s_Output_Roi_byparc.mat',prfix);
        load(outfile_full);
        sse_h(kq,:,:) = permute( out.opt_stats, [3 2 1]);
    else
        disp('not on list!');
    end
end

kq=0;
for i=1:numel(e)
    if ~strcmpi(e(i).name,'sub-NCOV1F006_ses-V01')
        kq=kq+1;
        prfix = e(i).name;    
        prfix,
        outfile_full = sprintf('texturial_matfiles_kerno/%s_Output_Roi_byparc.mat',prfix);
        load(outfile_full);
        sse_k(kq,:,:) = permute( out.opt_stats, [3 2 1]);
    else
        disp('not on list!');
    end
end

return;
%%
%% STARTASRARET
%%

% fixed BW, BASE SETTING
for i=1:numel(e)
    if ~strcmpi(e(i).name,'sub-NCOV1F006_ses-V01')
        prfix = e(i).name;    
        prfix,
        t1_file = sprintf('texturial_data/%s_t1.nii.gz',prfix);
        mskfile = sprintf('texturial_data/%s_seg.nii.gz',prfix);
        outfile = sprintf('texturial_matfiles_all/rebin_h/%s_fix',prfix);
        texmap_wrapper3( t1_file, mskfile, [masklab(subusblist)], 0, outfile, [], [], [], [], 'HDE', 0.048 );
    else
        disp('not on list!');
    end
end
for i=1:numel(e)
    if ~strcmpi(e(i).name,'sub-NCOV1F006_ses-V01')
        prfix = e(i).name;    
        prfix,
        t1_file = sprintf('texturial_data/%s_t1.nii.gz',prfix);
        mskfile = sprintf('texturial_data/%s_seg.nii.gz',prfix);
        outfile = sprintf('texturial_matfiles_all/rebin_t/%s_fix',prfix);
        texmap_wrapper3( t1_file, mskfile, [masklab(subusblist)], 0, outfile, [], [], [], [], 'KDE-GT', 0.023 );
    else
        disp('not on list!');
    end
end


% fixed BW, narrow-bw
for i=1:numel(e)
    if ~strcmpi(e(i).name,'sub-NCOV1F006_ses-V01')
        prfix = e(i).name;    
        prfix,
        t1_file = sprintf('texturial_data/%s_t1.nii.gz',prfix);
        mskfile = sprintf('texturial_data/%s_seg.nii.gz',prfix);
        outfile = sprintf('texturial_matfiles_all/rebin_h/%s_bw3Lo',prfix);
        texmap_wrapper3( t1_file, mskfile, [masklab(subusblist)], 0, outfile, [], [], [], [], 'HDE', 0.048/3 );
    else
        disp('not on list!');
    end
end
for i=1:numel(e) % [9 14 16 26
    if ~strcmpi(e(i).name,'sub-NCOV1F006_ses-V01')
        prfix = e(i).name;    
        prfix,
        t1_file = sprintf('texturial_data/%s_t1.nii.gz',prfix);
        mskfile = sprintf('texturial_data/%s_seg.nii.gz',prfix);
        outfile = sprintf('texturial_matfiles_all/rebin_t/%s_bw3Lo',prfix);
        texmap_wrapper3( t1_file, mskfile, [masklab(subusblist)], 0, outfile, [], [], [], [], 'KDE-GT', 0.023/3 );
    else
        disp('not on list!');
    end
end
% fixed BW, wide-bw
for i=1:numel(e)
    if ~strcmpi(e(i).name,'sub-NCOV1F006_ses-V01')
        prfix = e(i).name;    
        prfix,
        t1_file = sprintf('texturial_data/%s_t1.nii.gz',prfix);
        mskfile = sprintf('texturial_data/%s_seg.nii.gz',prfix);
        outfile = sprintf('texturial_matfiles_all/rebin_h/%s_bw3Hi',prfix);
        texmap_wrapper3( t1_file, mskfile, [masklab(subusblist)], 0, outfile, [], [], [], [], 'HDE', 0.048*3 );
    else
        disp('not on list!');
    end
end
for i=1:numel(e)
    if ~strcmpi(e(i).name,'sub-NCOV1F006_ses-V01')
        prfix = e(i).name;    
        prfix,
        t1_file = sprintf('texturial_data/%s_t1.nii.gz',prfix);
        mskfile = sprintf('texturial_data/%s_seg.nii.gz',prfix);
        outfile = sprintf('texturial_matfiles_all/rebin_t/%s_bw3Hi',prfix);
        texmap_wrapper3( t1_file, mskfile, [masklab(subusblist)], 0, outfile, [], [], [], [], 'KDE-GT', 0.023*3 );
    else
        disp('not on list!');
    end
end


% fixed BW, narrow-rng
for i=1:numel(e)
    if ~strcmpi(e(i).name,'sub-NCOV1F006_ses-V01')
        prfix = e(i).name;    
        prfix,
        t1_file = sprintf('texturial_data/%s_t1.nii.gz',prfix);
        mskfile = sprintf('texturial_data/%s_seg.nii.gz',prfix);
        outfile = sprintf('texturial_matfiles_all/rebin_h/%s_rngLo',prfix);
        texmap_wrapper3( t1_file, mskfile, [masklab(subusblist)], 0, outfile, [], [], ['PctRange-5'], [], 'HDE', 0.048 );
    else
        disp('not on list!');
    end
end
for i=1:numel(e)
    if ~strcmpi(e(i).name,'sub-NCOV1F006_ses-V01')
        prfix = e(i).name;    
        prfix,
        t1_file = sprintf('texturial_data/%s_t1.nii.gz',prfix);
        mskfile = sprintf('texturial_data/%s_seg.nii.gz',prfix);
        outfile = sprintf('texturial_matfiles_all/rebin_t/%s_rngLo',prfix);
        texmap_wrapper3( t1_file, mskfile, [masklab(subusblist)], 0, outfile, [], [], ['PctRange-5'], [], 'KDE-GT', 0.023 );
    else
        disp('not on list!');
    end
end
% fixed BW, wide-rng
for i=1:numel(e)
    if ~strcmpi(e(i).name,'sub-NCOV1F006_ses-V01')
        prfix = e(i).name;    
        prfix,
        t1_file = sprintf('texturial_data/%s_t1.nii.gz',prfix);
        mskfile = sprintf('texturial_data/%s_seg.nii.gz',prfix);
        outfile = sprintf('texturial_matfiles_all/rebin_h/%s_rngHi',prfix);
        texmap_wrapper3( t1_file, mskfile, [masklab(subusblist)], 0, outfile, [], [], ['PctRange-20'], [], 'HDE', 0.048 );
    else
        disp('not on list!');
    end
end
for i=1:numel(e)
    if ~strcmpi(e(i).name,'sub-NCOV1F006_ses-V01')
        prfix = e(i).name;    
        prfix,
        t1_file = sprintf('texturial_data/%s_t1.nii.gz',prfix);
        mskfile = sprintf('texturial_data/%s_seg.nii.gz',prfix);
        outfile = sprintf('texturial_matfiles_all/rebin_t/%s_rngHi',prfix);
        texmap_wrapper3( t1_file, mskfile, [masklab(subusblist)], 0, outfile, [], [], ['PctRange-20'], [], 'KDE-GT', 0.023 );
    else
        disp('not on list!');
    end
end


% fixed BW, +noise
for i=1:numel(e)
    if ~strcmpi(e(i).name,'sub-NCOV1F006_ses-V01')
        prfix = e(i).name;    
        prfix,
        t1_file = sprintf('texturial_data/%s_t1.nii.gz',prfix);
        mskfile = sprintf('texturial_data/%s_seg.nii.gz',prfix);
        outfile = sprintf('texturial_matfiles_all/rebin_h/%s_shf05pct',prfix);
        texmap_wrapper3( t1_file, mskfile, [masklab(subusblist)], 0, outfile, [], [], [], [], 'HDE', 0.048 );
    else
        disp('not on list!');
    end
end
for i=1:numel(e)
    if ~strcmpi(e(i).name,'sub-NCOV1F006_ses-V01')
        prfix = e(i).name;    
        prfix,
        t1_file = sprintf('texturial_data/%s_t1.nii.gz',prfix);
        mskfile = sprintf('texturial_data/%s_seg.nii.gz',prfix);
        outfile = sprintf('texturial_matfiles_all/rebin_t/%s_shf05pct',prfix);
        texmap_wrapper3( t1_file, mskfile, [masklab(subusblist)], 0, outfile, [], [], [], [], 'KDE-GT', 0.023 );
    else
        disp('not on list!');
    end
end

%% interrogate?


for i=1%:numel(e)
    if ~strcmpi(e(i).name,'sub-NCOV1F006_ses-V01')
        prfix = e(i).name;    
        prfix,
        t1_file = sprintf('texturial_data/%s_t1.nii.gz',prfix);
        mskfile = sprintf('texturial_data/%s_seg.nii.gz',prfix);
        outfile = sprintf('texturial_matfiles_all/rebin_t/%s_xx',prfix);
        texmap_wrapper3( t1_file, mskfile, [masklab(subusblist)], 0, outfile, [], [], [], [], 'KDE-GT', 0.023,1 );
        texmap_wrapper3( t1_file, mskfile, [masklab(subusblist)], 0, outfile, [], [], ['PctRange-20'], [], 'KDE-GT', 0.023,1 );
    else
        disp('not on list!');
    end
end

%collecting datas2
kq=0;
for i=1:numel(e)
    if ~strcmpi(e(i).name,'sub-NCOV1F006_ses-V01')
        kq=kq+1;
        prfix = e(i).name; 
        prfix,
        outfile_full = sprintf('texturial_matfiles_all/rebin_h/%s_fix_Output_Roi_byparc.mat',prfix);
        load(outfile_full);
        met_h{1}(kq,:,:) = permute( out.metrics_av, [3 2 1]);
        outfile_full = sprintf('texturial_matfiles_all/rebin_h/%s_bw3Lo_Output_Roi_byparc.mat',prfix);
        load(outfile_full);
        met_h{2}(kq,:,:) = permute( out.metrics_av, [3 2 1]);
        outfile_full = sprintf('texturial_matfiles_all/rebin_h/%s_bw3Hi_Output_Roi_byparc.mat',prfix);
        load(outfile_full);
        met_h{3}(kq,:,:) = permute( out.metrics_av, [3 2 1]);
%         outfile_full = sprintf('texturial_matfiles_all/rebin_h/%s_rngLo_Output_Roi_byparc.mat',prfix);
%         load(outfile_full);
%         met_h{4}(kq,:,:) = permute( out.metrics_av, [3 2 1]);
%         outfile_full = sprintf('texturial_matfiles_all/rebin_h/%s_rngHi_Output_Roi_byparc.mat',prfix);
%         load(outfile_full);
%         met_h{5}(kq,:,:) = permute( out.metrics_av, [3 2 1]);
%         outfile_full = sprintf('texturial_matfiles_all/rebin_h/%s_noi10pct_Output_Roi_byparc.mat',prfix);
%         load(outfile_full);
%         met_h{6}(kq,:,:) = permute( out.metrics_av, [3 2 1]);
%         outfile_full = sprintf('texturial_matfiles_all/rebin_h/%s_noi20pct_Output_Roi_byparc.mat',prfix);
%         load(outfile_full);
%         met_h{7}(kq,:,:) = permute( out.metrics_av, [3 2 1]);
%         outfile_full = sprintf('texturial_matfiles_all/rebin_h/%s_noi30pct_Output_Roi_byparc.mat',prfix);
%         load(outfile_full);
%         met_h{8}(kq,:,:) = permute( out.metrics_av, [3 2 1]);
%         outfile_full = sprintf('texturial_matfiles_all/rebin_h/%s_shf05pct_Output_Roi_byparc.mat',prfix);
%         load(outfile_full);
%         met_h{9}(kq,:,:) = permute( out.metrics_av, [3 2 1]);
    else
        disp('not on list!');
    end
end

kq=0;
for i=1:numel(e)
    if ~strcmpi(e(i).name,'sub-NCOV1F006_ses-V01')
        kq=kq+1;
        prfix = e(i).name;    
        prfix,
        outfile_full = sprintf('texturial_matfiles_all/rebin_t/%s_fix_Output_Roi_byparc.mat',prfix);
        load(outfile_full);
        met_k{1}(kq,:,:) = permute( out.metrics_av, [3 2 1]);
        outfile_full = sprintf('texturial_matfiles_all/rebin_t/%s_bw3Lo_Output_Roi_byparc.mat',prfix);
        load(outfile_full);
        met_k{2}(kq,:,:) = permute( out.metrics_av, [3 2 1]);
        outfile_full = sprintf('texturial_matfiles_all/rebin_t/%s_bw3Hi_Output_Roi_byparc.mat',prfix);
        load(outfile_full);
        met_k{3}(kq,:,:) = permute( out.metrics_av, [3 2 1]);
%         outfile_full = sprintf('texturial_matfiles_all/rebin_t/%s_rngLo_Output_Roi_byparc.mat',prfix);
%         load(outfile_full);
%         met_k{4}(kq,:,:) = permute( out.metrics_av, [3 2 1]);
%         outfile_full = sprintf('texturial_matfiles_all/rebin_t/%s_rngHi_Output_Roi_byparc.mat',prfix);
%         load(outfile_full);
%         met_k{5}(kq,:,:) = permute( out.metrics_av, [3 2 1]);
%         outfile_full = sprintf('texturial_matfiles_all/rebin_t/%s_noi10pct_Output_Roi_byparc.mat',prfix);
%         load(outfile_full);
%         met_k{6}(kq,:,:) = permute( out.metrics_av, [3 2 1]);
%         outfile_full = sprintf('texturial_matfiles_all/rebin_t/%s_noi20pct_Output_Roi_byparc.mat',prfix);
%         load(outfile_full);
%         met_k{7}(kq,:,:) = permute( out.metrics_av, [3 2 1]);
%         outfile_full = sprintf('texturial_matfiles_all/rebin_t/%s_noi30pct_Output_Roi_byparc.mat',prfix);
%         load(outfile_full);
%         met_k{8}(kq,:,:) = permute( out.metrics_av, [3 2 1]);
%         outfile_full = sprintf('texturial_matfiles_all/rebin_t/%s_shf05pct_Output_Roi_byparc.mat',prfix);
%         load(outfile_full);
%         met_k{9}(kq,:,:) = permute( out.metrics_av, [3 2 1]);
    else
        warning('not on list!');
    end
end
idlist = 1:70;
% k1=1;
% k2=5;
% figure;
% for r=1:5
%     dev(:,1) = median( 100*abs(met_h{k2}(:,:,r)-met_h{k1}(:,:,r))./abs(met_h{k1}(:,:,r)) );
%     dev(:,2) = median( 100*abs(met_k{k2}(:,:,r)-met_k{k1}(:,:,r))./abs(met_k{k1}(:,:,r)) );
%     subplot(2,5,r); bar( (dev) );
% end

k1=1;
k2=2;
figure;
for r=1:5
    dev(:,1) = diag(corr(met_h{k1}(idlist,:,r),met_h{k2}(idlist,:,r),'type','Spearman'));
    dev(:,2) = diag(corr(met_k{k1}(idlist,:,r),met_k{k2}(idlist,:,r),'type','Spearman'));
    subplot(2,5,r); bar( (dev) ); ylim([0 1.1]);
end





