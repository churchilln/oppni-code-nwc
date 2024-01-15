close all;
clear;

% ==> sample scripts for running texture analysis locally

list_s = [10:13 17:18, 49:54, 16];
list_sp= [10:13 17:18 26 28, 49:54 58 60, 16];
list_c = [11101:11141 11143:11175, 12101:12141 12143:12175];
masklab = [list_s list_c];

%lists:

list1 = [2 41]; % white matter
list2 = [10:13 17:18, 49:54]; % subcortical
list3 = [11101:11141 11143:11175, 12101:12141 12143:12175]; % cortical
list4 = [11164 12164]; % iolfactory
list5 = [11106:11110 12106:12110]; % cingulate
list6 = [11117:11118 11148:11150, 12117:12118 12148:12150]; % insula
list7 = [11101:11141 11143:11175, 12101:12141 12143:12175]; % putamen


% hardcode absolute path
pathloc_fs = 'Neurocovid/NCOV_SUBS';
pathloc_af = 'Neurocovid/fmri_proc';
e = dir(sprintf('%s/sub-*-V01',pathloc_fs));
if isempty(e)
    error('no folders found!');
end

% % for a quick plot...
% for i=1:numel(e) %list1(WM):  RAD=5,B=0.175 / RAD=7,B=0.10 / RAD=9,B=0.070 / RAD=11, B=0.055
%     if ~strcmpi(e(i).name,'sub-NCOV1F006_ses-V01')
%         prfix = e(i).name;    
%         prfix,
% 
%         t1_file = sprintf('texturial_data/%s_t1.nii.gz',prfix);
%         mskfile = sprintf('texturial_data/%s_seg.nii.gz',prfix);
%         outfile = sprintf('texturial_matfiles_vox/%s_scx',prfix);
%         texmap_wrapper3( t1_file, mskfile, list_sp, 1, outfile, 7, 7, [], [], 'HDE', ['estim'],1,1 ); %[0.0425]
%     else
%         disp('not on list!');
%     end
% end

% -- bandwidth estimation loop below -- (0.101?

                 %list3(GM):  RAD=5,B=0.160 / RAD=7,B=0.10 / RAD=9,B=0.073 RAD=11, B=0.055 
for i=1:numel(e) %list1(WM):  RAD=5,B=0.175 / RAD=7,B=0.10 / RAD=9,B=0.070 / RAD=11, B=0.055
    if ~strcmpi(e(i).name,'sub-NCOV1F006_ses-V01')
        prfix = e(i).name;    
        prfix,

        t1_file = sprintf('texturial_data/%s_t1.nii.gz',prfix);
        mskfile = sprintf('texturial_data/%s_seg.nii.gz',prfix);
        outfile = sprintf('texturial_matfiles_moar2/%s_wm',prfix);

        texmap_wrapper3( t1_file, mskfile, list_sp, 1, outfile, 7, 7, [], [], 'HDE', 'estim',0,1 ); %[0.0425]
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
        outfile = sprintf('texturial_matfiles_moar2/%s_wm_Output_Vox_restrict.mat',prfix);
        load(outfile);
        optsz(:,kq) = out.opt_stats(:,3);
    else
        disp('not on list!');
    end
end
median(optsz(:)),

return;

% -- now running batching -- %

for i=1:numel(e) %list1(WM):  RAD=5,B=0.175 / RAD=7,B=0.10 / RAD=9,B=0.070 / RAD=11, B=0.055
    if ~strcmpi(e(i).name,'sub-NCOV1F006_ses-V01')
        prfix = e(i).name;    
        prfix,

        t1_file = sprintf('texturial_data/%s_t1.nii.gz',prfix);
        mskfile = sprintf('texturial_data/%s_seg.nii.gz',prfix);
%         outfile = sprintf('texturial_matfiles_vox/%s_sc',prfix);
%         texmap_wrapper3( t1_file, mskfile, list_s, 1, outfile, 7, 7, [], [], 'HDE', [0.0909],0,1 ); %[0.0425]
%         outfile = sprintf('texturial_matfiles_vox/%s_sc1',prfix);
%         texmap_wrapper3( t1_file, mskfile, list_sp, 1, outfile, 7, 7, [], [1], 'HDE', [0.10],0,1 ); %[0.0425]
        outfile = sprintf('texturial_matfiles_vox/%s_sc2',prfix);
        texmap_wrapper3( t1_file, mskfile, list_sp, 1, outfile, 7, 7, [], [2], 'HDE', [0.10],0,1 ); %[0.0425]
        outfile = sprintf('texturial_matfiles_vox/%s_sc3',prfix);
        texmap_wrapper3( t1_file, mskfile, list_sp, 1, outfile, 7, 7, [], [3], 'HDE', [0.10],0,1 ); %[0.0425]
    else
        disp('not on list!');
    end
end


return;

for i=1:numel(e)
    if ~strcmpi(e(i).name,'sub-NCOV1F006_ses-V01')
        prfix = e(i).name;    
        prfix,

        t1_file = sprintf('texturial_data/%s_t1.nii.gz',prfix);
        mskfile = sprintf('texturial_data/%s_seg.nii.gz',prfix);
        outfile = sprintf('texturial_matfiles_moar/%s_a1',prfix);

        texmap_wrapper3( t1_file, mskfile, list1, 0, outfile, [], [], [], [], 'HDE', 0.05,0,1 ); %[0.0425]
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
        outfile = sprintf('texturial_matfiles_moar/%s_a2',prfix);

        texmap_wrapper3( t1_file, mskfile, list2, 0, outfile, [], [], [], [], 'HDE', 0.05,0,1 ); %[0.0425]
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
        outfile = sprintf('texturial_matfiles_moar/%s_a3',prfix);

        texmap_wrapper3( t1_file, mskfile, list3, 0, outfile, [], [], [], [], 'HDE', 0.05,0,1 ); %[0.0425]
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
        outfile = sprintf('texturial_matfiles_moar/%s_a4',prfix);

        texmap_wrapper3( t1_file, mskfile, list4, 0, outfile, [], [], [], [], 'HDE', 0.05,0,0 ); %[0.0425]
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
        outfile = sprintf('texturial_matfiles_moar/%s_a5',prfix);

        texmap_wrapper3( t1_file, mskfile, list5, 0, outfile, [], [], [], [], 'HDE', 0.05,0,1 ); %[0.0425]
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
        outfile = sprintf('texturial_matfiles_moar/%s_a6',prfix);

        texmap_wrapper3( t1_file, mskfile, list6, 0, outfile, [], [], [], [], 'HDE', 0.05,0,1 ); %[0.0425]
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
        outfile = sprintf('texturial_matfiles_moar/%s_a7',prfix);

        texmap_wrapper3( t1_file, mskfile, list7, 0, outfile, [], [], [], [], 'HDE', 0.05,0,0 ); %[0.0425]
    else
        disp('not on list!');
    end
end

% > redo olfac (both sides)
for i=1:numel(e)
    if ~strcmpi(e(i).name,'sub-NCOV1F006_ses-V01')
        prfix = e(i).name;    
        prfix,

        t1_file = sprintf('texturial_data/%s_t1.nii.gz',prfix);
        mskfile = sprintf('texturial_data/%s_seg.nii.gz',prfix);
        outfile = sprintf('texturial_matfiles_moar/%s_a8',prfix);

        texmap_wrapper3( t1_file, mskfile, list4, 0, outfile, [], [], [], [], 'HDE', 0.05,0,0 ); %[0.0425]
    else
        disp('not on list!');
    end
end









% >merged --> 0.08 / separate --> 
for i=1:numel(e)
    if ~strcmpi(e(i).name,'sub-NCOV1F006_ses-V01')
        prfix = e(i).name;    
        prfix,

        t1_file = sprintf('texturial_data/%s_t1.nii.gz',prfix);
        mskfile = sprintf('texturial_data/%s_seg.nii.gz',prfix);
        outfile = sprintf('texturial_matfiles_moar/%s_a4opt',prfix);

        texmap_wrapper3( t1_file, mskfile, list4, 0, outfile, [], [], [], [], 'HDE', 'optim',0,0 ); %[0.0425]
    else
        disp('not on list!');
    end
end

%--plotteur
kq=0;
for i=1:numel(e)
    if ~strcmpi(e(i).name,'sub-NCOV1F006_ses-V01')
        prfix = e(i).name;    
        prfix,
        kq=kq+1;
        outfile = sprintf('texturial_matfiles_moar/%s_a4opt_Output_Roi_byparc.mat',prfix);
        load(outfile);
        Bset(kq,1) = out.Bvl_for_fitt;
        cvg(:,kq) = mean(out.CVfull,2);
    else
        disp('not on list!');
    end
end

%--plotteur
for i=1:numel(e)
    if ~strcmpi(e(i).name,'sub-NCOV1F006_ses-V01')
        prfix = e(i).name;    
        prfix,

        t1_file = sprintf('texturial_data/%s_t1.nii.gz',prfix);
        mskfile = sprintf('texturial_data/%s_seg.nii.gz',prfix);
        outfile = sprintf('texturial_matfiles_moar/%s_a4xxxxxx',prfix);

        texmap_wrapper3( t1_file, mskfile, list4, 0, outfile, [], [], [], [], 'HDE', 0.08,1,1 ); %[0.0425]
    else
        disp('not on list!');
    end
end


return;
































% sample run -- for plots
%running texture analysis, fixed BW value (cortical only
for i=1:numel(e)
    if ~strcmpi(e(i).name,'sub-NCOV1F006_ses-V01')
        prfix = e(i).name;    
        prfix,

        t1_file = sprintf('texturial_data/%s_t1.nii.gz',prfix);
        mskfile = sprintf('texturial_data/%s_seg.nii.gz',prfix);
        outfile = sprintf('texturial_matfiles_moar/%s_a1',prfix);

        texmap_wrapper3( t1_file, mskfile, [2 41], 0, outfile, [], [], [], [], 'HDE', 0.05,0 ); %[0.0425]
    else
        disp('not on list!');
    end
end
return;

%running texture analysis, optimization (subcortical only)

for i=1:numel(e)
    if ~strcmpi(e(i).name,'sub-NCOV1F006_ses-V01')
        prfix = e(i).name;    
        prfix,

        t1_file = sprintf('texturial_data/%s_t1.nii.gz',prfix);
        mskfile = sprintf('texturial_data/%s_seg.nii.gz',prfix);
        outfile = sprintf('texturial_matfiles_defrange/H_%s',prfix);
        texmap_wrapper3( t1_file, mskfile, [list_s([1 3 5])], 0, outfile, [], [], [], [], 'HDE', 'optim' );
        outfile = sprintf('texturial_matfiles_defrange/KT_%s',prfix);
        texmap_wrapper3( t1_file, mskfile, [list_s([1 3 5])], 0, outfile, [], [], [], [], 'KDE-GT', 'optim' );
        outfile = sprintf('texturial_matfiles_defrange/K_%s',prfix);
        texmap_wrapper3( t1_file, mskfile, [list_s([1 3])], 0, outfile, [], [], [], [], 'KDE-GU', 'optim' );
        return;
    else
        disp('not on list!');
    end
end

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

for i=1:numel(e)
    if ~strcmpi(e(i).name,'sub-NCOV1F006_ses-V01')
        prfix = e(i).name;    
        prfix,

        t1_file = sprintf('texturial_data/%s_t1.nii.gz',prfix);
        mskfile = sprintf('texturial_data/%s_seg.nii.gz',prfix);
        outfile = sprintf('texturial_matfiles_histo/%s',prfix);

        texmap_wrapper3( t1_file, mskfile, [masklab([1 40 80 120 160])], 0, outfile, [], [], [], [], 'HDE', 'optim' );
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
        outfile = sprintf('texturial_matfiles_kerno/%s',prfix);

        texmap_wrapper3( t1_file, mskfile, [masklab([1 40 80 120 160])], 0, outfile, [], [], [], [], 'KDE-G', 'optim' );
    else
        disp('not on list!');
    end
end

