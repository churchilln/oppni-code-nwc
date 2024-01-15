close all;
clear;

Xdemo = load('Neurocovid/datasheets_v2/demosheet.txt');    % old: id / group / days-vent / yrs-educ / age / sex / upscor / upcode
                                                           % old: id / ses / 3=days / 4=smys / 5=fever / 6=cough / 7=sorethroat / 8=shortbreath / 9=fatigue / 10=GI / 11=smell+taste/ 12=autre

% ==> scripts for pulling down full optimization runs -- generating summary stats etc,

list_s = [10:13 17:18, 49:54];
list_c = [11101:11141 11143:11175, 12101:12141 12143:12175];
masklab = [list_s list_c];
segl = [1*ones(13,1); 2*ones(148,1)];

% hardcode absolute path
pathloc_fs = 'Neurocovid/NCOV_SUBS';
pathloc_af = 'Neurocovid/fmri_proc';
e = dir(sprintf('%s/sub-*-V01',pathloc_fs));
if isempty(e)
    error('no folders found!');
end

return;
%% HDE (0.048)
kq=0;
for i=1:numel(e)
    if ~strcmpi(e(i).name,'sub-NCOV1F006_ses-V01')
        prfix = e(i).name;    
        prfix,
        kq=kq+1;

        ix = find( Xdemo(:,1) == str2num(prfix(11:13)) );
        group(kq,1) = Xdemo(ix,2);
    
        outfile = sprintf('texturial_matfiles_all/texturial_matfiles_optim_h/%s_Output_Roi_byparc.mat.mat',prfix);
        load(outfile); % (bw x dir x roi)

        for j=1:13 % indexing>direction
            for k=1:160 % indexing>roi
                [vx,ix] = min( out.CVfull(:,j,k) );
                szarr(k,kq,j) = out.SZfull(ix,j,k); % (roi x subj x direction)
                nzarr(k,kq,j) = szarr(k,kq,j) ./ sum( [0.70 0.30]' .* out.SZfull(14:15,j,k) );
                cvarr(k,kq,j) = vx;
                nvarr(k,kq,j) = out.sizecheck(k,6);
                ivarr(k,kq,j) = ix;
            end
        end
        prfix_arr{kq} = prfix;
    end
end
save('texturial_code/bwarr_h.mat','szarr','prfix_arr','group');
return;
%% KDE-GU (0.023)
kq=0;
for i=1:numel(e)
    if ~strcmpi(e(i).name,'sub-NCOV1F006_ses-V01')
        prfix = e(i).name;    
        prfix,
        kq=kq+1;

        ix = find( Xdemo(:,1) == str2num(prfix(11:13)) );
        group(kq,1) = Xdemo(ix,2);
    
        outfile = sprintf('texturial_matfiles_all/texturial_matfiles_optim_u/%s_Output_Roi_byparc.mat.mat',prfix);
        load(outfile); % (bw x dir x roi)

        for j=1:13 % indexing>direction
            for k=1:160 % indexing>roi
                [vx,ix] = min( out.CVfull(:,j,k) );
                szarr(k,kq,j) = out.SZfull(ix,j,k); % (roi x subj x direction)
                nzarr(k,kq,j) = szarr(k,kq,j) ./ sum( [0.70 0.30]' .* out.SZfull(14:15,j,k) );
                cvarr(k,kq,j) = vx;
                nvarr(k,kq,j) = out.sizecheck(k,6);
                ivarr(k,kq,j) = ix;
            end
        end
        prfix_arr{kq} = prfix;
    end
end
save('texturial_code/bwarr_u.mat','szarr','prfix_arr','group');
return;
%% KDE-GT (0.023)
kq=0;
for i=1:numel(e)
    if ~strcmpi(e(i).name,'sub-NCOV1F006_ses-V01')
        prfix = e(i).name;    
        prfix,
        kq=kq+1;

        ix = find( Xdemo(:,1) == str2num(prfix(11:13)) );
        group(kq,1) = Xdemo(ix,2);
    
        outfile = sprintf('texturial_matfiles_all/texturial_matfiles_optim_t/%s_Output_Roi_byparc.mat.mat',prfix);
        load(outfile); % (bw x dir x roi)

        for j=1:13 % indexing>direction
            for k=1:160 % indexing>roi
                [vx,ix] = min( out.CVfull(:,j,k) );
                szarr(k,kq,j) = out.SZfull(ix,j,k); % (roi x subj x direction)
                nzarr(k,kq,j) = szarr(k,kq,j) ./ sum( [0.70 0.30]' .* out.SZfull(14:15,j,k) );
                cvarr(k,kq,j) = vx;
                nvarr(k,kq,j) = out.sizecheck(k,6);
                ivarr(k,kq,j) = ix;
            end
        end
        prfix_arr{kq} = prfix;
    end
end
save('texturial_code/bwarr_t.mat','szarr','prfix_arr','group');

% plotting arrays by size
figure;
t = tiledlayout(4,4,'TileSpacing','Compact','Padding','Compact');
for i=1:13
nexttile
imagesc( szarr(:,:,i),prctile(szarr(:),[5 95]) ); colormap jet;
end


a = reshape(std(szarr,0,3),[],1);
b = reshape(std(szarr,0,2),[],1);
c = reshape(std(szarr,0,1),[],1);

xcat=[];
for i=1:13
    xcat = [xcat, reshape(szarr(:,:,i),[],1)];
end
figure,boxplot(xcat);

xcat=[];
for i=1:71
    xcat = [xcat, reshape(szarr(:,i,:),[],1)];
end
figure,boxplot(xcat);

xcat=[];
for i=1:160
    xcat = [xcat, reshape(szarr(i,:,:),[],1)];
end
figure,bar(median(xcat));









return;
kq=0;
missinge=zeros(numel(e),3)
for i=1:numel(e)
    if ~strcmpi(e(i).name,'sub-NCOV1F006_ses-V01')
        prfix = e(i).name;    
        prfix,
        kq=kq+1;
    
        outfile = sprintf('texturial_matfiles_all/texturial_matfiles_h_aD/%s_aDaS_Output_Roi_byparc.mat',prfix);
        if exist(outfile,'file')
        load(outfile); % (bw x dir x roi)
        mset_av{1}(:,:,kq) = out.metrics_av;
        else
        mset_av{1}(:,:,kq) = NaN*ones(160,6);
        missinge(kq,1) = 1;
        end
        
        outfile = sprintf('texturial_matfiles_all/texturial_matfiles_u_aD/%s_aDaS_Output_Roi_byparc.mat',prfix);
        if exist(outfile,'file')
        load(outfile); % (bw x dir x roi)
        mset_av{2}(:,:,kq) = out.metrics_av;
        else
        mset_av{2}(:,:,kq) = NaN*ones(160,6);
        missinge(kq,2) = 1;
        end

        outfile = sprintf('texturial_matfiles_all/texturial_matfiles_t_aD/%s_aDaS_Output_Roi_byparc.mat',prfix);
        if exist(outfile,'file')
        load(outfile); % (bw x dir x roi)
        mset_av{3}(:,:,kq) = out.metrics_av;
        else
        mset_av{3}(:,:,kq) = NaN*ones(160,6);
        missinge(kq,3) = 1;
        end

    end
end
figure,bar( diag(corr(permute(mset_av{1}(:,4,:),[3 1 2]),permute(mset_av{2}(:,4,:),[3 1 2]),'rows','pairwise')))

return;

szoptim = mean(szarr,3); % roi x subj
% 1.average across directions
kq=0;
for i=1:numel(e)
    if ~strcmpi(e(i).name,'sub-NCOV1F006_ses-V01')
        prfix = e(i).name;    
        prfix,
        kq=kq+1;

        t1_file = sprintf('texturial_data/%s_t1.nii.gz',prfix);
        mskfile = sprintf('texturial_data/%s_seg.nii.gz',prfix);
        outfile = sprintf('texturial_matfiles_aD/%s',prfix);

        texmap_wrapper3( t1_file, mskfile, [masklab], 0, outfile, [], [], [], [], 'KDE-G', szoptim(:,kq) );
    else
        disp('not on list!');
    end
end

szoptim = median(mean(szarr,3),2); % roi x 1
% 2.then average across subjects
kq=0;
for i=1:numel(e)
    if ~strcmpi(e(i).name,'sub-NCOV1F006_ses-V01')
        prfix = e(i).name;    
        prfix,
        kq=kq+1;

        t1_file = sprintf('texturial_data/%s_t1.nii.gz',prfix);
        mskfile = sprintf('texturial_data/%s_seg.nii.gz',prfix);
        outfile = sprintf('texturial_matfiles_aDaS/%s',prfix);

        texmap_wrapper3( t1_file, mskfile, [masklab], 0, outfile, [], [], [], [], 'KDE-G', szoptim(:) );
    else
        disp('not on list!');
    end
end


szoptim = median(median(mean(szarr,3),2));
% 3.then average across rois
kq=0;
for i=1:numel(e)
    if ~strcmpi(e(i).name,'sub-NCOV1F006_ses-V01')
        prfix = e(i).name;    
        prfix,
        kq=kq+1;

        t1_file = sprintf('texturial_data/%s_t1.nii.gz',prfix);
        mskfile = sprintf('texturial_data/%s_seg.nii.gz',prfix);
        outfile = sprintf('texturial_matfiles_aDaSaR/%s',prfix);

        texmap_wrapper3( t1_file, mskfile, [masklab], 0, outfile, [], [], [], [], 'KDE-G', szoptim );
    else
        disp('not on list!');
    end
end

