function [subject_list] = P2_qc_diffusion( inputfile, outpath, qc_subj_idxes )
%
% . this script: 
% . (a) checks whether input/pipelin/param/template/task/seed files can be read 
% . (b) generates output directory structure, 
% . (c) creates dir-specific InputStructs to protect against overwriting
% . --> should precede any majory script execution

% declaring path
CODE_PATH = fileparts(which('P1_diff_dataProcessing.m'));
if CODE_PATH(end)~='/'
    CODE_PATH = [CODE_PATH '/'];
end

if nargin<3
    qc_subj_idxes=[];
end

%% ========= PHASE ZERO GO ========= %%

outpath = fullfile(outpath,'diff_proc'); % subdir should be diff_proc
if ~exist(outpath,'dir') error('diff proc directory does not exist!'); end
e=dir([outpath,'*']); % dir to get absolute
if ~strcmpi( e(1).name, 'diff_proc') error('first dir should be "diff_proc"'); end
outpath = fullfile(e(1).folder,e(1).name); % convert to absolute

% basic file checks ... construct input/pipeline/param structure -> should throw error if inputs non-valid
InputStruct = Read_Input_File_diff(inputfile);

% now store all files to output
for(ns=1:numel(InputStruct))
    % now reconstruct full folder hierarchy for this subject
    mkdir_r( fullfile( outpath,InputStruct(ns).PREFIX, {'rawdata','diff_proc_p1','diff_proc_p2'}) );
    % store subject ids
    subject_list{ns} = InputStruct(ns).PREFIX;
end

disp('Step-0 Complete!');

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

mkdir_r([outpath,'/_group_level/QC/qc.quant']);
% store information about masking sublist...
subject_list_forqc = subject_list(qc_subj_idxes);
if exist([outpath,'/_group_level/QC/qc.quant/pipe_BASE_qc_subj_idxes.mat'],'file')
    x=load([outpath,'/_group_level/QC/qc.quant/pipe_BASE_qc_subj_idxes.mat']);
    if     ~isempty( setdiff(subject_list_forqc,x.subject_list_forqc) ) 
        error('custom subject list for qcing :: subjects in new list not present in old! delete group level folders if you want to update!')
    elseif ~isempty( setdiff(x.subject_list_forqc,subject_list_forqc) )
        error('custom subject list for qcing :: subjects not in new list that are present in old! delete group level folders if you want to update!')
    else
        disp('custom subject list for qcing :: list is consistent with old one ... continuing without modification!')
    end
else    
    save([outpath,'/_group_level/QC/qc.quant/pipe_BASE_qc_subj_idxes.mat'],'qc_subj_idxes','subject_list_forqc');
end

if ~exist(  fullfile(outpath,'_group_level','QC','qc.quant',['QCStruct_quant_pipe_BASE.mat'])  ,'file')
    
    for ni=1:numel(qc_subj_idxes)

        nr=1;
    
        fprintf('=== PHASE 1, subject %u/%u ===\n',ni,numel(qc_subj_idxes)),
    
        % *** NO SUBJECT SPECIFIC STRUCT FILE YET ***
        % run thru inputstruct and find element matching ni'th ID
        nix=[];
        for ij=1:numel(InputStruct)
            if strcmpi(InputStruct(ij).PREFIX,subject_list_forqc{ni})
                nix = [nix ij];
            end
        end
        if numel(nix)==0 || numel(nix)>2
            error('something went wrong with your indexing')
        end
        InputStruct_ssa = InputStruct(nix); %% single subject
    
        opath0 = fullfile(outpath,InputStruct_ssa.PREFIX,'rawdata');
        opath1 = fullfile(outpath,InputStruct_ssa.PREFIX,'diff_proc_p1/pre_eddy');
        opath2 = fullfile(outpath,InputStruct_ssa.PREFIX,'diff_proc_p1/eddy');
        opath3 = fullfile(outpath,InputStruct_ssa.PREFIX,'diff_proc_p2/dti');
        opath4 = fullfile(outpath,InputStruct_ssa.PREFIX,'diff_proc_p2/noddi');
        opath5 = fullfile(outpath,InputStruct_ssa.PREFIX,'diff_proc_p2/dki');
    
        QCStruct_quant(ni).PREFIX = InputStruct_ssa.PREFIX;
        QCStruct_quant(ni).N_anat = InputStruct_ssa.N_anat;
        QCStruct_quant(ni).N_diff_fwd = InputStruct_ssa.N_diff_fwd;
        QCStruct_quant(ni).N_diff_rev = InputStruct_ssa.N_diff_rev;
    
        mpe = load(sprintf('%s/eddy_unwarp.eddy_parameters',opath2));
        mpe = [mpe(:,1:3) mpe(:,4:6) * (180/pi)];
    
        % absolute [motion] disp
        rmsd=[];
        for i=1:size(mpe,1)-1
            rmsd = [rmsd; mean(  (mpe(i+1:end,:)-mpe(i,:)).^2, 2  ) ];
        end
        QCStruct_quant(ni).drun(nr).mot_avg_tot = mean(rmsd);
        QCStruct_quant(ni).drun(nr).mot_max_tot = max(rmsd);
        % relative [motion] framewise disp
        rmsd = mean( (mpe(2:end,:)-mpe(1:end-1,:)).^2,2);
        QCStruct_quant(ni).drun(nr).mot_avg_rel = mean(rmsd);
        QCStruct_quant(ni).drun(nr).mot_max_rel = max(rmsd);
    
    
        M=load_untouch_niiz(sprintf('%s/refavg_brain_mask.nii.gz',opath1));
        V=load_untouch_niiz(sprintf('%s/eddy_unwarp.eddy_outlier_free_data.nii.gz',opath2));
        Xbval  = load(sprintf('%s/diff_fwd_cat.bval',opath0));
        volmat = nifti_to_mat(V,M);
        % --> OLD VERSION USED TO LUMP ALL VOLUMES INTO RMS ESTIMATE, NOW SPLIT I 
        % --> FOR NOW, OUTLIER CHECKS IMPLICITLY COMBINE ALL DIFF WEIGHTINGS, MAY NEED TO TWEAK LATER! 
        volmat = volmat(:,Xbval>0);

        % dvars [motion] disp
        rmsd=[];
        for i=1:size(volmat,2)-1
            rmsd = [rmsd; mean(  (volmat(:,i+1:end)-volmat(:,i)).^2, 1  )' ];
        end
        QCStruct_quant(ni).drun(nr).dvr_avg_tot = mean(rmsd);
        QCStruct_quant(ni).drun(nr).dvr_max_tot = max(rmsd);
        % dvars [motion] framewise disp
        rmsd = mean( (volmat(:,2:end)-volmat(:,1:end-1)).^2,1);
        QCStruct_quant(ni).drun(nr).dvr_avg_rel = mean(rmsd);
        QCStruct_quant(ni).drun(nr).dvr_max_rel = max(rmsd);
    
    %     % piecewise??
    %     load(sprintf('%s/nvol_fwd.mat',opath0));
    %     Xbval = load(sprintf('%s/diff_fwd_cat.bval',opath0));
    %     Xbvec = load(sprintf('%s/eddy_unwarp.eddy_rotated_bvecs',opath2));
    %     for i=1:numel(nvol_fwd)
    %         ist=sum(nvol_fwd(1:(i-1)))+1;
    %         iln=nvol_fwd(i);
    %         ied=ist+iln-1;
    %         % nifti
    %         unix(sprintf('fslroi %s/eddy_unwarp.eddy_outlier_free_data.nii.gz %s/datachunk.nii.gz %u %u',opath2,opath3,ist-1,iln)); % zero-rel
    %         % bval
    %         dlmwrite(sprintf('%s/datachunk.bval',opath3),Xbval(:,ist:ied),'delimiter',' ')
    %         % bvec
    %         dlmwrite(sprintf('%s/datachunk.bvec',opath3),Xbvec(:,ist:ied),'delimiter',' ','precision',10);
    %         % then delete the "datachunks" + MO map (I don't tend to use it for much...
    %         unix(sprintf('rm %s/datachunk* %s/dtifit_*_MO.nii.gz',opath3,opath3));
    %     end
    
    end

    disp('exporting quantitative qc files.')
    
    mkdir_r(fullfile(outpath,'_group_level','QC','qc.quant'));
    % export the .mat structure to folder
    save( fullfile(outpath,'_group_level','QC','qc.quant',['QCStruct_quant_pipe_BASE.mat']), 'QCStruct_quant' );
else
    disp('qcfile found, reloading...');
    %
    load( fullfile(outpath,'_group_level','QC','qc.quant',['QCStruct_quant_pipe_BASE.mat']));
end


% now go through array and export to arrays, print to file
ytmp = {'mot_avg_tot','mot_avg_rel','mot_max_tot','mot_max_rel','dvr_avg_tot','dvr_avg_rel','dvr_max_tot','dvr_max_rel'};
ytst = {'gam',        'gam',        'gam',        'gam',        'gam',        'gam',        'gam',        'gam',       };
yid  = {'ID'};
% functional first
kq=0; clear catid catmat;
for ni=1:numel(QCStruct_quant)
    for nr=1
        kq=kq+1;
        xtmp   = QCStruct_quant(ni).drun(nr);
        catid{kq} = strcat(QCStruct_quant(ni).PREFIX,'_run(',num2str(nr),')');
        cattmp = [];
        for iu=1:numel(ytmp)
            cattmp = [cattmp xtmp.(ytmp{iu})];
        end
        catmat(kq,:) = [cattmp];
    end
end
catid = pad(catid); yid{1}=pad('ID',numel(catid{1})); % pad out the IDs
filo = fopen( fullfile(outpath,'_group_level','QC','qc.quant',['table_func_stats_pipe_BASE.txt']),'w');
fprintf(filo,'%s | %15s %15s %15s %15s %15s %15s %15s %15s\n',yid{1},ytmp{:});
for i=1:size(catmat,1)
    fprintf(filo,'%s | %15.2f %15.2f %15.2f %15.2f %15.2f %15.2f %15.2f %15.2f\n',catid{i},catmat(i,:));
end
fclose(filo);

out = batch_outlier_testing( catmat, ytst, 0.05, 'FDR' );
out.thr(out.thr<0)=0;

clear abar; legcell = {'mot-mpe','mot-bold'};
abar(:,1) = sum(out.thr(:,[ 1: 4]),2);
abar(:,2) = sum(out.thr(:,[ 5: 7]),2);
figure, bar( abar,'stacked' ); ylim([0 10]);
title('outlier counts - structural data');
legend(legcell);

ix = find( sum(abar,2)>0 );

fprintf('\n\nTotal of %u diffusion runs with outlier values (FDR=0.05):\n',numel(ix)),
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
