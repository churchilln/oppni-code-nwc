function P1_fmri_prepareData( inputfile, pipelinefile, paramlist, outpath )
%
% =========================================================================
% P1_FMRI_PREPAREDATA: this script should be run after the "P0" step in the bold
% pipeline. It does the following:
% . (a) checks for existence of raw data files
% . (b) ports the raw data over to dir-struct (with basic proc; drop n z-clip)
% . (c) generates QC0 (compatibility checker) + files for QC1 (manual checker)
% =========================================================================
%
% Syntax:
%
%     P1_fmri_prepareData( inputfile, pipelinefile, paramlist, outpath )
%
% Inputs:
%
%      inputfile : string, giving name of input textfile listing data to process
%      pipelinefiles : string, giving name of pipeline textfile specifying which pipeline steps to apply to the data
%      paramlist : string, giving name of parameter file specifying special arguments to use in pipelines
%      outpath : string, specifying destination directory for processed outputs
%

% declaring path
CODE_PATH = fileparts(which('P0_fmri_populateDirectories.m'));
if CODE_PATH(end)~='/'
    CODE_PATH = [CODE_PATH '/'];
end

% initializing structure and running checkes
[subject_list, InputStruct_aug, ~, ~] = P0_fmri_populateDirectories( inputfile, pipelinefile, paramlist, outpath );
% now augmenting outpath... do this after P0!
outpath = fullfile(outpath,'fmri_proc'); % subdir should be fmri_proc
if ~exist(outpath,'dir') error('fmri proc directory does not exist!'); end
e=dir([outpath,'*']); % dir to get absolute
if ~strcmpi( e(1).name, 'fmri_proc') error('first dir should be "fmri_proc"'); end
outpath = fullfile(e(1).folder,e(1).name); % convert to absolute
% now check for expected raw data files
File_Existence_Checker_fmri(InputStruct_aug,outpath,0); clear InputStruct_aug;

%% Copy over raw files and construct qc0/qc1 sheets

% create raw data folders, then copy over func,physio,anat data
for ns=1:numel(subject_list)
    
    % check existence of subject specific struct file
    if exist(fullfile( outpath,subject_list{ns},'InputStruct_ssa.mat'),'file')
        load( fullfile( outpath,subject_list{ns},'InputStruct_ssa.mat') );
    else
        error('cannot find Input struct file for subject: %s \n');
    end

    opath0 = fullfile(outpath,InputStruct_ssa.PREFIX,'rawdata');

    % preparing compatibility qc structure
    QCStruct_compat(ns).PREFIX = InputStruct_ssa.PREFIX;
    QCStruct_compat(ns).N_anat = InputStruct_ssa.N_anat;
    QCStruct_compat(ns).N_func = InputStruct_ssa.N_func;

    % go through anatomical data runs
    for nr=1:InputStruct_ssa.N_anat
        % port over, minimally proc and check the anatomical scan
        if ~exist(sprintf('%s/anat%u.nii',opath0,nr),'file')
            if contains(InputStruct_ssa.arun(nr).ANAT_filename,'.nii.gz')
                unix(sprintf('cp %s %s/anat%u.nii.gz',InputStruct_ssa.arun(nr).ANAT_filename,opath0,nr));
                unix(sprintf('gunzip %s/anat%u.nii.gz',opath0,nr)); % unzipt
            elseif contains(InputStruct_ssa.arun(nr).ANAT_filename,'.nii')
                unix(sprintf('rm %s/anat%u.nii.gz',opath0,nr));
                unix(sprintf('cp %s %s/anat%u.nii',InputStruct_ssa.arun(nr).ANAT_filename,opath0,nr));
            else
                error('unrecognized input file type for:\n\t%s\n',InputStruct_ssa.arun(nr).ANAT_filename);
            end
        end
        % gathering information for qc compatibility stats
        hdr = load_nii_hdr(sprintf('%s/anat%u.nii',opath0,nr));
        % now zippit
        unix(sprintf('gzip %s/anat%u.nii',opath0,nr));
        % orientation info
        orient_RASO = sign([ hdr.hist.srow_x(1) hdr.hist.srow_y(2) hdr.hist.srow_z(3) ]);
        if    ( prod( orient_RASO == [ 1 1 1] )==1 ) orient_RASO(4) =  1; % RAS=radiologic
        elseif( prod( orient_RASO == [-1 1 1] )==1 ) orient_RASO(4) = -1; % LAS=neurologic
        else orient_RASO(4) = 0; % other/unidentified
        end
        % > minimal proc: z-axis clipping
        if ischar(InputStruct_ssa.arun(nr).ZCLIP_thr) && strcmpi(InputStruct_ssa.arun(nr).ZCLIP_thr,'AUTO')
            if ~exist(sprintf('%s/anat%u_zclip.nii.gz',opath0,nr),'file')
                zval = autoclipper( sprintf('%s/anat%u.nii.gz',opath0,nr) );
                unix(sprintf('@clip_volume -input %s/anat%u.nii.gz -below %.02f -prefix %s/anat%u_zclip.nii.gz', ...
                    opath0,nr, [zval],opath0,nr));
            end
        elseif isnumeric(InputStruct_ssa.arun(nr).ZCLIP_thr) && isfinite(InputStruct_ssa.arun(nr).ZCLIP_thr)
            if ~exist(sprintf('%s/anat%u_zclip.nii.gz',opath0,nr),'file')
                unix(sprintf('@clip_volume -input %s/anat%u.nii.gz -below %.02f -prefix %s/anat%u_zclip.nii.gz', ...
                    opath0,nr, [InputStruct_ssa.arun(nr).ZCLIP_thr],opath0,nr));
            end
        else
            disp('skipping z-clipping...');
        end
        %--- store to qc compatibility struct
        QCStruct_compat(ns).arun(nr).voxres_xyzt  = hdr.dime.pixdim(2:5);
        QCStruct_compat(ns).arun(nr).voxdim_xyzt  = hdr.dime.dim(2:5);
        QCStruct_compat(ns).arun(nr).voxdim_t_adj = 1;
        QCStruct_compat(ns).arun(nr).bitpix       = hdr.dime.bitpix;
        QCStruct_compat(ns).arun(nr).orient_RASO  =orient_RASO;
    end

    % go through functional data runs
    for nr=1:InputStruct_ssa.N_func
        % port over, minimally proc and check the functional run
        if ~exist(sprintf('%s/func%u.nii',opath0,nr),'file')
            if contains(InputStruct_ssa.frun(nr).FUNC_filename,'.nii.gz')            
                unix(sprintf('cp %s %s/func%u.nii.gz',InputStruct_ssa.frun(nr).FUNC_filename,opath0,nr));
                unix(sprintf('gunzip %s/func%u.nii.gz',opath0,nr)); % unzipt
            elseif contains(InputStruct_ssa.frun(nr).FUNC_filename,'.nii')
                unix(sprintf('rm %s/func%u.nii.gz',opath0,nr));
                unix(sprintf('cp %s %s/func%u.nii',InputStruct_ssa.frun(nr).FUNC_filename,opath0,nr));
            else
                error('unrecognized input file type for:\n\t%s\n',InputStruct_ssa.frun(nr).FUNC_filename);
            end
        end
        % gathering information for qc compatibility stats
        hdr      = load_nii_hdr(sprintf('%s/func%u.nii',opath0,nr));
        % now zippit
        unix(sprintf('gzip %s/func%u.nii',opath0,nr));
        % time info
        Nt_raw   = hdr.dime.dim(5);
        Nt_adj   = hdr.dime.dim(5) - InputStruct_ssa.frun(nr).DROP_first - InputStruct_ssa.frun(nr).DROP_last;
        % orientation info
        orient_RASO = sign([ hdr.hist.srow_x(1) hdr.hist.srow_y(2) hdr.hist.srow_z(3) ]);
        if    ( prod( orient_RASO == [ 1 1 1] )==1 ) orient_RASO(4) =  1; % RAS=radiologic
        elseif( prod( orient_RASO == [-1 1 1] )==1 ) orient_RASO(4) = -1; % LAS=neurologic
        else orient_RASO(4) = 0; % other/unidentified
        end
        % > minimal proc: scan dropping from runs
        if InputStruct_ssa.frun(nr).DROP_first>0 || InputStruct_ssa.frun(nr).DROP_last>0
            newSTART_afniIdx = InputStruct_ssa.frun(nr).DROP_first; % afni uses zero-rel index
            newEND_afniIdx   = Nt_raw - InputStruct_ssa.frun(nr).DROP_last - 1;  
            if ~exist(sprintf('%s/func%u_drop.nii.gz',opath0,nr),'file')
                unix(sprintf('3dTcat -prefix %s/func%u_drop.nii.gz ''%s/func%u.nii.gz[%u..%u]''',opath0,nr, opath0,nr, newSTART_afniIdx,newEND_afniIdx));
            end
            % making epi mean vol, for later visual qc purposes
            if ~exist( sprintf('%s/func%u_drop_tav.nii.gz',opath0,nr), 'file')
                unix(sprintf('3dTstat -mean -prefix %s/func%u_drop_tav.nii.gz %s/func%u_drop.nii.gz',opath0,nr,opath0,nr));
            end
            if InputStruct_ssa.frun(nr).DROP_first>0 && ~exist(sprintf('%s/func%u_drop_magsat.nii.gz',opath0,nr),'file')
                unix(sprintf('3dcalc -a ''%s/func%u_drop.nii.gz[0]'' -b ''%s/func%u.nii.gz[0]'' -expr ''b-a'' -prefix %s/func%u_drop_magsat.nii.gz',opath0,nr,opath0,nr,opath0,nr));
            end
        else
            disp('skipping drop step...');
            if ~exist( sprintf('%s/func%u_tav.nii.gz',opath0,nr), 'file')
                unix(sprintf('3dTstat -mean -prefix %s/func%u_tav.nii.gz %s/func%u.nii.gz',opath0,nr,opath0,nr));
            end
        end
        %--- store to qc compatibility struct
        QCStruct_compat(ns).frun(nr).voxres_xyzt  = hdr.dime.pixdim(2:5);
        QCStruct_compat(ns).frun(nr).voxdim_xyzt  = hdr.dime.dim(2:5);
        QCStruct_compat(ns).frun(nr).voxdim_t_adj = Nt_adj;
        QCStruct_compat(ns).frun(nr).bitpix       = hdr.dime.bitpix;
        QCStruct_compat(ns).frun(nr).orient_RASO  = orient_RASO;

        % if physio data exists, port over, minimally proc and check the physio run too
        if ~isempty(InputStruct_ssa.frun(nr).PHYSIO_filename)
            unix(sprintf('cp %s.resp.1D %s/physio%u.resp.1D',InputStruct_ssa.frun(nr).PHYSIO_filename,opath0,nr));
            unix(sprintf('cp %s.puls.1D %s/physio%u.puls.1D',InputStruct_ssa.frun(nr).PHYSIO_filename,opath0,nr));
        end
        % physio data dropping to conform to functionale!!
        %%% ===> TO DO!! <=== %%%
        % physio data compatibility checking!!
        %%% ===> TO DO!! <=== %%%
    end
end

disp('exporting compatibility check files.');

% export the .mat structure to folder
save( fullfile(outpath,'_group_level','QC','qc.compat','QCStruct_compat.mat'), 'QCStruct_compat' );

% now go through array and export to arrays, print to file
ytmp = {'ID','res-x','res-y','res-z','res-t','dim-x','dim-y','dim-z','dim-t','dim-t-adj','dtype','R>L','A>P','S>I','RAS'};

% functional first
kq=0; clear catid catmat;
for ns=1:numel(QCStruct_compat)
    for nr=1:QCStruct_compat(ns).N_func
        kq=kq+1;
        xtmp   = QCStruct_compat(ns).frun(nr);
        catid{kq} = strcat(QCStruct_compat(ns).PREFIX,'_run(',num2str(nr),')');
        catmat(kq,:) = [xtmp.voxres_xyzt(:); xtmp.voxdim_xyzt(:); xtmp.voxdim_t_adj; xtmp.bitpix; xtmp.orient_RASO(:)]';
    end
end
catid = pad(catid); ytmp{1}=pad('ID',numel(catid{1})); % pad out the IDs
filo = fopen( fullfile(outpath,'_group_level','QC','qc.compat','table_func_stats.txt'),'w');
fprintf(filo,'%s | %5s  %5s  %5s  %5s  %5s  %5s  %5s  %5s  %9s  %5s  %5s  %5s  %5s  %5s\n',ytmp{:});
for i=1:size(catmat,1)
    fprintf(filo,'%s | %5.2f  %5.2f  %5.2f  %5.2f  %5u  %5u  %5u  %5u  %9u  %5u  %5d  %5d  %5d  %5d\n',catid{i},catmat(i,:));
end
fclose(filo);
%
filo = fopen( fullfile(outpath,'_group_level','QC','qc.compat','outlier_func_stats.txt'),'w');
ko=0;
for i=1:size(catmat,2)
    xmed = median(catmat(:,i));
    ix = find( abs(catmat(:,i)-xmed) > 0.1*abs(xmed) );
    if ~isempty(ix)
        for j = 1:numel(ix)
            ko=ko+1;
            fprintf(filo,'%u. Line=%u/%u, ID=%s, metric=%s:  value is %.2f, median is %.2f\n',ko,ix(j),size(catmat,1),catid{ix(j)},ytmp{i}, catmat(ix(j),i),xmed);
        end
    end
end
fprintf(filo,'\n==> A total of %u functional outlier instances found!\n',ko);
fclose(filo);

% anatomical nirst
kq=0; clear catid catmat;
for ns=1:numel(QCStruct_compat)
    for nr=1:QCStruct_compat(ns).N_anat
        kq=kq+1;
        xtmp   = QCStruct_compat(ns).arun(nr);
        catid{kq} = strcat(QCStruct_compat(ns).PREFIX,'_run(',num2str(nr),')');
        catmat(kq,:) = [xtmp.voxres_xyzt(:); xtmp.voxdim_xyzt(:); xtmp.voxdim_t_adj; xtmp.bitpix; xtmp.orient_RASO(:)]';
    end
end
catid = pad(catid); ytmp{1}=pad('ID',numel(catid{1})); % pad out the IDs
filo = fopen( fullfile(outpath,'_group_level','QC','qc.compat','table_anat_stats.txt'),'w');
fprintf(filo,'%s | %5s  %5s  %5s  %5s  %5s  %5s  %5s  %5s  %9s  %5s  %5s  %5s  %5s  %5s\n',ytmp{:});
for i=1:size(catmat,1)
    fprintf(filo,'%s | %5.2f  %5.2f  %5.2f  %5.2f  %5u  %5u  %5u  %5u  %9u  %5u  %5d  %5d  %5d  %5d\n',catid{i},catmat(i,:));
end
fclose(filo);
%
filo = fopen( fullfile(outpath,'_group_level','QC','qc.compat','outlier_anat_stats.txt'),'w');
ko=0;
for i=1:size(catmat,2)
    xmed = median(catmat(:,i));
    ix = find( abs(catmat(:,i)-xmed) > 0.1*abs(xmed) );
    if ~isempty(ix)
        for j = 1:numel(ix)
            ko=ko+1;
            fprintf(filo,'%u. Line=%u/%u, ID=%s, metric=%s:  value is %.2f, median is %.2f\n',ko,ix(j),size(catmat,1),catid{ix(j)},ytmp{i}, catmat(ix(j),i),xmed);
        end
    end
end
fprintf(filo,'\n==> A total of %u anatomical outlier instances found!\n',ko);
fclose(filo);

disp('Step-1 Complete!');
