function P1_perf_prepareData( inputfile, pipelinefile, paramlist, outpath )
%
% =========================================================================
% P1_PERF_PREPAREDATA: this script should be run after the "P0" step in the perf
% pipeline. It does the following:
% . (a) checks for existence of raw data files
% . (b) ports the raw data over to dir-struct (with basic proc; drop n z-clip)
% . (c) generates QC0 (compatibility checker) + files for QC1 (manual checker)
% =========================================================================
%
% Syntax:
%
%     P1_perf_prepareData( inputfile, pipelinefile, paramlist, outpath )
%
% Inputs:
%
%      inputfile : string, giving name of input textfile listing data to process
%      pipelinefiles : string, giving name of pipeline textfile specifying which pipeline steps to apply to the data
%      paramlist : string, giving name of parameter file specifying special arguments to use in pipelines
%      outpath : string, specifying destination directory for processed outputs
%

% declaring path
CODE_PATH = fileparts(which('P0_perf_populateDirectories.m'));
if CODE_PATH(end)~='/'
    CODE_PATH = [CODE_PATH '/'];
end

% initializing structure and running checkes
[subject_list, InputStruct_aug, ~, ~] = P0_perf_populateDirectories( inputfile, pipelinefile, paramlist, outpath );
% now augmenting outpath... do this after P0!
outpath = fullfile(outpath,'perf_proc'); % subdir should be fmri_proc
% now check for expected raw data files
File_Existence_Checker_perf(InputStruct_aug,outpath,0); clear InputStruct_aug;

%% Copy over raw files and construct qc0/qc1 sheets

% create raw data folders, then copy over perf,physio,anat data
for ns=1:numel(subject_list)
    
    % check existence of subject specific struct file
    if exist(fullfile( outpath,subject_list{ns},'InputStruct_ssa.mat'),'file')
        load( fullfile( outpath,subject_list{ns},'InputStruct_ssa.mat') )
    else
        error('cannot find Input struct file for subject: %s \n');
    end

    opath0 = fullfile(outpath,InputStruct_ssa.PREFIX,'rawdata');

    % preparing compatibility qc structure
    QCStruct_compat(ns).PREFIX = InputStruct_ssa.PREFIX;
    QCStruct_compat(ns).N_anat = InputStruct_ssa.N_anat;
    QCStruct_compat(ns).N_perf = InputStruct_ssa.N_perf;
    QCStruct_compat(ns).N_m0ref = InputStruct_ssa.N_m0ref;

    % go through anatomical data runs
    for nr=1:InputStruct_ssa.N_anat
        % port over, minimally proc and check the anatomical scan
        if ~exist(sprintf('%s/anat%u.nii',opath0,nr),'file')
            if contains(InputStruct_ssa.arun(nr).ANAT_filename,'.nii.gz')
                unix(sprintf('cp %s %s/anat%u.nii.gz',InputStruct_ssa.arun(nr).ANAT_filename,opath0,nr))
                unix(sprintf('gunzip %s/anat%u.nii.gz',opath0,nr)); % unzipt
            elseif contains(InputStruct_ssa.arun(nr).ANAT_filename,'.nii')
                unix(sprintf('cp %s %s/anat%u.nii',InputStruct_ssa.arun(nr).ANAT_filename,opath0,nr))
            else
                error('unrecognized input file type for:\n\t%s\n',InputStruct_ssa.arun(nr).ANAT_filename)
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

    % go through perfusion data runs
    for nr=1:InputStruct_ssa.N_perf
        % port over, minimally proc and check the perfusion run
        if ~exist(sprintf('%s/perf%u.nii',opath0,nr),'file')

            if contains(InputStruct_ssa.prun(nr).PERF_filename,'.nii.gz')
                if strcmpi(InputStruct_ssa.m0loc,'separate')
                    unix(sprintf('cp %s %s/perf%u.nii.gz',InputStruct_ssa.prun(nr).PERF_filename,opath0,nr));
                elseif strcmpi(InputStruct_ssa.m0loc,'infile')
                    unix(sprintf('3dTcat -prefix %s/perf%u.nii.gz ''%s[%u..$]''',opath0,nr, InputStruct_ssa.prun(nr).PERF_filename, (InputStruct_ssa.mrun(nr).M0REF_filename)))
                end
                unix(sprintf('gunzip %s/perf%u.nii.gz',opath0,nr)); % unzipt
            elseif contains(InputStruct_ssa.prun(nr).PERF_filename,'.nii')
                if strcmpi(InputStruct_ssa.m0loc,'separate')
                    unix(sprintf('cp %s %s/perf%u.nii',InputStruct_ssa.prun(nr).PERF_filename,opath0,nr))
                elseif strcmpi(InputStruct_ssa.m0loc,'infile')
                    unix(sprintf('3dTcat -prefix %s/perf%u.nii ''%s[%u..$]''',opath0,nr, InputStruct_ssa.prun(nr).PERF_filename, (InputStruct_ssa.mrun(nr).M0REF_filename)))
                end
            else
                error('unrecognized input file type for:\n\t%s\n',InputStruct_ssa.prun(nr).PERF_filename)
            end

        end
        % gathering information for qc compatibility stats
        hdr      = load_nii_hdr(sprintf('%s/perf%u.nii',opath0,nr));
        % now zippit
        unix(sprintf('gzip %s/perf%u.nii',opath0,nr));
        % time info
        Nt_raw   = hdr.dime.dim(5);
        Nt_adj   = hdr.dime.dim(5) - mod(Nt_raw,2);
        % orientation info
        orient_RASO = sign([ hdr.hist.srow_x(1) hdr.hist.srow_y(2) hdr.hist.srow_z(3) ]);
        if    ( prod( orient_RASO == [ 1 1 1] )==1 ) orient_RASO(4) =  1; % RAS=radiologic
        elseif( prod( orient_RASO == [-1 1 1] )==1 ) orient_RASO(4) = -1; % LAS=neurologic
        else orient_RASO(4) = 0; % other/unidentified
        end
        % > minimal proc: --
        if mod(Nt_raw,2)>0
            if strcmpi( InputStruct_ssa.PWDROP, 'NONE' )
                error('Odd number of Tag/control scans identified -- need to drop one to proceed!')
            else
                if strcmpi( InputStruct_ssa.PWDROP, 'FIRST' )
                    newSTART_afniIdx = 1; % afni uses zero-rel index
                    newEND_afniIdx   = (Nt_raw - 1); % saa, ending on N-1th scan
                elseif strcmpi( InputStruct_ssa.PWDROP, 'LAST' )
                    newSTART_afniIdx = 0; % afni uses zero-rel index
                    newEND_afniIdx   = (Nt_raw - 1) - 1; % saa, ending on N-1th scan
                else
                    error('unrecognized drop scheme! should be NONE, FIRST or LAST!')
                end
                if ~exist(sprintf('%s/perf%u_drop.nii.gz',opath0,nr),'file')
                    unix(sprintf('3dTcat -prefix %s/perf%u_drop.nii.gz ''%s/perf%u.nii.gz[%u..%u]''',opath0,nr, opath0,nr, newSTART_afniIdx,newEND_afniIdx));
                end
                % making epi mean vol, for later visual qc purposes
                if ~exist( sprintf('%s/perf%u_drop_tav.nii.gz',opath0,nr), 'file')
                    unix(sprintf('3dTstat -mean -prefix %s/perf%u_drop_tav.nii.gz %s/perf%u_drop.nii.gz',opath0,nr,opath0,nr));
                end
                % making meandiff vol, to verify label order
                if ~exist( sprintf('%s/perf%u_drop_TC_dif.nii.gz',opath0,nr), 'file')
                    if strcmpi( InputStruct_ssa.LAB_ORDER, 'TC' )
                        unix(sprintf('3dTstat -mean -prefix %s/perf%u_drop_T_tav.nii.gz ''%s/perf%u_drop.nii.gz[0..$(2)]''',opath0,nr,opath0,nr));
                        unix(sprintf('3dTstat -mean -prefix %s/perf%u_drop_C_tav.nii.gz ''%s/perf%u_drop.nii.gz[1..$(2)]''',opath0,nr,opath0,nr));
                    elseif strcmpi( InputStruct_ssa.LAB_ORDER, 'CT' )
                        unix(sprintf('3dTstat -mean -prefix %s/perf%u_drop_T_tav.nii.gz ''%s/perf%u_drop.nii.gz[1..$(2)]''',opath0,nr,opath0,nr));
                        unix(sprintf('3dTstat -mean -prefix %s/perf%u_drop_C_tav.nii.gz ''%s/perf%u_drop.nii.gz[0..$(2)]''',opath0,nr,opath0,nr));
                    else
                        error('unrecognized labelling order, should be TC or CT');
                    end
                    unix(sprintf('3dcalc -a %s/perf%u_drop_T_tav.nii.gz -b %s/perf%u_drop_C_tav.nii.gz -expr ''a-b'' -prefix %s/perf%u_drop_TC_dif.nii.gz',opath0,nr,opath0,nr,opath0,nr));
                    unix(sprintf('rm %s/perf%u_drop_T_tav.nii.gz %s/perf%u_drop_C_tav.nii.gz',opath0,nr,opath0,nr))
                end    
            end
        else
            if strcmpi( InputStruct_ssa.PWDROP, 'FIRST' ) || strcmpi( InputStruct_ssa.PWDROP, 'LAST' )
                error('Even number of Tag/control scans identified -- no need to drop one!')
            elseif strcmpi( InputStruct_ssa.PWDROP, 'NONE' )
                disp('skipping drop step...');
                % making epi mean vol, for later visual qc purposes
                if ~exist( sprintf('%s/perf%u_tav.nii.gz',opath0,nr), 'file')
                    unix(sprintf('3dTstat -mean -prefix %s/perf%u_tav.nii.gz %s/perf%u.nii.gz',opath0,nr,opath0,nr));
                end
                % making meandiff vol, to verify label order
                if ~exist( sprintf('%s/perf%u_TC_dif.nii.gz',opath0,nr), 'file')
                    if strcmpi( InputStruct_ssa.LAB_ORDER, 'TC' )
                        unix(sprintf('3dTstat -mean -prefix %s/perf%u_T_tav.nii.gz ''%s/perf%u.nii.gz[0..$(2)]''',opath0,nr,opath0,nr));
                        unix(sprintf('3dTstat -mean -prefix %s/perf%u_C_tav.nii.gz ''%s/perf%u.nii.gz[1..$(2)]''',opath0,nr,opath0,nr));
                    elseif strcmpi( InputStruct_ssa.LAB_ORDER, 'CT' )
                        unix(sprintf('3dTstat -mean -prefix %s/perf%u_T_tav.nii.gz ''%s/perf%u.nii.gz[1..$(2)]''',opath0,nr,opath0,nr));
                        unix(sprintf('3dTstat -mean -prefix %s/perf%u_C_tav.nii.gz ''%s/perf%u.nii.gz[0..$(2)]''',opath0,nr,opath0,nr));
                    else
                        error('unrecognized labelling order, should be TC or CT');
                    end
                    unix(sprintf('3dcalc -a %s/perf%u_T_tav.nii.gz -b %s/perf%u_C_tav.nii.gz -expr ''a-b'' -prefix %s/perf%u_TC_dif.nii.gz',opath0,nr,opath0,nr,opath0,nr));
                    unix(sprintf('rm %s/perf%u_T_tav.nii.gz %s/perf%u_C_tav.nii.gz',opath0,nr,opath0,nr))
                end     
            else
                error('unrecognized drop scheme! should be NONE, FIRST or LAST!')
            end
        end
        %--- store to qc compatibility struct
        QCStruct_compat(ns).prun(nr).voxres_xyzt  = hdr.dime.pixdim(2:5);
        QCStruct_compat(ns).prun(nr).voxdim_xyzt  = hdr.dime.dim(2:5);
        QCStruct_compat(ns).prun(nr).voxdim_t_adj = Nt_adj;
        QCStruct_compat(ns).prun(nr).bitpix       = hdr.dime.bitpix;
        QCStruct_compat(ns).prun(nr).orient_RASO  = orient_RASO;
    end

    % go through calibration data runs
    catlist = [];
    m0run_idx = [];
    for nr=1:InputStruct_ssa.N_m0ref

        % port over, minimally proc and check the perfusion run
        if ~exist(sprintf('%s/m0ref%u.nii',opath0,nr),'file')

            %--> in case of numeric m0 field, get these from start of perf-run
    
            if strcmpi(InputStruct_ssa.m0loc,'separate')
                if contains(InputStruct_ssa.mrun(nr).M0REF_filename,'.nii.gz')            
                    unix(sprintf('cp %s %s/m0ref%u.nii.gz',InputStruct_ssa.mrun(nr).M0REF_filename,opath0,nr))
                    unix(sprintf('gunzip %s/m0ref%u.nii.gz',opath0,nr)); % unzipt
                elseif contains(InputStruct_ssa.mrun(nr).M0REF_filename,'.nii')
                    unix(sprintf('cp %s %s/m0ref%u.nii',InputStruct_ssa.mrun(nr).M0REF_filename,opath0,nr))
                else
                    error('unrecognized input file type for:\n\t%s\n',InputStruct_ssa.mrun(nr).M0REF_filename)
                end
            elseif strcmpi(InputStruct_ssa.m0loc,'infile')
                if contains(InputStruct_ssa.prun(nr).PERF_filename,'.nii.gz')            
                    unix(sprintf('3dTcat -prefix %s/m0ref%u.nii.gz ''%s[0..%u]''',opath0,nr, InputStruct_ssa.prun(nr).PERF_filename, (InputStruct_ssa.mrun(nr).M0REF_filename-1)))
                    unix(sprintf('gunzip %s/m0ref%u.nii.gz',opath0,nr)); % unzipt
                elseif contains(InputStruct_ssa.prun(nr).PERF_filename,'.nii')
                    unix(sprintf('3dTcat -prefix %s/m0ref%u.nii ''%s[0..%u]''',opath0,nr, InputStruct_ssa.prun(nr).PERF_filename, (InputStruct_ssa.mrun(nr).M0REF_filename-1)))
                else
                    error('unrecognized input file type for:\n\t%s\n',InputStruct_ssa.mrun(nr).M0REF_filename)
                end
            elseif strcmpi(InputStruct_ssa.m0loc,'none') % special case, pulling first "control" scan from each
                if     strcmpi( InputStruct_ssa.LAB_ORDER, 'TC' )
                    refo=2; % first contrl, seocnd vol
                elseif strcmpi( InputStruct_ssa.LAB_ORDER, 'CT' )
                    refo=1; % first control, first ovlume
                end
                if contains(InputStruct_ssa.prun(nr).PERF_filename,'.nii.gz')            
                    unix(sprintf('3dTcat -prefix %s/m0ref%u.nii.gz ''%s[%u]''',opath0,nr, InputStruct_ssa.mrun(nr).PERF_filename, (refo-1)))
                    unix(sprintf('gunzip %s/m0ref%u.nii.gz',opath0,nr)); % unzipt
                elseif contains(InputStruct_ssa.prun(nr).PERF_filename,'.nii')
                    unix(sprintf('3dTcat -prefix %s/m0ref%u.nii ''%s[%u]''',opath0,nr, InputStruct_ssa.mrun(nr).PERF_filename, (refo-1)))
                else
                    error('unrecognized input file type for:\n\t%s\n',InputStruct_ssa.mrun(nr).M0REF_filename)
                end
            end

        end
        % gathering information for qc compatibility stats
        hdr      = load_nii_hdr(sprintf('%s/m0ref%u.nii',opath0,nr));
        % now zippit
        unix(sprintf('gzip %s/m0ref%u.nii',opath0,nr));
        % time info
        Nt_raw   = hdr.dime.dim(5);
        % orientation info
        orient_RASO = sign([ hdr.hist.srow_x(1) hdr.hist.srow_y(2) hdr.hist.srow_z(3) ]);
        if    ( prod( orient_RASO == [ 1 1 1] )==1 ) orient_RASO(4) =  1; % RAS=radiologic
        elseif( prod( orient_RASO == [-1 1 1] )==1 ) orient_RASO(4) = -1; % LAS=neurologic
        else orient_RASO(4) = 0; % other/unidentified
        end
        % > minimal proc: --
        if ~exist( sprintf('%s/m0ref%u_tav.nii.gz',opath0,nr), 'file')
            unix(sprintf('3dTstat -mean -prefix %s/m0ref%u_tav.nii.gz %s/m0ref%u.nii.gz',opath0,nr,opath0,nr));
        end
        %--- store to qc compatibility struct
        QCStruct_compat(ns).prun(nr).voxres_xyzt_m0  = hdr.dime.pixdim(2:5);
        QCStruct_compat(ns).prun(nr).voxdim_xyzt_m0  = hdr.dime.dim(2:5);
        QCStruct_compat(ns).prun(nr).voxdim_t_adj_m0 = Nt_raw;
        QCStruct_compat(ns).prun(nr).bitpix_m0       = hdr.dime.bitpix;
        QCStruct_compat(ns).prun(nr).orient_RASO_m0  = orient_RASO;

        %%% get list for catting...
        catlist = [catlist, sprintf(' %s/m0ref%u.nii.gz',opath0,nr)];
        m0run_idx = [m0run_idx, nr*ones(1,Nt_raw)];
    end
    unix(sprintf('fslmerge -t %s/m0ref_cat.nii.gz%s',opath0,catlist));
    writematrix(m0run_idx,sprintf('%s/m0ref_run_index.txt',opath0));
    unix(sprintf('rm%s',catlist))
end

disp('exporting compatibility check files.')

% export the .mat structure to folder
save( fullfile(outpath,'_group_level','QC','qc.compat','QCStruct_compat.mat'), 'QCStruct_compat' );

% now go through array and export to arrays, print to file
ytmp = {'ID','res-x','res-y','res-z','res-t','dim-x','dim-y','dim-z','dim-t','dim-t-adj','dtype','R>L','A>P','S>I','RAS'};

% perfusion first
kq=0; clear catid catmat;
for ns=1:numel(QCStruct_compat)
    for nr=1:QCStruct_compat(ns).N_perf
        kq=kq+1;
        xtmp   = QCStruct_compat(ns).prun(nr);
        catid{kq} = strcat(QCStruct_compat(ns).PREFIX,'_run(',num2str(nr),')');
        catmat(kq,:) = [xtmp.voxres_xyzt(:); xtmp.voxdim_xyzt(:); xtmp.voxdim_t_adj; xtmp.bitpix; xtmp.orient_RASO(:)]';
    end
end
catid = pad(catid); ytmp{1}=pad('ID',numel(catid{1})); % pad out the IDs
filo = fopen( fullfile(outpath,'_group_level','QC','qc.compat','table_perf_stats.txt'),'w');
fprintf(filo,'%s | %5s  %5s  %5s  %5s  %5s  %5s  %5s  %5s  %9s  %5s  %5s  %5s  %5s  %5s\n',ytmp{:});
for i=1:size(catmat,1)
    fprintf(filo,'%s | %5.2f  %5.2f  %5.2f  %5.2f  %5u  %5u  %5u  %5u  %9u  %5u  %5d  %5d  %5d  %5d\n',catid{i},catmat(i,:));
end
fclose(filo);
%
filo = fopen( fullfile(outpath,'_group_level','QC','qc.compat','outlier_perf_stats.txt'),'w');
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
fprintf(filo,'\n==> A total of %u perfusion outlier instances found!\n',ko);
fclose(filo);

% m0-ref first
kq=0; clear catid catmat;
for ns=1:numel(QCStruct_compat)
    for nr=1:QCStruct_compat(ns).N_m0ref
        kq=kq+1;
        xtmp   = QCStruct_compat(ns).prun(nr);
        catid{kq} = strcat(QCStruct_compat(ns).PREFIX,'_run(',num2str(nr),')');
        catmat(kq,:) = [xtmp.voxres_xyzt_m0(:); xtmp.voxdim_xyzt_m0(:); xtmp.voxdim_t_adj_m0; xtmp.bitpix_m0; xtmp.orient_RASO_m0(:)]';
    end
end
catid = pad(catid); ytmp{1}=pad('ID',numel(catid{1})); % pad out the IDs
filo = fopen( fullfile(outpath,'_group_level','QC','qc.compat','table_m0ref_stats.txt'),'w');
fprintf(filo,'%s | %5s  %5s  %5s  %5s  %5s  %5s  %5s  %5s  %9s  %5s  %5s  %5s  %5s  %5s\n',ytmp{:});
for i=1:size(catmat,1)
    fprintf(filo,'%s | %5.2f  %5.2f  %5.2f  %5.2f  %5u  %5u  %5u  %5u  %9u  %5u  %5d  %5d  %5d  %5d\n',catid{i},catmat(i,:));
end
fclose(filo);
%
filo = fopen( fullfile(outpath,'_group_level','QC','qc.compat','outlier_m0ref_stats.txt'),'w');
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
fprintf(filo,'\n==> A total of %u m0-ref outlier instances found!\n',ko);
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
