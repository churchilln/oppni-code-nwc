function [subject_list] = P1_diff_dataProcessing( inputfile, outpath, bigskip, list )
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
    bigskip = 0;
end
if nargin<4
    list=[];
end

PipeStruct.BREF = 'PROX1'; % just uses first-in-run (fwd) and last-in-run (rev)
% PipeStruct.TOPUP   = 'ON'; % apply topup
PipeStruct.TOPUP   = 'OFF';
%% ========= PHASE ZERO GO ========= %%

outpath = fullfile(outpath,'diff_proc'); % subdir should be diff_proc

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

if bigskip==0

%% ========= PHASE ONE GO ========= %%


for ns=1:numel(subject_list)

    fprintf('=== PHASE 1, subject %u/%u ===\n',ns,numel(subject_list)),

    % *** NO SUBJECT SPECIFIC STRUCT FILE YET ***
    InputStruct_ssa = InputStruct(ns); %% single subject

    opath0 = fullfile(outpath,InputStruct_ssa.PREFIX,'rawdata');

    % preparing compatibility qc structure
    QCStruct_compat(ns).PREFIX = InputStruct_ssa.PREFIX;
    QCStruct_compat(ns).N_anat = InputStruct_ssa.N_anat;
    QCStruct_compat(ns).N_diff_fwd = InputStruct_ssa.N_diff_fwd;
    QCStruct_compat(ns).N_diff_rev = InputStruct_ssa.N_diff_rev;

    % preparing unprocced qc structure
    QCStruct_qunproc(ns).PREFIX = InputStruct_ssa.PREFIX;
    QCStruct_qunproc(ns).N_anat = InputStruct_ssa.N_anat;
    QCStruct_qunproc(ns).N_diff_fwd = InputStruct_ssa.N_diff_fwd;
    QCStruct_qunproc(ns).N_diff_rev = InputStruct_ssa.N_diff_rev;
    
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

    % go through diffusion-fwd data runs
    if ~exist(sprintf('%s/diff_fwd_cat.nii.gz',opath0),'file')
        for nr=1:InputStruct_ssa.N_diff_fwd
            % NII port over,
            if ~exist(sprintf('%s/diff_fwd%u.nii',opath0,nr),'file')
                if contains(InputStruct_ssa.dfrun(nr).DIFF_FWD_filename,'.nii.gz')            
                    unix(sprintf('cp %s %s/diff_fwd%u.nii.gz',InputStruct_ssa.dfrun(nr).DIFF_FWD_filename,opath0,nr))
                    unix(sprintf('gunzip %s/diff_fwd%u.nii.gz',opath0,nr)); % unzipt
                elseif contains(InputStruct_ssa.dfrun(nr).DIFF_FWD_filename,'.nii')
                    unix(sprintf('cp %s %s/diff_fwd%u.nii',InputStruct_ssa.dfrun(nr).DIFF_FWD_filename,opath0,nr))
                else
                    error('unrecognized input file type for:\n\t%s\n',InputStruct_ssa.dfrun(nr).DIFF_FWD_filename)
                end
            end
            % BVAL port over,
            if ~exist(sprintf('%s/diff_fwd%u.bval',opath0,nr),'file')
                [inpth,insub,~]=fileparts( InputStruct_ssa.dfrun(nr).DIFF_FWD_filename );
                [~,insub2,~]=fileparts(insub);
                unix(sprintf('cp %s/%s.bval %s/diff_fwd%u.bval',inpth,insub2,opath0,nr))
            end
            % BVEC port over,
            if ~exist(sprintf('%s/diff_fwd%u.bvec',opath0,nr),'file')
                [inpth,insub,~]=fileparts( InputStruct_ssa.dfrun(nr).DIFF_FWD_filename );
                [~,insub2,~]=fileparts(insub);
                unix(sprintf('cp %s/%s.bvec %s/diff_fwd%u.bvec',inpth,insub2,opath0,nr))
            end
    
        end
        % ## and concattify
        if InputStruct_ssa.N_diff_fwd>1
            catlist=[];
            for nr=1:InputStruct_ssa.N_diff_fwd
                catlist = [catlist, sprintf(' %s/diff_fwd%u.nii',opath0,nr)];
                V_hdr=load_nii_hdr(sprintf('%s/diff_fwd%u.nii',opath0,nr));
                nvol_fwd(nr) = V_hdr.dime.dim(5);
            end
            unix(sprintf('fslmerge -t %s/diff_fwd_cat.nii.gz%s',opath0,catlist));
            % concat bval files...
            bvcat=cell(1);
            for nr=1:InputStruct_ssa.N_diff_fwd
                fidbv  = fopen(sprintf('%s/diff_fwd%u.bval',opath0,nr),'rt'); % open file
                tline  = fgetl(fidbv); kq=0; % get 1st line
                while ischar(tline) 
                    kq=kq+1;
                    bvcat{kq} = [bvcat{kq} tline ' ']; % read in lines, concat to cells (resp each. line)
                    tline     = fgetl(fidbv);
                end
                fclose(fidbv);
            end
            fidOUT = fopen(sprintf('%s/diff_fwd_cat.bval',opath0),'wt'); % create bval
            for(i=1:length(bvcat))
                fprintf(fidOUT, '%s\n', bvcat{i});
            end
            fclose(fidOUT);
            % concat bvec files...
            bvcat=cell(3);
            for nr=1:InputStruct_ssa.N_diff_fwd
                fidbv  = fopen(sprintf('%s/diff_fwd%u.bvec',opath0,nr),'rt'); % open file
                tline  = fgetl(fidbv); kq=0; % get 1st line
                while ischar(tline) 
                    kq=kq+1;
                    bvcat{kq} = [bvcat{kq} tline ' ']; % read in lines, concat to cells (resp each. line)
                    tline     = fgetl(fidbv);
                end
                fclose(fidbv);
            end            
            fidOUT = fopen(sprintf('%s/diff_fwd_cat.bvec',opath0),'wt'); % create bval
            for(i=1:length(bvcat))
                fprintf(fidOUT, '%s\n', bvcat{i});
            end
            fclose(fidOUT);
            if (exist(sprintf('%s/diff_fwd_cat.nii',opath0),'file') || exist(sprintf('%s/diff_fwd_cat.nii.gz',opath0),'file')) && exist(sprintf('%s/diff_fwd_cat.bval',opath0),'file') && exist(sprintf('%s/diff_fwd_cat.bvec',opath0),'file')
                for nr=1:InputStruct_ssa.N_diff_fwd
                    unix(sprintf('rm %s/diff_fwd%u.nii %s/diff_fwd%u.bval %s/diff_fwd%u.bvec',opath0,nr,opath0,nr,opath0,nr))
                end
            else
                error('failed to correctly cat forward data!')
            end
        else
            V_hdr=load_nii_hdr(sprintf('%s/diff_fwd1.nii',opath0));
            nvol_fwd(1) = V_hdr.dime.dim(5);
            unix(sprintf('mv %s/diff_fwd1.nii  %s/diff_fwd_cat.nii' ,opath0,opath0))
            unix(sprintf('gzip %s/diff_fwd_cat.nii',opath0))
            unix(sprintf('mv %s/diff_fwd1.bval %s/diff_fwd_cat.bval',opath0,opath0))
            unix(sprintf('mv %s/diff_fwd1.bvec %s/diff_fwd_cat.bvec',opath0,opath0))
        end
        save(sprintf('%s/nvol_fwd.mat',opath0),'nvol_fwd');
    end

    % go through diffusion-rev data runs
    if strcmpi(InputStruct_ssa.REV_MODE,'NONE')
        disp('skipping any reverse encoding stuff');
    else
        if ~exist(sprintf('%s/diff_rev_cat.nii.gz',opath0),'file')
            for nr=1:InputStruct_ssa.N_diff_rev
                % NII port over,
                if ~exist(sprintf('%s/diff_rev%u.nii',opath0,nr),'file')
                    if contains(InputStruct_ssa.drrun(nr).DIFF_REV_filename,'.nii.gz')            
                        unix(sprintf('cp %s %s/diff_rev%u.nii.gz',InputStruct_ssa.drrun(nr).DIFF_REV_filename,opath0,nr))
                        unix(sprintf('gunzip %s/diff_rev%u.nii.gz',opath0,nr)); % unzipt
                    elseif contains(InputStruct_ssa.drrun(nr).DIFF_REV_filename,'.nii')
                        unix(sprintf('cp %s %s/diff_rev%u.nii',InputStruct_ssa.drrun(nr).DIFF_REV_filename,opath0,nr))
                    else
                        error('unrecognized input file type for:\n\t%s\n',InputStruct_ssa.drrun(nr).DIFF_REV_filename)
                    end
                end
                if strcmpi(InputStruct_ssa.REV_MODE,'FULL') %% case where full diffusion encoding data collected
                    % BVAL port over,
                    if ~exist(sprintf('%s/diff_rev%u.bval',opath0,nr),'file')
                        [inpth,insub,~]=fileparts( InputStruct_ssa.drrun(nr).DIFF_REV_filename );
                        [~,insub2,~]=fileparts(insub);
                        unix(sprintf('cp %s/%s.bval %s/diff_rev%u.bval',inpth,insub2,opath0,nr))
                    end
                    % BVEC port over,
                    if ~exist(sprintf('%s/diff_rev%u.bvec',opath0,nr),'file')
                        [inpth,insub,~]=fileparts( InputStruct_ssa.drrun(nr).DIFF_REV_filename );
                        [~,insub2,~]=fileparts(insub);
                        unix(sprintf('cp %s/%s.bvec %s/diff_rev%u.bvec',inpth,insub2,opath0,nr))
                    end
                elseif strcmpi(InputStruct_ssa.REV_MODE,'REF') %% case where only b=0 data collected
                    % BVAL port over,
                    if ~exist(sprintf('%s/diff_rev%u.bval',opath0,nr),'file')
                        V_hdr=load_nii_hdr(sprintf('%s/diff_rev%u.nii',opath0,nr));
                        str = sprintf('%u ',zeros(1,V_hdr.dime.dim(5)));
                        fid=fopen(sprintf('%s/diff_rev%u.bval',opath0,nr),'w');
                        fprintf(fid,'%s\n',str(1:end-1));
                        fclose(fid);
                    end
                    % BVEC port over,
                    if ~exist(sprintf('%s/diff_rev%u.bvec',opath0,nr),'file')
                        V_hdr=load_nii_hdr(sprintf('%s/diff_rev%u.nii',opath0,nr));
                        str = sprintf('%u ',zeros(1,V_hdr.dime.dim(5)));
                        fid=fopen(sprintf('%s/diff_rev%u.bvec',opath0,nr),'w');
                        for i=1:3
                        fprintf(fid,'%s\n',str(1:end-1));
                        end
                        fclose(fid);
                    end
                end
            end
            % ## and concattify
            if InputStruct_ssa.N_diff_rev>1
                catlist=[];
                for nr=1:InputStruct_ssa.N_diff_rev
                    catlist = [catlist, sprintf(' %s/diff_rev%u.nii',opath0,nr)];
                    V_hdr=load_nii_hdr(sprintf('%s/diff_rev%u.nii',opath0,nr));
                    nvol_rev(nr) = V_hdr.dime.dim(5);
                end
                unix(sprintf('fslmerge -t %s/diff_rev_cat.nii.gz%s',opath0,catlist));
                % concat bval files...
                bvcat=cell(1);
                for(nr=1:InputStruct_ssa.N_diff_rev)
                    fidbv  = fopen(sprintf('%s/diff_rev%u.bval',opath0,nr),'rt'); % open file
                    tline  = fgetl(fidbv); kq=0; % get 1st line
                    while ischar(tline) 
                        kq=kq+1;
                        bvcat{kq} = [bvcat{kq} tline ' ']; % read in lines, concat to cells (resp each. line)
                        tline     = fgetl(fidbv);
                    end
                    fclose(fidbv);
                end
                fidOUT = fopen(sprintf('%s/diff_rev_cat.bval',opath0),'wt'); % create bval
                for i=1:length(bvcat)
                    fprintf(fidOUT, '%s\n', bvcat{i});
                end
                fclose(fidOUT);
                % concat bvec files...
                bvcat=cell(3);
                for nr=1:InputStruct_ssa.N_diff_rev
                    fidbv  = fopen(sprintf('%s/diff_rev%u.bvec',opath0,nr),'rt'); % open file
                    tline  = fgetl(fidbv); kq=0; % get 1st line
                    while ischar(tline) 
                        kq=kq+1;
                        bvcat{kq} = [bvcat{kq} tline ' ']; % read in lines, concat to cells (resp each. line)
                        tline     = fgetl(fidbv);
                    end
                    fclose(fidbv);
                end            
                fidOUT = fopen(sprintf('%s/diff_rev_cat.bvec',opath0),'wt'); % create bval
                for i=1:length(bvcat)
                    fprintf(fidOUT, '%s\n', bvcat{i});
                end
                fclose(fidOUT);
                if exist(sprintf('%s/diff_rev_cat.nii',opath0),'file') && exist(sprintf('%s/diff_rev_cat.bval',opath0),'file') && exist(sprintf('%s/diff_rev_cat.bvec',opath0),'file')
                    for nr=1:InputStruct_ssa.N_diff_rev
                        unix(sprintf('rm %s/diff_rev%u.nii %s/diff_rev%u.bval %s/diff_rev%u.bvec',opath0,nr,opath0,nr,opath0,nr));
                    end
                else
                    error('failed to correctly cat forward data!')
                end
            else
                V_hdr=load_nii_hdr(sprintf('%s/diff_rev1.nii',opath0));
                nvol_rev(1) = V_hdr.dime.dim(5);
                unix(sprintf('mv %s/diff_rev1.nii  %s/diff_rev_cat.nii' ,opath0,opath0));
                unix(sprintf('gzip %s/diff_rev_cat.nii',opath0));
                unix(sprintf('mv %s/diff_rev1.bval %s/diff_rev_cat.bval',opath0,opath0));
                unix(sprintf('mv %s/diff_rev1.bvec %s/diff_rev_cat.bvec',opath0,opath0));
            end
            save(sprintf('%s/nvol_rev.mat',opath0),'nvol_rev');
        end
    end
end

%% ========= PHASE TWO GO ========= %%

for ns=1:numel(subject_list)

    fprintf('=== PHASE 2, subject %u/%u ===\n',ns,numel(subject_list)),

    InputStruct_ssa = InputStruct(ns); %% single subject

    opath0 = fullfile(outpath,InputStruct_ssa.PREFIX,'rawdata');
    opath1 = fullfile(outpath,InputStruct_ssa.PREFIX,'diff_proc_p1/pre_eddy');
    opath2 = fullfile(outpath,InputStruct_ssa.PREFIX,'diff_proc_p1/eddy');
    %
    mkdir(opath1);
    mkdir(opath2);

    if ~exist(sprintf('%s/b0ref_fwd.nii.gz',opath1),'file')
        xbfwd = load(sprintf('%s/diff_fwd_cat.bval',opath0));
        ixfwd = find(xbfwd<1E-9);
        % extract first b0 of catt'd run -- this is our reference in general!
        unix(sprintf('fslroi %s/diff_fwd_cat.nii.gz %s/b0ref_fwd.nii.gz %u 1',opath0,opath1,ixfwd(1)-1));
    end

    % phase encoding string
    if     strcmpi(InputStruct_ssa.PE_FWD,'A>>P') pef_str='0 -1 0';
    elseif strcmpi(InputStruct_ssa.PE_FWD,'P>>A') pef_str='0 1 0';
    elseif strcmpi(InputStruct_ssa.PE_FWD,'L>>R') pef_str='0 -1 0';
    elseif strcmpi(InputStruct_ssa.PE_FWD,'R>>L') pef_str='0 1 0';
    else   error('unrecognized PE_FWD option')
    end
    
    % NB "stretched" front = P->A / "squashed" front = A->P
    % rot = 1E-3 * SE * FEPI, where SE = echo spacing in ms, FEPI = epi factor
    % in siemens protocol, find echo spacing in msec (~0.8) / and epi factor (~100), 
    %
    % Effective Echo Spacing (s) = 1/(BandwidthPerPixelPhaseEncode * MatrixSizePhase)
    % Total readout time (FSL) = (number of echoes - 1) * echo spacing
    % Total Readout Time (SPM) = 1/(BandwidthPerPixelPhaseEncode)
    %
    % FEPI = "ReconMatrixPE": 96
    %
    % SE(s) = 1/(96 * 41.667)
    % TR(s) = SE(s) * (96)*SE(s)
    if ~exist(sprintf('%s/refavg_brain.nii.gz',opath1),'file')
        if strcmpi(PipeStruct.TOPUP,'ON')
            
            if ~exist(sprintf('%s/b0ref_rev.nii.gz',opath1),'file')
                xbrev = load(sprintf('%s/diff_rev_cat.bval',opath0));
                ixrev = find(xbrev<1E-9);
                % extract last b0 of catt'd run -- also for topup reference
                unix(sprintf('fslroi %s/diff_rev_cat.nii.gz %s/b0ref_rev.nii.gz %u 1',opath0,opath1,ixrev(end)-1));
            end
    
            % phase encoding string
            if     strcmpi(InputStruct_ssa.PE_REV,'A>>P') per_str='0 -1 0';
            elseif strcmpi(InputStruct_ssa.PE_REV,'P>>A') per_str='0 1 0';
            elseif strcmpi(InputStruct_ssa.PE_REV,'L>>R') per_str='0 -1 0';
            elseif strcmpi(InputStruct_ssa.PE_REV,'R>>L') per_str='0 1 0';
            else   error('unrecognized PE_FWD option')
            end
            % concat b0 files on time ... for topup!
            unix(sprintf('fslmerge -t %s/b0ref_fwd_rev.nii %s/b0ref_fwd.nii %s/b0ref_rev.nii',opath1,opath1,opath1));
            % construct acqfile
            fid=fopen(sprintf('%s/acq_fwd_rev.txt',opath1),'w');
            fprintf(fid,'%s %.05f\n',pef_str,InputStruct_ssa.TRO_MSEC/1000); % A->P is -1 
            fprintf(fid,'%s %.05f\n',per_str,InputStruct_ssa.TRO_MSEC/1000); % P->A is +1, last is total readout time (s)
            fclose(fid);
    
            % doin thems topsups
            if ~exist(sprintf('%s/b0ref_fwd_rev_undist',opath1),'file')
                disp('  running topup...')
                unix(sprintf('topup --imain=%s/b0ref_fwd_rev.nii --datain=%s/acq_fwd_rev.txt --config=b02b0.cnf --out=%s/topup_fwd_rev --iout=%s/b0ref_fwd_rev_undist --fout=%s/topup_fwd_rev_fout',...
                    opath1,opath1,opath1,opath1,opath1)); %12:40
                disp('  ...topup done!');
            end
            
            % now getting mean b0 map from aligned oes 
            unix(sprintf('fslmaths %s/b0ref_fwd_rev_undist -Tmean %s/b0ref_avg_undist',opath1,opath1));
            % and masking out brain tissue
            unix(sprintf('bet %s/b0ref_avg_undist %s/refavg_brain -m -f 0.2',opath1,opath1));
        else
            % construct acqfile
            fid=fopen(sprintf('%s/acq_fwd.txt',opath1),'w');
            fprintf(fid,'%s %.05f\n',pef_str,InputStruct_ssa.TRO_MSEC/1000); % A->P is -1 
            fclose(fid);
            % and masking just using selected b0 volume
            unix(sprintf('bet %s/b0ref_fwd.nii %s/refavg_brain -m -f 0.2',opath1,opath1));
        end
    end

    % and construct "index" file for eddy later ... fwd only for now
    load(sprintf('%s/nvol_fwd.mat',opath0));
    fid = fopen(sprintf('%s/idx_fwd.txt',opath1),'w');
    fprintf(fid,'%u ',1*ones(1,sum(nvol_fwd)));
    fclose(fid);

    if strcmpi(InputStruct_ssa.REV_MODE,'FULL')
        error('full diffusion encoding with 2 PE directions currently unsupported');
        % >> still need to write code that cats'em together for eddy processing, and appropriately runs DTI/DKI/NODDI processing 
    end
    
    % eddy correction 2: repollin
    disp('  running eddy...');
    if ~exist(sprintf('%s/eddy_unwarp.nii.gz',opath2),'file')
        if strcmpi(PipeStruct.TOPUP,'ON')
        unix(sprintf('eddy --imain=%s/diff_fwd_cat.nii --mask=%s/refavg_brain_mask --index=%s/idx_fwd.txt --acqp=%s/acq_fwd_rev.txt --bvecs=%s/diff_fwd_cat.bvec --bvals=%s/diff_fwd_cat.bval --fwhm=0 --topup=%s/topup_fwd_rev --flm=quadratic --out=%s/eddy_unwarp --data_is_shelled --repol',...
            opath0,opath1,opath1,opath1,opath0,opath0,opath1,opath2));
        else
        unix(sprintf('eddy --imain=%s/diff_fwd_cat.nii --mask=%s/refavg_brain_mask --index=%s/idx_fwd.txt --acqp=%s/acq_fwd.txt --bvecs=%s/diff_fwd_cat.bvec --bvals=%s/diff_fwd_cat.bval --fwhm=0 --flm=quadratic --out=%s/eddy_unwarp --data_is_shelled --repol',...
            opath0,opath1,opath1,opath1,opath0,opath0,opath2));
        end
    end
    disp('  ...eddy done!');
end

end % finish bigskip

%% ========= PHASE THREE GO ========= %%

for ns=1:numel(subject_list)

    if( isempty(list) || sum(list==ns)>0 )
    
        fprintf('=== PHASE 3, subject %u/%u ===\n',ns,numel(subject_list)),
    
        InputStruct_ssa = InputStruct(ns); %% single subject
    
        opath0 = fullfile(outpath,InputStruct_ssa.PREFIX,'rawdata');
        opath1 = fullfile(outpath,InputStruct_ssa.PREFIX,'diff_proc_p1/pre_eddy');
        opath2 = fullfile(outpath,InputStruct_ssa.PREFIX,'diff_proc_p1/eddy');
        opath3 = fullfile(outpath,InputStruct_ssa.PREFIX,'diff_proc_p2/dti');
        opath4 = fullfile(outpath,InputStruct_ssa.PREFIX,'diff_proc_p2/noddi');
        opath5 = fullfile(outpath,InputStruct_ssa.PREFIX,'diff_proc_p2/dki');
        %
        mkdir(opath3);
        mkdir(opath4);
        mkdir(opath5);
    
        % dti proc first
        load(sprintf('%s/nvol_fwd.mat',opath0));
        Xbval = load(sprintf('%s/diff_fwd_cat.bval',opath0));
        Xbvec = load(sprintf('%s/eddy_unwarp.eddy_rotated_bvecs',opath2));
        for i=1:numel(nvol_fwd)
            if ~exist(sprintf('%s/dtifit_%u_FA.nii.gz',opath3,i),'file')
                ist=sum(nvol_fwd(1:(i-1)))+1;
                iln=nvol_fwd(i);
                ied=ist+iln-1;
                % nifti
                unix(sprintf('fslroi %s/eddy_unwarp.eddy_outlier_free_data.nii.gz %s/datachunk.nii.gz %u %u',opath2,opath3,ist-1,iln)); % zero-rel
                % bval
                dlmwrite(sprintf('%s/datachunk.bval',opath3),Xbval(:,ist:ied),'delimiter',' ')
                % bvec
                dlmwrite(sprintf('%s/datachunk.bvec',opath3),Xbvec(:,ist:ied),'delimiter',' ','precision',10);
                % now dti fitting...
                unix(sprintf('dtifit -k %s/datachunk.nii.gz -o %s/dtifit_%u -m %s/refavg_brain_mask.nii.gz -r %s/datachunk.bvec -b %s/datachunk.bval',...
                    opath3,opath3,i,opath1,opath3,opath3));
                % then delete the "datachunks" + MO map (I don't tend to use it for much...
                unix(sprintf('rm %s/datachunk* %s/dtifit_*_MO.nii.gz',opath3,opath3));
            else
                disp('skipping dti fitting! already done...')
            end
        end
        if  numel(nvol_fwd)>1 && ~exist([opath4,'/NODDI_fit_odi.nii'],'file') %noddi
            unix(sprintf('cp %s/eddy_unwarp.eddy_outlier_free_data.nii.gz %s/tmpnii.nii.gz',opath2,opath4));
            unix(sprintf('gunzip %s/tmpnii.nii.gz',opath4));
            unix(sprintf('cp %s/refavg_brain_mask.nii.gz %s/tmpmsk.nii.gz',opath1,opath4));
            unix(sprintf('gunzip %s/tmpmsk.nii.gz',opath4));
            % creating .mat ROI set
            CreateROI( sprintf('%s/tmpnii.nii',opath4), sprintf('%s/tmpmsk.nii',opath4), [opath4,'/DTI_Multi_ROI.mat'] );
            % defining the NODDI protocol
            protocol = FSL2Protocol(sprintf('%s/diff_fwd_cat.bval',opath0), sprintf('%s/eddy_unwarp.eddy_rotated_bvecs',opath2));
            % make model
            noddi = MakeModel('WatsonSHStickTortIsoV_B0');
            % batch fitting, all brain voxels
            batch_fitting_single([opath4,'/DTI_Multi_ROI.mat'], protocol, noddi, [opath4,'/paramFit.mat']);
            % store as nifti file
            SaveParamsAsNIfTI([opath4,'/paramFit.mat'],[opath4,'/DTI_Multi_ROI.mat'],sprintf('%s/tmpmsk.nii',opath4),[opath4,'/NODDI_fit']);
        else
            disp('skipping noddi fitting! already done...')
        end
        if numel(nvol_fwd)>1 && ~exist([opath5,'/dtifit_ms_kurt.nii.gz'],'file') % dki
            unix(sprintf('dtifit -k %s/eddy_unwarp.eddy_outlier_free_data.nii.gz -o %s/dtifit_ms -m %s/refavg_brain_mask.nii.gz -r %s/eddy_unwarp.eddy_rotated_bvecs -b %s/diff_fwd_cat.bval --kurt --kurtdir',...
                opath2,opath5,opath1,opath2,opath0));
            unix(sprintf('rm %s/dtifit_*_MO.nii.gz %s/dtifit_*_FA.nii.gz %s/dtifit_*_MD.nii.gz %s/dtifit_*_V*.nii.gz %s/dtifit_*_L*.nii.gz',...
                opath5,opath5,opath5,opath5,opath5));
        else
            disp('skipping dki fitting! already done...')
        end
    
    end
end
