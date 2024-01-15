function P2_fmri_dataProcessing( inputfile, pipelinefile, paramlist, outpath, manual_qc, bigskip )
%
% . this script:
% . does processing!!! f
%

if nargin<6
    bigskip=0;
end

% initializing structure and running checkes
[subject_list, InputStruct_aug, ParamStruct_aug] = P0_fmri_populateDirectories( inputfile, pipelinefile, paramlist, outpath );
% now augmenting outpath... do this after P0!
outpath = fullfile(outpath,'fmri_proc'); % subdir should be fmri_proc
% check for missing files from previous step too
File_Existence_Checker(InputStruct_aug,outpath,1); 

% checking consistency of pipelines across participants -- just a precaution
for ns=1:numel(subject_list)

    % check existence of subject specific struct file
    if exist(fullfile( outpath,subject_list{ns},'InputStruct_ssa.mat'),'file')  
        load( fullfile( outpath,subject_list{ns},'InputStruct_ssa.mat') )
    else
        error('cannot find Input struct file for subject: %s \n');
    end
    if ns==1
        PIPELINE_ref = InputStruct_ssa.PIPELINES;
    else
        if ~isequal( PIPELINE_ref, InputStruct_ssa.PIPELINES)
            error('inconsistent pipelines between %s and %s\n\tCheck your code!\n',subject_list{1}, subject_list{ns});
        end
    end
end
% check validity of analysis model etc. --> carry forward information about the analysis too!
analysis_struct = check_analysis_model( ParamStruct_aug.ANALYSIS );
% also, modify "number of components" field if more than one contrast is specified
if( ~isempty(strfind(ParamStruct_aug.CONTRAST,',' )) )
    analysis_struct.num_comp = 'multi_component';
end
% % check if pipeline steps can be executed --> **currently deactivated, handeded adaptively in pipelines
% Pipe_Compatibility_Checker( InputStruct_aug, PIPELINE_ref, analysis_struct);

%% Checking on man. qc before proceeding...

% check the filetype
if iscell(manual_qc)

    if numel(manual_qc)~=2
        error('if man qc sheets provided, expects two: {anat-qc, func-qc}')
    elseif ~exist(manual_qc{1},'file')
        error('anat qc sheet not found')
    elseif ~exist(manual_qc{2},'file')
        error('func qc sheet not found')
    else
        QCStruct_manual = Read_ManQC_File(inputfile,manual_qc{1},manual_qc{2});
    end
    % export the .mat structure to folder
    save( fullfile(outpath,'_group_level','QC','qc.manual','QCStruct_manual.mat'), 'QCStruct_manual' );

    % now go through array and export to arrays, print to file
    ytmp = {'ID','BCOV','WMGM','SHRP','ALGH','ZSIP','SUSC','RING','INHO','Warn','Fail','Reject?','Notes & Initials'};
    
    % functional first
    kq=0; clear catid catin catno catmat;
    for ns=1:numel(QCStruct_manual)
        for nr=1:QCStruct_manual(ns).N_func
            kq=kq+1;
            xtmp   = QCStruct_manual(ns).frun(nr);
            catid{kq} = strcat(QCStruct_manual(ns).PREFIX,'_run(',num2str(nr),')');
            catin{kq} = QCStruct_manual(ns).INIT_func;
            catno{kq} = xtmp.NOTE;
            catmat(kq,:) = [xtmp.BCOV, xtmp.WMGM xtmp.SHRP xtmp.ALGH xtmp.ZSIP xtmp.SUSC xtmp.SUMWARN xtmp.SUMFAIL xtmp.REJECT]; %xtmp.RING xtmp.INHO
        end
    end
    catid = pad(catid); ytmp{1}=pad('ID',numel(catid{1})); % pad out the IDs
    filo = fopen( fullfile(outpath,'_group_level','QC','qc.manual','table_func_stats.txt'),'w');
    fprintf(filo,'%s | %5s  %5s  %5s  %5s  %5s  %5s  %5s  %5s | %5s %5s %8s   %s\n',ytmp{:});
    for i=1:size(catmat,1)
        fprintf(filo,'%s | %5u  %5u  %5u  %5u  %5u  %5u      -      - | %5u %5u %8u   %s (%s) \n',catid{i},catmat(i,:),catno{i},catin{i});
    end
    fclose(filo);
    %
    filo = fopen( fullfile(outpath,'_group_level','QC','qc.manual','outlier_func_stats.txt'),'w');
    ko=0;
    ix = find( catmat(:,end)>0 );
    if ~isempty(ix)
        for j = 1:numel(ix)
            ko=ko+1;
            fprintf(filo,'%u. Line=%u/%u, ID=%s, rejected:  with %u warns and %u fails\n',ko,ix(j),size(catmat,1),catid{ix(j)}, catmat(ix(j),end-2), catmat(ix(j),end-1) );
        end
    end
    fprintf(filo,'\n==> A total of %u functional scan rejections found!\n',ko);
    fclose(filo);
    
    % anatomical nirst
    kq=0; clear catid catin catno catmat;
    for ns=1:numel(QCStruct_manual)
        for nr=1:QCStruct_manual(ns).N_anat
            kq=kq+1;
            xtmp   = QCStruct_manual(ns).arun(nr);
            catid{kq} = strcat(QCStruct_manual(ns).PREFIX,'_run(',num2str(nr),')');
            catin{kq} = QCStruct_manual(ns).INIT_anat;
            catno{kq} = xtmp.NOTE;
            catmat(kq,:) = [xtmp.BCOV, xtmp.WMGM xtmp.SHRP xtmp.ALGH xtmp.ZSIP xtmp.SUSC xtmp.RING xtmp.INHO xtmp.SUMWARN xtmp.SUMFAIL xtmp.REJECT]; %xtmp.RING xtmp.INHO
        end
    end
    catid = pad(catid); ytmp{1}=pad('ID',numel(catid{1})); % pad out the IDs
    filo = fopen( fullfile(outpath,'_group_level','QC','qc.manual','table_anat_stats.txt'),'w');
    fprintf(filo,'%s | %5s  %5s  %5s  %5s  %5s  %5s  %5s  %5s | %5s %5s %8s   %s\n',ytmp{:});
    for i=1:size(catmat,1)
        fprintf(filo,'%s | %5u  %5u  %5u  %5u  %5u  %5u  %5u  %5u | %5u %5u %8u   %s (%s) \n',catid{i},catmat(i,:),catno{i},catin{i});
    end
    fclose(filo);
    %
    filo = fopen( fullfile(outpath,'_group_level','QC','qc.manual','outlier_anat_stats.txt'),'w');
    ko=0;
    ix = find( catmat(:,end)>0 );
    if ~isempty(ix)
        for j = 1:numel(ix)
            ko=ko+1;
            fprintf(filo,'%u. Line=%u/%u, ID=%s, rejected:  with %u warns and %u fails\n',ko,ix(j),size(catmat,1),catid{ix(j)}, catmat(ix(j),end-2), catmat(ix(j),end-1) );
        end
    end
    fprintf(filo,'\n==> A total of %u anatomical scan rejections found!\n',ko);
    fclose(filo);

    % now export values back to the local directories...
    for ns=1:numel(QCStruct_manual)
        opath0 = fullfile( outpath, QCStruct_manual(ns).PREFIX,'rawdata' );
        file_man = fopen([opath0,'/manual_qc_rating.txt'],'w');

        for nr=1:QCStruct_manual(ns).N_anat
		    catvec = [QCStruct_manual(ns).arun(nr).BCOV QCStruct_manual(ns).arun(nr).WMGM QCStruct_manual(ns).arun(nr).SHRP QCStruct_manual(ns).arun(nr).ALGH QCStruct_manual(ns).arun(nr).ZSIP QCStruct_manual(ns).arun(nr).SUSC QCStruct_manual(ns).arun(nr).RING QCStruct_manual(ns).arun(nr).INHO];
            fprintf(file_man, ['\n===> ANATOMICAL RUN %u\n\n',...
             '\t     Brain coverage (BCOV): %u\n', ...
             '\t     WM/GM contrast (WMGM): %u\n', ...
             '\t Sharpness at edges (SHRP): %u\n', ...
             '\t            Ringing (RING): %u\n', ...
             '\t        Alias/ghost (ALGH): %u\n', ...
             '\t       Zipper/spike (ZISP): %u\n', ...
             '\t      Inhomogeneity (INHO): %u\n', ...
             '\t     Susceptibility (SUSC): %u\n', ...
             '\t      * Misc. notes (NOTE): %s\n\n'], ...
             nr,catvec,QCStruct_manual(ns).arun(nr).NOTE);
        end
        for nr=1:QCStruct_manual(ns).N_func
		    catvec = [QCStruct_manual(ns).frun(nr).BCOV QCStruct_manual(ns).frun(nr).WMGM QCStruct_manual(ns).frun(nr).SHRP QCStruct_manual(ns).frun(nr).ALGH QCStruct_manual(ns).frun(nr).ZSIP QCStruct_manual(ns).frun(nr).SUSC];
            fprintf(file_man, ['\n===> FUNCTIONAL RUN %u\n\n',...
             '\t     Brain coverage (BCOV): %u\n', ...
             '\t     WM/GM contrast (WMGM): %u\n', ...
             '\t Sharpness at edges (SHRP): %u\n', ...
             '\t        Alias/ghost (ALGH): %u\n', ...
             '\t       Zipper/spike (ZISP): %u\n', ...
             '\t     Susceptibility (SUSC): %u\n', ...
             '\t      * Misc. notes (NOTE): %s\n\n'], ...
             nr,catvec,QCStruct_manual(ns).frun(nr).NOTE);
        end
        fprintf(file_man, '\n\n* Reviewer initials (INIT), anatomical scans: %s\n* Reviewer initials (INIT), functional scans: %s\n',QCStruct_manual(ns).INIT_anat,QCStruct_manual(ns).INIT_func);
        fclose(file_man);
    end

elseif strcmpi(manual_qc,'SKIP')
    warning('You chose to skip any manual QCing. Proceed with caution!')
else

    for ns=1:numel(subject_list)
        opath0 = fullfile( outpath, subject_list{ns},'rawdata' );
        QCStruct_manual(ns) =  Read_ManQC_1subj([opath0,'/manual_qc_rating.txt'], N_anat,N_func);
        %% **************** TO BECOMPLETED ***************** %%
    end
end


%% ANAT Processing ...
for ns=1:numel(subject_list)

    % check existence of subject specific struct file
    if exist(fullfile( outpath,subject_list{ns},'InputStruct_ssa.mat'),'file')  
        load( fullfile( outpath,subject_list{ns},'InputStruct_ssa.mat') )
    else
        error('cannot find Input struct file for subject: %s \n');
    end

    fprintf('\n===> anat-proc. now on subj %u/%u: %s...\n',ns,numel(subject_list),subject_list{ns}),
    
    % populating structural alignment dir
    mkdir_r( fullfile(outpath,InputStruct_ssa.PREFIX,'anat_proc',ParamStruct_aug.WarpPrefix{1},'seg' ) );

    nr=1;
    % quick formatting stuff
    opath0 = fullfile( outpath, InputStruct_ssa.PREFIX,'rawdata');
    opath1 = fullfile( outpath, InputStruct_ssa.PREFIX,'anat_proc');
    opath2 = fullfile(outpath,InputStruct_ssa.PREFIX,'anat_proc',ParamStruct_aug.WarpPrefix{1});

    if (ischar(InputStruct_ssa.arun(nr).ZCLIP_thr) && strcmpi(InputStruct_ssa.arun(nr).ZCLIP_thr,'AUTO')) || (isnumeric(InputStruct_ssa.arun(nr).ZCLIP_thr) && isfinite(InputStruct_ssa.arun(nr).ZCLIP_thr))
        zclip_tag = '_zclip';
    else
        zclip_tag = '';
    end

    if ~exist(sprintf('%s/anat%u%s_deob.nii',opath0,nr,zclip_tag),'file')
        unix(sprintf('3dWarp -oblique2card -prefix %s/anat%u_deob.nii -wsinc5 %s/anat%u%s.nii', opath1,nr, opath0,nr, zclip_tag));
    end
    if ~exist( sprintf('%s/anat%u_2std.nii', opath1,nr), 'file')
        unix(sprintf('fslreorient2std %s/anat%u_deob.nii %s/anat%u_2std.nii', opath1,nr, opath1,nr));        
        unix(sprintf('gunzip %s/anat%u_2std.nii.gz', opath1,nr));
    end

    % construct the mask next, for skull stripping in alignment

    % actual alignment part...
    if strcmpi(InputStruct_ssa.PIPELINES.WARP{1},'AFNI')

        if ~exist( [opath2,'/anatQQ.x.nii'], 'file' )
            %-- alignment to template for structural
            sswarper_script( sprintf('%s/anat%u_2std.nii', opath1,nr), 'x', ParamStruct_aug.TEMPLATE_loc, opath2);
            %unix(sprintf('@SSwarper -input %s/anat%u_2std.nii -base %s -subid x -odir %s', opath1,nr, ParamStruct_aug.TEMPLATE_loc, opath2));
        else
            disp('warp found! skipping to next step...')
        end
        % create renames copy of warped image for compatibility w other stuff
        if ~exist(sprintf('%s/anat_warped.nii',opath2),'file')
            unix( sprintf('cp %s/anatQQ.x.nii %s/anat_warped.nii',opath2,opath2) );
        end
        % copy "clean" masked file into seg folder for tissue segmentation
        if ~exist(sprintf('%s/seg/anat_stripref.nii',opath2),'file')
            unix( sprintf('cp %s/anatSS.x.nii %s/seg/anat_stripref.nii',opath2,opath2) );
        end

    elseif strcmpi(InputStruct_ssa.PIPELINES.WARP{1},'AFNI1') %% --sswarper + epimasking

        if ~exist( [opath2,'/amskBrainExtractionBrain.nii.gz'], 'file' )
            unix(sprintf('antsBrainExtraction.sh -d 3 -a %s/anat%u_2std.nii -e MNI152_2009_template_SSW_0.nii.gz -m MNI152_2009_template_SSW_mask.nii.gz -o %s/amsk -k 1',...
                opath1,nr,opath2));
        end

    elseif strcmpi(InputStruct_ssa.PIPELINES.WARP{1},'ANTS')
        %%% *** Need fixing for compatibility with param struct
        % creating bias-corrected, skull-stripped anats
        unix(sprintf('antsBrainExtraction.sh -d 3 -a %s/anat%u_2std.nii -e MNI152_2009_template_SSW_0.nii -m MNI152_2009_template_SSW_mask.nii -o %s/amsk -k 1',opath1,nr, opath2));
        % params from "newAntsExample.sh"+
        fixd = 'MNI152_2009_template_SSW_brain.nii';
        movn = sprintf('%s/amskBrainExtractionBrain.nii.gz',opath2);
        nomen = sprintf('%s/aln',opath2);
        unix(sprintf(['antsRegistration ', ...
        '-d 3 ', ...
        '-r [ %s, %s , 1 ] ', ...
        '-m mattes[  %s, %s , 1 , 32, regular, 0.3 ] ', ...
        '-t translation[ 0.1 ] ', ...
        '-c [ 10000x111110x11110,1.e-8,20 ] ', ...
        '-s 4x2x1vox  ', ...
        '-f 6x4x2 -l 1 ', ...
        '-m mattes[  %s, %s , 1 , 32, regular, 0.3 ] ', ...
        '-t rigid[ 0.1 ] ', ...
        '-c [ 10000x111110x11110,1.e-8,20 ] ', ...
        '-s 4x2x1vox  ', ...
        '-f 3x2x1 -l 1 ', ...
        '-m mattes[  %s, %s , 1 , 32, regular, 0.3 ] ', ...
        '-t affine[ 0.1 ] ', ...
        '-c [ 10000x111110x11110,1.e-8,20 ] ', ...
        '-s 4x2x1vox  ', ...
        '-f 3x2x1 -l 1 ', ...
        '-m mattes[  %s, %s , 0.5 , 32 ] ', ...
        '-m cc[  %s, %s , 0.5 , 4 ] ', ...
        '-t SyN[ .20, 3, 0 ] ', ...
        '-c [ 100x100x50,-0.01,5 ] ', ...
        '-s 1x0.5x0vox ', ...
        '-f 4x2x1 -l 1 -u 1 -z 1 ', ...
        '-o [ %s,%s_diff.nii.gz,%s_inv.nii.gz]'],...
        fixd, movn, fixd, movn, fixd, movn, fixd, movn, fixd, movn, fixd, movn, nomen, nomen, nomen ));
        %
        unix(sprintf('antsApplyTransforms -d 3 -i %s -r %s -n linear -t %s1Warp.nii.gz -t %s0GenericAffine.mat -o %s_warped.nii.gz',movn,fixd,nomen,nomen,nomen));
    end

    disp('t1 mask n warp done! Now for some tissue segmentation...');
    
    %-- get segmentation
    if ~exist(sprintf('%s/seg/anat_seg_GM.nii',opath2),'file')
        unix(sprintf('fast -R 0.3 -H 0.1 -t 1 %s/seg/anat_stripref.nii',opath2));
        unix(sprintf('gunzip %s/seg/anat_stripref_pve_*.nii.gz',opath2));
        unix(sprintf('rm %s/seg/anat_stripref_mixeltype.nii.gz',opath2));
        unix(sprintf('rm %s/seg/anat_stripref_pveseg.nii.gz',opath2));
        unix(sprintf('rm %s/seg/anat_stripref_seg.nii.gz',opath2));
        % renaming 'em
        VX = load_untouch_nii(sprintf('%s/seg/anat_stripref.nii',opath2));
        for i=1:3
            VS = load_untouch_nii(sprintf('%s/seg/anat_stripref_pve_%u.nii',opath2,i-1));
            sval(i,1) = mean(double(VX.img(VS.img>0.5)));
        end
        isort = sortrows([(1:3)',sval],2,'ascend');
        isort = isort(:,1);           % pves indexed by increasing mean intensity
        tisslist = {'CSF','GM','WM'}; % tissues in increasing order of T1 intensity
        for i=1:3
            unix(sprintf('mv %s/seg/anat_stripref_pve_%u.nii %s/seg/anat_seg_%s.nii',opath2, (isort(i)-1), opath2, tisslist{i}));
        end
    end
    
    % warping segmentations
    tisslist = {'CSF','GM','WM'}; % tissues in increasing order of T1 intensity
    for i=1:3
        if ~exist( sprintf('%s/seg/anat_seg_%s_warped.nii',opath2,tisslist{i}),'file')
            unix(sprintf('3dNwarpApply -master %s/anatQQ.x.nii -source %s/seg/anat_seg_%s.nii -nwarp "%s/anatQQ.x_WARP.nii %s/anatQQ.x.aff12.1D" -prefix %s.x.nlin',opath2,opath2,tisslist{i},opath2,opath2,tisslist{i}));
            unix(sprintf('3dAFNItoNIFTI %s.x.nlin+tlrc -prefix %s.x.nlin.nii',tisslist{i},tisslist{i}));
            unix(sprintf('rm %s.x.nlin+tlrc.BRIK %s.x.nlin+tlrc.BRIK.gz %s.x.nlin+tlrc.HEAD',tisslist{i},tisslist{i},tisslist{i}));
            unix(sprintf('mv %s.x.nlin.nii %s/seg/anat_seg_%s_warped.nii',tisslist{i},opath2,tisslist{i}));
        else
            disp('skipping tissue seg warping...')
        end
    end

    % create individual masks -- for later qc
    if ~exist(sprintf('%s/seg/anat_mask.nii',opath2),'file') || ~exist(sprintf('%s/anat_warped_mask.nii',opath2),'file')
        unix(sprintf('3dmask_tool -dilate_input 5 -5 -fill_holes -input %s/seg/anat_stripref.nii -prefix %s/seg/anat_mask.nii',opath2,opath2))
        unix(sprintf('3dmask_tool -dilate_input 5 -5 -fill_holes -input %s/anat_warped.nii -prefix %s/anat_warped_mask.nii',opath2,opath2))
    end
end

return;

%% ANATOMICAL INTERMEZZO

% warp-specific path
catpath = [ParamStruct_aug.WarpPrefix{1},'.',ParamStruct_aug.VOXRES,'mm'];

if ~exist( [outpath,'/_group_level/masks/',catpath,'/anat_brain_mask_grp.nii'], 'file') % Brain mask

    for ns=1:numel(subject_list)
    
        % check existence of subject specific struct file
        if exist(fullfile( outpath,subject_list{ns},'InputStruct_ssa.mat'),'file')  
            load( fullfile( outpath,subject_list{ns},'InputStruct_ssa.mat') )
        else
            error('cannot find Input struct file for subject: %s \n');
        end
        opath2 = fullfile(outpath,InputStruct_ssa.PREFIX,'anat_proc',ParamStruct_aug.WarpPrefix{1});
        Vw = load_untouch_nii(sprintf('%s/anat_warped_mask.nii',opath2));
        if ns==1
            volref = double(Vw.img);
        else
            volref = volref + double(Vw.img);
        end
    end
    volref = volref./numel(subject_list);
    Vw.img = double(volref>0.50); % included in majority of individuals
    save_untouch_nii(Vw,[outpath,'/_group_level/masks/',catpath,'/anat_brain_mask_grp.nii']);
end

if ~exist( [outpath,'/_group_level/brain_maps/',catpath,'/anat_brain_grp.nii'], 'file') % Mean image

    for ns=1:numel(subject_list)
    
        % check existence of subject specific struct file
        if exist(fullfile( outpath,subject_list{ns},'InputStruct_ssa.mat'),'file')  
            load( fullfile( outpath,subject_list{ns},'InputStruct_ssa.mat') )
        else
            error('cannot find Input struct file for subject: %s \n');
        end
        opath2 = fullfile(outpath,InputStruct_ssa.PREFIX,'anat_proc',ParamStruct_aug.WarpPrefix{1});
        Vw = load_untouch_nii(sprintf('%s/anat_warped.nii',opath2));
        if ns==1
            volref = double(Vw.img);
        else
            volref = volref + double(Vw.img);
        end
    end
    Vw.img = volref./numel(subject_list);
    save_untouch_nii(Vw,[outpath,'/_group_level/brain_maps/',catpath,'/anat_brain_grp.nii']);
end

tisslist = {'CSF','GM','WM'};
for i=1:3
    if ~exist( [outpath,'/_group_level/brain_maps/',catpath,'/anat_',tisslist{i},'_grp.nii'], 'file') % Mean image
    
        for ns=1:numel(subject_list)
        
            % check existence of subject specific struct file
            if exist(fullfile( outpath,subject_list{ns},'InputStruct_ssa.mat'),'file')  
                load( fullfile( outpath,subject_list{ns},'InputStruct_ssa.mat') )
            else
                error('cannot find Input struct file for subject: %s \n');
            end
            opath2 = fullfile(outpath,InputStruct_ssa.PREFIX,'anat_proc',ParamStruct_aug.WarpPrefix{1});
            Vw = load_untouch_nii(sprintf('%s/seg/anat_seg_%s_warped.nii',opath2,tisslist{i}));
            if ns==1
                volref = double(Vw.img);
            else
                volref = volref + double(Vw.img);
            end
        end
        Vw.img = volref./numel(subject_list);
        save_untouch_nii(Vw,[outpath,'/_group_level/brain_maps/',catpath,'/anat_',tisslist{i},'_grp.nii']);
    end
end

% rough "ventricular/sulcal map" for alignment checking
Va = load_untouch_nii([outpath,'/_group_level/brain_maps/',catpath,'/anat_brain_grp.nii']);
Vb = load_untouch_nii([outpath,'/_group_level/masks/',catpath,'/anat_brain_mask_grp.nii']);
vsm = double(Va.img);
vsi = 1 - (vsm - min(vsm(:)))./(max(vsm(:))-min(vsm(:)));
vsi  = vsi .* double( Vb.img );
V.img = double( vsi > 0.5);% prctile(vsi(vsi>0),75));
save_untouch_nii(V,[outpath,'/_group_level/brain_maps/',catpath,'/anat_sulcal_grp.nii'])


%% PHYSIO PROCESSING
for ns=1:numel(subject_list)

    % check existence of subject specific struct file
    if exist(fullfile( outpath,subject_list{ns},'InputStruct_ssa.mat'),'file') 
        load( fullfile( outpath,subject_list{ns},'InputStruct_ssa.mat') )
    else
        error('cannot find Input struct file for subject: %s \n');
    end

    fprintf('\n===> phys-proc. now on subj %u/%u: %s...\n',ns,numel(subject_list),subject_list{ns}),
    
    disp('nothin for physio-proc so far. put in soon!')
end

%% FUNC PROCESSING...BLOCK1...
for ns=1:numel(subject_list)

    % check existence of subject specific struct file
    if exist(fullfile( outpath,subject_list{ns},'InputStruct_ssa.mat'),'file')
        load( fullfile( outpath,subject_list{ns},'InputStruct_ssa.mat') )
    else
        error('cannot find Input struct file for subject: %s \n');
    end

    fprintf('\n===> func-proc. now on subj %u/%u: %s...\n',ns,numel(subject_list),subject_list{ns}),

    %%%%%========== compatibility adjustments for BLOCK1 ... RICOR only
    missing_physio=0;
    for nr=1:InputStruct_ssa.N_func
        if isempty(InputStruct_ssa.frun(nr).PHYSIO_filename)
            missing_physio=missing_physio+1;
        end
    end
    if missing_physio>0 && ~(strcmpi(InputStruct_ssa.PIPELINES.RICOR{1},'0') || strcmpi(InputStruct_ssa.PIPELINES.RICOR{1},'OFF'))
        warning('subject %s has %u/%u runs with missing physio. Turning RICOR off!',InputStruct_ssa.PREFIX, missing_physio,InputStruct_ssa.N_func)
        InputStruct_ssa.PIPELINES.RICOR{1}='0';
    end
    %%%%%========== compatibility adjustments, done

    % populating functional alignment dir
    mkdir_r( fullfile(outpath,InputStruct_ssa.PREFIX,'func_proc_p1',[ParamStruct_aug.WarpPrefix{1},'.',ParamStruct_aug.VOXRES,'mm'],'seg' ) );
    mkdir_r( fullfile(outpath,InputStruct_ssa.PREFIX,'func_proc_p2',[ParamStruct_aug.WarpPrefix{1},'.',ParamStruct_aug.VOXRES,'mm'] ) );

    % quick formatting stuff
    opath0 = fullfile(outpath,InputStruct_ssa.PREFIX,'rawdata');
    opath1 = fullfile(outpath,InputStruct_ssa.PREFIX,'func_proc_p1');
    opath2 = fullfile(outpath,InputStruct_ssa.PREFIX,'phys_proc');
    opath3a= fullfile(outpath,InputStruct_ssa.PREFIX,'anat_proc',ParamStruct_aug.WarpPrefix{1});
    opath3f= fullfile(outpath,InputStruct_ssa.PREFIX,'func_proc_p1',[ParamStruct_aug.WarpPrefix{1},'.',ParamStruct_aug.VOXRES,'mm']);
    opath4f= fullfile(outpath,InputStruct_ssa.PREFIX,'func_proc_p2',[ParamStruct_aug.WarpPrefix{1},'.',ParamStruct_aug.VOXRES,'mm']);
    
    for nr=1:InputStruct_ssa.N_func

        if InputStruct_ssa.frun(nr).DROP_first>0 || InputStruct_ssa.frun(nr).DROP_last>0
            drop_tag = '_drop';
        else
            drop_tag = '';
        end

    %%  preliminary estimates of motion/motion-like displacement...

        if ~exist(sprintf('%s/func%u_motref.txt',opath1,nr),'file') || ~exist(sprintf('%s/func%u_outlier_dat.mat',opath1,nr),'file')
    
            if ~exist(sprintf('%s/init_mot_estim/func%u_smo6_mask.nii',opath1,nr),'file')
                % preparing data for first estimates of displacement...smooth 'n' mask
                unix(sprintf('3dmerge -prefix %s/init_mot_estim/func%u_smo6.nii -doall -1blur_fwhm 6 %s/func%u%s.nii',opath1,nr,opath0,nr,drop_tag));
                unix(sprintf('3dAutomask -prefix %s/init_mot_estim/func%u_smo6_mask.nii %s/init_mot_estim/func%u_smo6.nii',opath1,nr,opath1,nr));
            else
                disp('skipping dataprep...')
            end
    
            %--- a. mindisp estimation for alignment... (we can keep this, since it probably won't change after slice-proc)
    
            VF     = load_untouch_nii( sprintf('%s/init_mot_estim/func%u_smo6.nii',opath1,nr) );
            MF     = load_untouch_nii( sprintf('%s/init_mot_estim/func%u_smo6_mask.nii',opath1,nr) );
            % convert to 4D fMRI data volume (VV.img) into 2D matrix (voxels x time) for analysis
            epimat = nifti_to_mat( VF,MF ); 
            epimat = bsxfun(@minus,epimat,mean(epimat,2)); % mean-centering
            [u,l,v]=svd( epimat,'econ' );
            Qdat = (v(:,1:end-1)*l(1:end-1,1:end-1))';
            Qmed = median( Qdat,2 );
            Dist = sqrt(sum(bsxfun(@minus,Qdat,Qmed).^2));
            % special cost-function vector (optimal is minimum)
            dmat(:,1) = (Dist-min(Dist))./(max(Dist)-min(Dist)); % standardized distance from medioid [0,1]
            dmat(:,2) = double( ([1; dmat(1:end-1,1)] > 0.5) | ([1;1; dmat(1:end-2,1)] > 0.5) | ([dmat(2:end,1) ;1] > 0.5) | ([dmat(3:end,1) ;1;1] > 0.5) ); % +1 if possible neighbouring spikes +- 2TR
            dmat(:,3) = zeros(numel(Dist),1); dmat(1:3,3)=1; dmat(end-2:end,3)=1; % +1 if start or end of the run
            [v imed1] = min( sum(dmat,2) ); % minimizer of 1, subject to costs 2 & 3
            %
            imed_afni = imed1-1; % zero-relative indexing
            dlmwrite( sprintf('%s/func%u_motref.txt',opath1,nr), [imed_afni] ); % ** found min-disp brick
    
            if ~exist(sprintf('%s/init_mot_estim/func%u_mc_smo6_mask_dil.nii',opath1,nr),'file')
                %--- b. displacement estimation, for qc and despiking purposes only, since we will be processing a bit before MC'ing afterwards...
                % motion correction for gettting mpes, and then smooth, then re-get (better?) brain mask
                unix(sprintf('3dvolreg -prefix %s/init_mot_estim/func%u_mc.nii -1Dfile %s/init_mot_estim/func%u_mpe -maxdisp1D %s/init_mot_estim/func%u_maxdisp -base %u %s/func%u%s.nii',...
                    opath1,nr, opath1,nr, opath1,nr, imed_afni,opath0,nr,drop_tag));
                unix(sprintf('3dmerge -prefix %s/init_mot_estim/func%u_mc_smo6.nii -doall -1blur_fwhm 6 %s/init_mot_estim/func%u_mc.nii',...
                    opath1,nr,opath1,nr));
                unix(sprintf('3dAutomask -prefix %s/init_mot_estim/func%u_mc_smo6_mask.nii %s/init_mot_estim/func%u_mc_smo6.nii',...
                    opath1,nr,opath1,nr)); % ** gets a decent functional mask
                % dilate the new brain mask for identifying outliers in despiking
                unix(sprintf('fslmaths %s/init_mot_estim/func%u_mc_smo6_mask.nii -dilD %s/init_mot_estim/func%u_mc_smo6_mask_dil.nii',...
                    opath1,nr,opath1,nr));
                unix(sprintf('gunzip %s/init_mot_estim/func%u_mc_smo6_mask_dil.nii.gz',...
                    opath1,nr));
            else
                disp('skipping initial disp estimation...')
            end
    
            %-- a. collecting outlier estimates
    
            % done on non-motcorred data; want to find the *really* extreme bad cases
            outlier_dat = spike_estimator( sprintf('%s/init_mot_estim/func%u_smo6.nii',opath1,nr), sprintf('%s/init_mot_estim/func%u_mc_smo6_mask_dil.nii',opath1,nr), sprintf('%s/init_mot_estim/func%u_mpe',opath1,nr), ['--unused--'], 7,1 );
            save(sprintf('%s/func%u_outlier_dat.mat',opath1,nr),'outlier_dat');

            unix(sprintf('rm %s/init_mot_estim/*.nii',opath1)); % clear nifti data to save on space -- keeps only motion disp data
        else
            imed_afni = load( sprintf('%s/func%u_motref.txt',opath1,nr) );
            load(sprintf('%s/func%u_outlier_dat.mat',opath1,nr)); % load "outlier_dat" structure
        end

    %%  slice-based processing...

        if ~exist( sprintf('%s/func%u_despike.nii',opath1,nr), 'file')
            %-- b. trimmming the outliers
            % replacing volumes that are both extreme BOLD values and correspond to high-estimated-motion timepoints
            fmri_interpolator( sprintf('%s/func%u%s.nii',opath0,nr,drop_tag), outlier_dat, InputStruct_ssa.PIPELINES.DESPIKE{1}, sprintf('%s/func%u_despike.nii',opath1,nr) )
        else
            disp('skipping despike interpolation...')
        end

        if ~exist( sprintf('%s/func%u_ricor.nii',opath1,nr), 'file')
            if strcmpi(InputStruct_ssa.PIPELINES.RICOR{1},'0') || strcmpi(InputStruct_ssa.PIPELINES.RICOR{1},'OFF')
                unix(sprintf('cp %s/func%u_despike.nii %s/func%u_ricor.nii',opath1,nr, opath1,nr));
            elseif strcmpi(InputStruct_ssa.PIPELINES.RICOR{1},'1') || strcmpi(InputStruct_ssa.PIPELINES.RICOR{1},'ON')
                %-- c. slice-based physio correction
                %
                % get the stripped-down physio data matrix
                fid = fopen( sprintf('%s.slibase.1D',ostr4) );
                tline = fgetl(fid);
                kq=0; ise=0;
                while ischar(tline) 
                    if    ( contains(tline,'# >') ) ise=1;
                    elseif( contains(tline,'# <') ) ise=2;
                    elseif( ise==1 ) % only if currently flagged on
                        kq=kq+1;
                        ricormat(kq,:) = str2num( tline );
                    end
                    tline = fgetl(fid);
                end
                fclose(fid);
                % take the matrix of slicewise physio covariates and regress slice-by-slice from the despiked data
                ricor_regress( sprintf('%s/func%u_despike.nii',opath1,nr), ricormat, sprintf('%s/func%u_ricor.nii',opath1,nr) );
            end
        else
            disp('skipping ricor regression...')
        end

        if ~exist(sprintf('%s/func%u_tshift.nii',opath1,nr),'file')
            if strcmpi(InputStruct_ssa.PIPELINES.TSHIFT{1},'0') || strcmpi(InputStruct_ssa.PIPELINES.TSHIFT{1},'OFF')
                unix(sprintf('cp %s/func%u_ricor.nii %s/func%u_tshift.nii',opath1,nr, opath1,nr));
            elseif strcmpi(InputStruct_ssa.PIPELINES.TSHIFT{1},'1') || strcmpi(InputStruct_ssa.PIPELINES.TSHIFT{1},'ON')
                %-- d. slice timing correction
                % adjusting for axial slice offsets, using interpolation, based on specified tpatt field
                unix(sprintf('3dTshift -prefix %s/func%u_tshift.nii -tpattern %s %s/func%u_ricor.nii',opath1,nr,InputStruct_ssa.TPATTERN,opath1,nr));
            end
        else
            disp('skipping slice time correction...')
        end

    %% alignment processing...

        % NB: we get mindisp reference for each run, even though we align all
        % to the run-1 reference atm. potentially useful for quantifying in-run vs between-run displacements

        if ~exist(sprintf('%s/func%u_2std_motref_masked.nii',opath1,nr),'file')
            % preparing to align... first reorient and deoblique
            unix(sprintf('3dWarp -oblique2card -prefix %s/func%u_deob.nii -wsinc5 %s/func%u_tshift.nii' , opath1,nr, opath1,nr));
            unix(sprintf('fslreorient2std %s/func%u_deob.nii %s/func%u_2std.nii',opath1,nr,opath1,nr));
            unix(sprintf('gunzip %s/func%u_2std.nii.gz',opath1,nr));
            if exist(sprintf('%s/func%u_2std.nii',opath1,nr),'file')
                unix(sprintf('rm %s/func%u_deob.nii', opath1,nr)); % not really useful intermediate - delete it
            else
                error('failed to create deob/reoriented functional file')
            end
            % final, best functional mask
            unix(sprintf('3dAutomask -prefix %s/func%u_2std_mask.nii %s/func%u_2std.nii',opath1,nr,opath1,nr)); % ** gets a decent functional mask
            % extract zero-motion "reference" volume and also mask it 
            unix(sprintf('3dTcat -prefix %s/func%u_2std_motref.nii ''%s/func%u_2std.nii[%u]''',opath1,nr,opath1,nr,imed_afni));
            unix(sprintf('3dcalc -prefix %s/func%u_2std_motref_masked.nii -a %s/func%u_2std_motref.nii -b %s/func%u_2std_mask.nii -expr ''a*b''',opath1,nr,opath1,nr,opath1,nr));
        else
            disp('skipping alignment prep...')
        end

        if strcmpi(InputStruct_ssa.PIPELINES.WARP{1},'AFNI')

            disp('AFNI-style alignment');

            if ~exist( sprintf('%s/func%u_motcor.nii',opath3f,nr),'file')
                %-- a. rigid within-run alignment of all func volumes to func run-1 refbrick ( 
                unix(sprintf('3dvolreg -zpad 1 -base %s/func1_2std_motref.nii -1Dfile %s/func%u_mpe -prefix %s/func%u_motcor.nii -cubic -1Dmatrix_save %s/func%u_motmat.1D %s/func%u_2std.nii', opath1,   opath3f,nr,opath3f,nr,opath3f,nr,opath1,nr));
            else
                disp('skipping motcor...')
            end

            if nr==1 && ~exist( sprintf('%s/anatSS.x_alj.nii',opath3f),'file')        
                %-- b. rigid alignment of run-1 func refbrick to run-1 t1 anat (stripped) --> ONLY FOR BASEBRICK OF FIRST RUN 
                skulstrip_anatomic = sprintf('%s/anatSS.x.nii',opath3a);
                unix(sprintf('align_epi_anat.py -anat2epi -anat %s -suffix _alj -epi %s/func1_2std_motref_masked.nii -epi_base 0  -epi_strip None  -anat_has_skull no  -ginormous_move -deoblique off -cost lpc+ZZ -volreg off -tshift off',skulstrip_anatomic,opath1));
                % convert output to .nii format, push to correct directory
                unix(sprintf('3dAFNItoNIFTI anatSS.x_alj+orig -prefix anatSS.x_alj.nii'));
                unix(sprintf('rm anatSS.x_alj+orig.BRIK anatSS.x_alj+orig.BRIK.gz anatSS.x_alj+orig.HEAD'));
                unix(sprintf('mv anatSS.x_alj* %s',opath3f));
            else
                disp('skipping t1 to epi...')
            end

            if ~exist( sprintf('%s/alg%u.x.affwarp.1D',opath3f,nr),'file')        
                %-- c. concatenate volreg/epi2anat/tlrc xforms
                unix(sprintf('cat_matvec -ONELINE %s/anatQQ.x.aff12.1D %s/anatSS.x_alj_mat.aff12.1D -I %s/func%u_motmat.1D > %s/alg%u.x.affwarp.1D', opath3a,opath3f,opath3f,nr,opath3f,nr));
            else
                disp('skipping warp concat...')
            end

            if ~exist( sprintf('%s/func%u_warped.nii',opath4f,nr),'file')        
                %-- d. apply concatenated xform: volreg/epi2anat/tlrc/NLtlrc; then apply non-linear standard-space warp
                unix(sprintf('3dNwarpApply -master %s/anatQQ.x.nii -dxyz 3 -source %s/func%u_2std.nii -nwarp "%s/anatQQ.x_WARP.nii %s/alg%u.x.affwarp.1D" -prefix ren%u.x.nlin',opath3a,opath1,nr,opath3a,opath3f,nr,nr));
                unix(sprintf('3dAFNItoNIFTI ren%u.x.nlin+tlrc -prefix ren%u.x.nlin.nii',nr,nr));
                unix(sprintf('rm ren%u.x.nlin+tlrc.BRIK ren%u.x.nlin+tlrc.BRIK.gz ren%u.x.nlin+tlrc.HEAD',nr,nr,nr));
                unix(sprintf('mv ren%u.x.nlin.nii %s/func%u_warped.nii',nr,opath4f,nr));
            else
                disp('skipping applywarp...')
            end

            % make mask, mean, sd maps --> for runs >1, this is just kept for qc purposes
            unix(sprintf('3dAutomask -prefix %s/func%u_warped_mask.nii %s/func%u_warped.nii',opath4f,nr,opath4f,nr));
            unix(sprintf('3dTstat -mean  -prefix %s/func%u_warped_tav.nii %s/func%u_warped.nii',opath4f,nr,opath4f,nr));
            unix(sprintf('3dTstat -stdev -prefix %s/func%u_warped_tsd.nii %s/func%u_warped.nii',opath4f,nr,opath4f,nr));
            % tidied up mask - for group masking
            if nr==1 && ~exist( sprintf('%s/func%u_warped_mask_clean.nii',opath4f,nr),'file')
                % tidied up functional mask in new space
                unix(sprintf('3dresample -master %s/func%u_warped_mask.nii -input %s/anatQQ.x.nii -prefix %s/anat_resam.nii',opath4f,nr,opath3a,opath3f));
                unix(sprintf('3dmask_tool -dilate_input 5 -5 -fill_holes -input %s/anat_resam.nii -prefix %s/anat_resam_mask.nii',opath3f,opath3f))
                unix(sprintf('3dmask_tool -input %s/func%u_warped_mask.nii %s/anat_resam_mask.nii -inter -prefix %s/func%u_warped_mask_clean.nii',opath4f,nr,opath3f,opath4f,nr))
            else
                disp('skipping newspace masking...')
            end

            tisslist = {'CSF','GM','WM'}; % tissues in increasing order of T1 intensity
            for i=1:3
                if nr==1 && ~exist( sprintf('%s/seg/anat_seg_%s_resam.nii',opath3f,tisslist{i}),'file')
                    unix(sprintf('3dresample -master %s/func1_warped_mask_clean.nii -input %s/seg/anat_seg_%s_warped.nii -prefix %s/seg/anat_seg_%s_resam.nii',opath4f,opath3a,tisslist{i},opath3f,tisslist{i}));
                else
                    disp('skipping tissue seg warping...')
                end
            end

        elseif strcmpi(PipeStruct.WARP{1},'ANTS')

            error('ants func to template warping is not yet available -- still in alpha testing ):');

        end
    end

end

%% INTERMEZZO A : constructing group-level tissue maps

% warp-specific path
catpath = [ParamStruct_aug.WarpPrefix{1},'.',ParamStruct_aug.VOXRES,'mm'];

if ~exist( [outpath,'/_group_level/masks/',catpath,'/func_brain_mask_grp.nii'], 'file') % Brain mask
    % consensus brain mask
    for ns=1:numel(subject_list)
        % check existence of subject specific struct file
        if exist(fullfile( outpath,subject_list{ns},'InputStruct_ssa.mat'),'file')
            load( fullfile( outpath,subject_list{ns},'InputStruct_ssa.mat') )
        else
            error('cannot find Input struct file for subject: %s \n');
        end
    
        opath4f= fullfile(outpath,InputStruct_ssa.PREFIX,'func_proc_p2',catpath);
        Vw=load_untouch_nii(sprintf('%s/func1_warped_mask_clean.nii',opath4f));
        if ns==1
            volref = double(Vw.img);
        else
            volref = volref + double(Vw.img);
        end
    end
    volref = volref./numel(subject_list);
    Vw.img = double(volref>0.50); % included in majority of individuals
    save_untouch_nii(Vw,[outpath,'/_group_level/masks/',catpath,'/func_brain_mask_grp.nii']);
else
    Vw = load_untouch_nii([outpath,'/_group_level/masks/',catpath,'/func_brain_mask_grp.nii']);
end
mask_brain = double(Vw.img);

if ~exist( [outpath,'/_group_level/masks/',catpath,'/func_CSF_mask_grp.nii'], 'file') % CSF mask
    % consensus CSF mask
    for ns=1:numel(subject_list)
        % check existence of subject specific struct file
        if exist(fullfile( outpath,subject_list{ns},'InputStruct_ssa.mat'),'file')
            load( fullfile( outpath,subject_list{ns},'InputStruct_ssa.mat') )
        else
            error('cannot find Input struct file for subject: %s \n');
        end
    
        opath3f= fullfile(outpath,InputStruct_ssa.PREFIX,'func_proc_p1',catpath);
        Vw=load_untouch_nii(sprintf('%s/seg/anat_seg_CSF_resam.nii',opath3f));
        if ns==1
            volref = double(Vw.img);
        else
            volref = volref + double(Vw.img);
        end
    end
    volref = volref./numel(subject_list);
    Vw.img = clust_up( double( volref > 0.90 ) ,20);
    save_untouch_nii(Vw,[outpath,'/_group_level/masks/',catpath,'/func_CSF_mask_grp.nii']);
else
    Vw = load_untouch_nii([outpath,'/_group_level/masks/',catpath,'/func_CSF_mask_grp.nii']);
end
mask_csf = double(Vw.img);

if ~exist( [outpath,'/_group_level/masks/',catpath,'/func_WM_mask_grp.nii'], 'file') % WM mask
    % consensus WM mask
    for ns=1:numel(subject_list)
        % check existence of subject specific struct file
        if exist(fullfile( outpath,subject_list{ns},'InputStruct_ssa.mat'),'file')
            load( fullfile( outpath,subject_list{ns},'InputStruct_ssa.mat') )
        else
            error('cannot find Input struct file for subject: %s \n');
        end
    
        opath3f= fullfile(outpath,InputStruct_ssa.PREFIX,'func_proc_p1',catpath);
        Vw=load_untouch_nii(sprintf('%s/seg/anat_seg_WM_resam.nii',opath3f));
        if ns==1
            volref = double(Vw.img);
        else
            volref = volref + double(Vw.img);
        end
    end
    volref = volref./numel(subject_list);
    volref = smooth3( volref, 'gaussian',[7 7 7], 0.85);
    Vw.img = clust_up( double( volref > 0.90 ) ,20);
    save_untouch_nii(Vw,[outpath,'/_group_level/masks/',catpath,'/func_WM_mask_grp.nii']);
else
    Vw = load_untouch_nii([outpath,'/_group_level/masks/',catpath,'/func_WM_mask_grp.nii']);
end
mask_wm = double(Vw.img);

if ~exist( [outpath,'/_group_level/masks/',catpath,'/func_tSD_mask_grp.nii'], 'file') % t-std mask
    % consensus VAR mask
    for ns=1:numel(subject_list)
        ns,
        % check existence of subject specific struct file
        if exist(fullfile( outpath,subject_list{ns},'InputStruct_ssa.mat'),'file')
            load( fullfile( outpath,subject_list{ns},'InputStruct_ssa.mat') )
        else
            error('cannot find Input struct file for subject: %s \n');
        end
    
        opath4f= fullfile(outpath,InputStruct_ssa.PREFIX,'func_proc_p2',catpath);
        Vw=load_untouch_nii(sprintf('%s/func1_warped_tsd.nii',opath4f));
        if ns==1
            volref = double(Vw.img);
        else
            volref = volref + double(Vw.img);
        end
    end
    volref = volref./numel(subject_list);
    Vw.img = clust_up( double( volref > prctile(volref(mask_brain>0),90) ) .* mask_brain ,20);
    save_untouch_nii(Vw,[outpath,'/_group_level/masks/',catpath,'/func_tSD_mask_grp.nii']);
else
    Vw = load_untouch_nii([outpath,'/_group_level/masks/',catpath,'/func_tSD_mask_grp.nii']);
end
mask_tsd = double(Vw.img);

if ~exist( [outpath,'/_group_level/masks/',catpath,'/func_GM_mask_grp.nii'], 'file') % GM inclusive mask
    % liberal GM mask
    tisslist = {'CSF','GM','WM'};
    for i=1:3
        for ns=1:numel(subject_list)
            % check existence of subject specific struct file
            if exist(fullfile( outpath,subject_list{ns},'InputStruct_ssa.mat'),'file')
                load( fullfile( outpath,subject_list{ns},'InputStruct_ssa.mat') )
            else
                error('cannot find Input struct file for subject: %s \n');
            end
        
            opath3f= fullfile(outpath,InputStruct_ssa.PREFIX,'func_proc_p1',catpath);
            Vw=load_untouch_nii(sprintf('%s/seg/anat_seg_%s_resam.nii',opath3f,tisslist{i}));
            if ns==1
                volref = double(Vw.img);
            else
                volref = volref + double(Vw.img);
            end
        end
        volref = volref./numel(subject_list);
        volref = smooth3( volref, 'gaussian',[7 7 7], 0.85);
        volvec(:,i) = volref(mask_brain>0);
    end
    tmp = mask_brain; tmp(tmp>0)= double( (volvec(:,2) ./ sum(volvec,2)) >=0.33 );
    Vw.img = tmp;
    save_untouch_nii(Vw,[outpath,'/_group_level/masks/',catpath,'/func_GM_mask_grp.nii']);
else
    Vw = load_untouch_nii([outpath,'/_group_level/masks/',catpath,'/func_GM_mask_grp.nii']);
end
mask_gm = double(Vw.img);

if ~exist( [outpath,'/_group_level/brain_maps/',catpath,'/func_tAV_grp.nii'], 'file') % group mean tav epi map
    for ns=1:numel(subject_list)
        ns,
        % check existence of subject specific struct file
        if exist(fullfile( outpath,subject_list{ns},'InputStruct_ssa.mat'),'file')
            load( fullfile( outpath,subject_list{ns},'InputStruct_ssa.mat') )
        else
            error('cannot find Input struct file for subject: %s \n');
        end
    
        opath4f= fullfile(outpath,InputStruct_ssa.PREFIX,'func_proc_p2',catpath);
        Vw=load_untouch_nii(sprintf('%s/func1_warped_tav.nii',opath4f));
        if ns==1
            volref = double(Vw.img);
        else
            volref = volref + double(Vw.img);
        end
    end
    Vw.img = volref./numel(subject_list);
    save_untouch_nii(Vw,[outpath,'/_group_level/brain_maps/',catpath,'/func_tAV_grp.nii']); % group mean sd epi map
end

if ~exist( [outpath,'/_group_level/brain_maps/',catpath,'/func_tSD_grp.nii'], 'file')
    for ns=1:numel(subject_list)
        ns,
        % check existence of subject specific struct file
        if exist(fullfile( outpath,subject_list{ns},'InputStruct_ssa.mat'),'file')
            load( fullfile( outpath,subject_list{ns},'InputStruct_ssa.mat') )
        else
            error('cannot find Input struct file for subject: %s \n');
        end
    
        opath4f= fullfile(outpath,InputStruct_ssa.PREFIX,'func_proc_p2',catpath);
        Vw=load_untouch_nii(sprintf('%s/func1_warped_tsd.nii',opath4f));
        if ns==1
            volref = double(Vw.img);
        else
            volref = volref + double(Vw.img);
        end
    end
    Vw.img = volref./numel(subject_list);
    save_untouch_nii(Vw,[outpath,'/_group_level/brain_maps/',catpath,'/func_tSD_grp.nii']);
end

%% INTERMEZZO B : constructing group-level spatial weighting on tissue maps

masklist = {'CSF','WM','tSD'};
parcdir  = [outpath,'/_group_level/parcellations/',catpath];
if ~exist([parcdir,'/Ugrp_wm.nii'], 'file') || ~exist([parcdir,'/Ugrp_csf.nii'], 'file') || ~exist([parcdir,'/Ugrp_var.nii'], 'file')
    %
    for i=1:numel(masklist)
        mkdir_r([parcdir,'/tmp']);
        maskname = [outpath,'/_group_level/masks/',catpath,'/func_',masklist{i},'_mask_grp.nii'];
        Mtmp = load_untouch_nii(maskname);
        ucat=[];
        unix(sprintf('rm %s/tmp/volblur.nii',parcdir))
        for ns=1:numel(subject_list)
            [ns numel(subject_list)],
            % check existence of subject specific struct file
            if exist(fullfile( outpath,subject_list{ns},'InputStruct_ssa.mat'),'file')
                load( fullfile( outpath,subject_list{ns},'InputStruct_ssa.mat') )
            else
                error('cannot find Input struct file for subject: %s \n');
            end
            % NB: only uses run-1 per subject
            opath4f= fullfile(outpath,InputStruct_ssa.PREFIX,'func_proc_p2',catpath);
            unix(sprintf('3dBlurInMask -input %s/func1_warped.nii -prefix %s/tmp/volblur.nii -mask %s -fwhm 6',opath4f,parcdir,maskname));
            Vtmp = load_untouch_nii(sprintf('%s/tmp/volblur.nii',parcdir));
            volmat = nifti_to_mat(Vtmp,Mtmp); npc_10 = floor(0.1*size(volmat,2));
            volmat = bsxfun(@rdivide, bsxfun(@minus,volmat,mean(volmat,2)), std(volmat,0,2)+eps);
            [u,~,~]=svd( volmat,'econ');
            ucat = [ucat, u(:,1:npc_10)];
            npc_10_list(ns,1) = npc_10;
            unix(sprintf('rm %s/tmp/volblur.nii',parcdir))
        end
        [Ugrp,~,~] = svd( ucat,'econ');
        Ugrp = Ugrp(:,1:round(median(npc_10_list)));
        save([parcdir,'/Ugrp_',masklist{i},'.mat'],'Ugrp','npc_10_list');
    
        clear TMPVOL;
        for(p=1:size(Ugrp,2) )
            tmp=double(Mtmp.img);
            tmp(tmp>0)= Ugrp(:,p);
            TMPVOL(:,:,:,p) = tmp;
        end
        nii=Vtmp;
        nii.img = TMPVOL;
        nii.hdr.dime.datatype = 16;
        nii.hdr.hist = Vtmp.hdr.hist;
        nii.hdr.dime.dim(5) = size(Ugrp,2);
        save_untouch_nii(nii,[parcdir,'/Ugrp_',masklist{i},'.nii']); 
    
        unix(sprintf('rm %s/*.nii', [parcdir,'/tmp']))
        unix(sprintf('rmdir %s', [parcdir,'/tmp']))
    end
end

%% FUNC PROCESSING...BLOCK2...
for ns=1:numel(subject_list)

    % check existence of subject specific struct file
    if exist(fullfile( outpath,subject_list{ns},'InputStruct_ssa.mat'),'file')
        load( fullfile( outpath,subject_list{ns},'InputStruct_ssa.mat') )
    else
        error('cannot find Input struct file for subject: %s \n');
    end

    fprintf('\n===> func-proc. round-2 now on subj %u/%u: %s...\n',ns,numel(subject_list),subject_list{ns}),

    %%%%%========== compatibility adjustments for BLOCK2 ... TASKREG only
    missing_physio=0;
    for nr=1:InputStruct_ssa.N_func
        if isempty(InputStruct_ssa.frun(nr).PHYSIO_filename)
            missing_physio=missing_physio+1;
        end
    end
    if analysis_struct.uses_taskfile==0  && ~(strcmpi(InputStruct_ssa.PIPELINES.TASKREG{1},'0') || strcmpi(InputStruct_ssa.PIPELINES.TASKREG{1},'OFF'))
        warning('your analysis model does not specify a task design. Turning TASKREG off!')
        InputStruct_ssa.PIPELINES.TASKREG{1}='0';
    end
    %%%%%========== compatibility adjustments, done

    % appending info about run-length (for regressor consrtuction and interp-contrast-list)
    for nr = 1:InputStruct_ssa.N_func
        %
        opath0   = fullfile(outpath,InputStruct_ssa.PREFIX,'rawdata');
        hdr      = load_nii_hdr(sprintf('%s/func%u.nii',opath0,nr));
        InputStruct_ssa.frun(nr).Nt_raw   = hdr.dime.dim(5);
        InputStruct_ssa.frun(nr).Nt_adj   = hdr.dime.dim(5) - InputStruct_ssa.frun(nr).DROP_first - InputStruct_ssa.frun(nr).DROP_last;
    end

    % extract contrasts specified for analysis --> appended into InputStruct, delete unformatted split-info field 
    InputStruct_ssa = interpret_contrast_list( InputStruct_ssa, analysis_struct, ParamStruct_aug.CONTRAST); % generate contrast list for each subject and run    
    %InputStruct_ssa = rmfield(InputStruct_ssa,'task_unformat');
    % extract seed contrasts SAA
    InputStruct_ssa = interpret_seed_list( InputStruct_ssa, analysis_struct, ParamStruct_aug.CONTRAST); % generate contrast list for each subject and run 
    %InputStruct_ssa = rmfield(InputStruct_ssa,'seed_unformat');

    % populating functional alignment dir
    mkdir_r( fullfile(outpath,InputStruct_ssa.PREFIX,'func_proc_p2',[ParamStruct_aug.WarpPrefix{1},'.',ParamStruct_aug.VOXRES,'mm'] ) );

    nr=1;
    %
    opathe  = fullfile(outpath,InputStruct_ssa.PREFIX,'func_proc_p2');
    opath3a = fullfile(outpath,InputStruct_ssa.PREFIX,'anat_proc',ParamStruct_aug.WarpPrefix{1});
    opath3f = fullfile(outpath,InputStruct_ssa.PREFIX,'func_proc_p1',[ParamStruct_aug.WarpPrefix{1},'.',ParamStruct_aug.VOXRES,'mm']);
    opath4f = fullfile(outpath,InputStruct_ssa.PREFIX,'func_proc_p2',[ParamStruct_aug.WarpPrefix{1},'.',ParamStruct_aug.VOXRES,'mm']);

    % spatial smoothing first
    unix(sprintf('3dmerge -prefix %s/func%u_warped_smo%s.nii -doall -1blur_fwhm %s %s/func%u_warped.nii',opath4f,nr,InputStruct_ssa.PIPELINES.SMOOTH{1},InputStruct_ssa.PIPELINES.SMOOTH{1},opath4f,nr));

    % loading data into mats:
    VSstr = sprintf('%s/func%u_warped_smo%s.nii',opath4f,nr,InputStruct_ssa.PIPELINES.SMOOTH{1});
    MBstr = ([outpath,'/_group_level/masks/',[ParamStruct_aug.WarpPrefix{1},'.',ParamStruct_aug.VOXRES,'mm'],'/func_brain_mask_grp.nii']);
    VS = load_untouch_nii(VSstr);
    MB = load_untouch_nii(MBstr);
    % doanload, store mean, and subtract out
    maskS = double(MB.img);
    volmatS = nifti_to_mat(VS,MB); 

    % regressione processingue ...
    
    % === construct legendre poly
    if strcmpi(InputStruct_ssa.PIPELINES.DETREND{1},'A')
        do = 1+floor( (InputStruct_ssa.TR_MSEC./1000) * InputStruct_ssa.frun(nr).Nt_adj ./ 150 );
    else
        do = str2num( InputStruct_ssa.PIPELINES.DETREND{1} );
    end
    xdet = det_regressor_builder( do, InputStruct_ssa.frun(nr).Nt_adj ); % we explicitly calculate these regressors, for examining design matrix later

    % === construct mpe stuff
    xmpe = mpe_regressor_builder( fullfile(outpath,InputStruct_ssa.PREFIX,'func_proc_p1',[ParamStruct_aug.WarpPrefix{1},'.',ParamStruct_aug.VOXRES,'mm'],sprintf('func%u_mpe',nr)), InputStruct_ssa.PIPELINES.MOTREG{1} );

    % === construct seed regs
    fstr = sprintf('%s/func%u_warped.nii',opath4f,nr);
    wcv_masks{1}   = [outpath,'/_group_level/masks/',[ParamStruct_aug.WarpPrefix{1},'.',ParamStruct_aug.VOXRES,'mm'],'/func_WM_mask_grp.nii'];
    wcv_masks{2}   = [outpath,'/_group_level/masks/',[ParamStruct_aug.WarpPrefix{1},'.',ParamStruct_aug.VOXRES,'mm'],'/func_CSF_mask_grp.nii'];
    wcv_masks{3}   = [outpath,'/_group_level/masks/',[ParamStruct_aug.WarpPrefix{1},'.',ParamStruct_aug.VOXRES,'mm'],'/func_tSD_mask_grp.nii'];
    wvc_weights{1} = [outpath,'/_group_level/parcellations/',[ParamStruct_aug.WarpPrefix{1},'.',ParamStruct_aug.VOXRES,'mm'],'/Ugrp_WM.nii'];
    wvc_weights{2} = [outpath,'/_group_level/parcellations/',[ParamStruct_aug.WarpPrefix{1},'.',ParamStruct_aug.VOXRES,'mm'],'/Ugrp_CSF.nii'];
    wvc_weights{3} = [outpath,'/_group_level/parcellations/',[ParamStruct_aug.WarpPrefix{1},'.',ParamStruct_aug.VOXRES,'mm'],'/Ugrp_tSD.nii'];
    %
    xseed = seed_regressor_builder( fstr, wcv_masks, wvc_weights, InputStruct_ssa.PIPELINES.ROIREG{1} );
    
    % === construct signal regs
    if strcmpi(InputStruct_ssa.PIPELINES.LOPASS{1},'0') || strcmpi(InputStruct_ssa.PIPELINES.LOPASS{1},'OFF')
        Xsignal = [];
    elseif strcmpi(InputStruct_ssa.PIPELINES.LOPASS{1},'1') || strcmpi(InputStruct_ssa.PIPELINES.LOPASS{1},'ON')
        Xsignal = InputStruct_ssa.task_contrast{1}.design_mat;
        Xsignal = bsxfun(@rdivide, Xsignal, sqrt(sum(Xsignal.^2)));
    else
        error('unknown taskreg option')
    end

    % === filter applied to regressors
    Xnoise = [xdet(:,2:end),xmpe,xseed]; % full noise matrix, normed to unit length
    Xnoise = bsxfun(@rdivide,Xnoise,sqrt(sum(Xnoise.^2)));
    [ output ] = GLM_model_fmri( volmatS, [0], [Xnoise], [Xsignal], 1, 1 );
    volmatF = output.vol_denoi;

    % === low-pass filtering, if requested
    if strcmpi(InputStruct_ssa.PIPELINES.LOPASS{1},'0') || strcmpi(InputStruct_ssa.PIPELINES.LOPASS{1},'OFF')
        disp('skipping lopass');
    elseif strcmpi(InputStruct_ssa.PIPELINES.LOPASS{1},'1') || strcmpi(InputStruct_ssa.PIPELINES.LOPASS{1},'ON')
        [volmatF] = quick_lopass( volmatF, (InputStruct_ssa.TR_MSEC./1000) );
    else
        error('unknown lopass option');
    end

    %% Transition to the analysis phase ... using "volmatF"

    oudirre = fullfile(outpath,InputStruct_ssa.PREFIX,'func_proc_p2',[ParamStruct_aug.WarpPrefix{1},'.',ParamStruct_aug.VOXRES,'mm'] );

    % first, export fully processed data to .nii and .mat formats
    if ~exist([oudirre,'/func',num2str(nr),'_fullproc.nii'],'file')
        clear TMPVOL;
        for(p=1:size(volmatF,2) )
            tmp=double(MB.img);
            tmp(tmp>0)= volmatF(:,p);
            TMPVOL(:,:,:,p) = tmp;
        end
        nii=VS;
        nii.img = TMPVOL;
        nii.hdr.dime.datatype = 16;
        nii.hdr.hist = VS.hdr.hist;
        nii.hdr.dime.dim(5) = size(volmatF,2);
        save_untouch_nii(nii,[oudirre,'/func',num2str(nr),'_fullproc.nii']); 
        save([oudirre,'/func',num2str(nr),'_fullproc.mat'],'volmatF');
    end

    % if seed-based analysis, make sure to produce appropriately resampled copies.....
    if analysis_struct.uses_roifile>0
        seedpath = fullfile(outpath,InputStruct_ssa.PREFIX,'func_seeds');
        opath4f= fullfile(outpath,InputStruct_ssa.PREFIX,'func_proc_p2',[ParamStruct_aug.WarpPrefix{1},'.',ParamStruct_aug.VOXRES,'mm']);
        mkdir_r([seedpath,'/',[ParamStruct_aug.WarpPrefix{1},'.',ParamStruct_aug.VOXRES,'mm']]);
        for i=1:numel(InputStruct_ssa.seed_unformat.seed)
            if ~exist(sprintf('%s/%s/%s.resam.nii',seedpath,[ParamStruct_aug.WarpPrefix{1},'.',ParamStruct_aug.VOXRES,'mm'], InputStruct_ssa.seed_unformat.seed(i).name),'file')
                unix(sprintf('cp %s %s/%s.nii',InputStruct_ssa.seed_unformat.seed(i).location, seedpath, InputStruct_ssa.seed_unformat.seed(i).name ));
                unix(sprintf('3dresample -master %s/func1_warped_mask_clean.nii -input  %s/%s.nii -prefix  %s/%s/%s.resam.nii',opath4f,seedpath, InputStruct_ssa.seed_unformat.seed(i).name,seedpath,[ParamStruct_aug.WarpPrefix{1},'.',ParamStruct_aug.VOXRES,'mm'], InputStruct_ssa.seed_unformat.seed(i).name));
            else
                disp('seedmap found! skipping to next');
            end
        end     
    end
    if analysis_struct.uses_taskfile>0
        %--> no special args atm
    end

    if ~strcmpi( ParamStruct_aug.ANALYSIS, 'NONE') && ~strcmpi( ParamStruct_aug.CONTRAST, 'NONE')

%         % then, execute the chosen contrast+analysis
%         disp('---analysis running---');
%         
%         params=[];
%         params.
%         if strcmpi(ParamStruct_aug.ANALYSIS,'CONN')
%             for i=1:numel(InputStruct_ssa.seed_info.contrast)
%                 V  = load_untouch_nii( sprintf('%s/%s/%s.resam.nii',seedpath,[ParamStruct_aug.WarpPrefix{1},'.',ParamStruct_aug.VOXRES,'mm'], InputStruct_ssa.seed_unformat.seed(i).name) );
%                 params.seedmat(:,i) = nifti_to_mat(V,MB);
%             end
%         end
%     
%         % get function handle for analysis model of interest
%         currPath=pwd;                   % get current path
%         cd(analysis_struct.filepath);    % jump to module directory
%         pfun= str2func(analysis_struct.model_name); % get function handle
%         cd(currPath);                   % jump back to current path
%         % call function
%         output = pfun( volmatF, InputStruct_ssa.task_info, params );  

    end

end


