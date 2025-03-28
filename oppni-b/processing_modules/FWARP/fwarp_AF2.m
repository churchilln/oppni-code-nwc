function fwarp_AF2( Funcfile_set, prefix_set, odir1, odir2, base_set, Anatloc, ParamCell )
%
% .fwarp_AF2:
% .functional warping using AFNI utilities
% .adapted from afni_proc.py script, using outputs from @SSwarper

% an expanded "optimizer" built in for initial steps
% > expects certain formatting in Anatloc folder!

if isempty(ParamCell) || isempty(ParamCell{1})
    anatType = 'T1';
else
    anatType = ParamCell{1};
end
if numel(ParamCell)<2 || isempty(ParamCell{2})
    searchSet = 'trans+move';
else
    searchSet = ParamCell{2};
end

%% setting param list

if strcmpi(searchSet,'trans+move') || strcmpi(searchSet,'move+trans')

    param_list = {'',               ' -rigid_body',...
                  ' -big_move',      ' -big_move -rigid_body',...
                  ' -giant_move',    ' -giant_move -rigid_body',...
                  ' -ginormous_move',' -ginormous_move -rigid_body'};

elseif strcmpi(searchSet,'trans')

    param_list = {'', ' -rigid_body'};

elseif strcmpi(searchSet,'move')

    param_list = {'', ' -big_move', ' -giant_move', ' -ginormous_move'};

end

%% to alignment...

% prefix for temp files
pref = [odir1,'/__opptmp_p2func_warp'];
% num func runs
N_func  = numel(Funcfile_set);
% basefile
Basefile_pref = sprintf('%s/motref',odir1); % *** NOTE WE ARE FIXED ON RUN-1 MOTREF FOR ALIGNMENT

for nr=1:N_func

    if ~exist( sprintf('%s/%s_warped.nii.gz',odir2,prefix_set{nr}) ,'file')
    
        % build directory struct recursively
        unix(sprintf('mkdir -p %s',pref));
    
        % extract zero-motion "reference" volume and also mask it (use this for motref and/or QC) 
        unix(sprintf('3dAutomask -prefix %s/%s_mask.nii.gz %s',odir1,prefix_set{nr},Funcfile_set{nr})); % ** gets a decent functional mask

        if nr==1 % Basefile constructed using only run-1 basebrick!
            unix(sprintf('3dTcat -prefix %s.nii.gz ''%s[%u]''',Basefile_pref,Funcfile_set{nr},base_set(nr)));
            unix(sprintf('3dcalc -prefix %s_masked.nii.gz -a %s.nii.gz -b %s/%s_mask.nii.gz -expr ''a*b''',Basefile_pref,Basefile_pref,odir1,prefix_set{nr}));
        end

        %-- a. rigid within-run alignment of all func volumes to func run-1 refbrick ( 
        unix(sprintf('3dvolreg -zpad 1 -base %s.nii.gz -1Dfile %s/%s_mpe -prefix %s/%s_motcor.nii.gz -cubic -1Dmatrix_save %s/%s_motmat.aff12.1D %s', ...
            Basefile_pref,odir1,prefix_set{nr},odir1,prefix_set{nr},odir1,prefix_set{nr},Funcfile_set{nr}));
    
        %-- b. rigid alignment of run-1 func refbrick to run-1 t1 anat (stripped) --> ONLY FOR BASEBRICK OF FIRST RUN 
        if nr==1

            currPath=pwd;   % get current path
            cd(pref);       % jump to temp directory --> because this script creates lots of local derivatives 

            if     strcmpi(anatType,'T1')
                %unix(sprintf('align_epi_anat.py -anat2epi -anat %s/anat_procss.nii.gz -suffix _alj -epi %s_masked.nii.gz -epi_base 0  -epi_strip None  -anat_has_skull no  -ginormous_move -deoblique off -cost lpc+ZZ -volreg off -tshift off', Anatloc,Basefile_pref));
                costtype = 'lpc+ZZ';
            elseif strcmpi(anatType,'FLAIR') || strmpi(anatType,'T2')
                %unix(sprintf('align_epi_anat.py -anat2epi -anat %s/anat_procss.nii.gz -suffix _alj -epi %s_masked.nii.gz -epi_base 0  -epi_strip None  -anat_has_skull no  -ginormous_move -deoblique off -cost lpa+ZZ -volreg off -tshift off', Anatloc,Basefile_pref));
                costtype = 'lpa+ZZ';
            else
                error('unrecognized anat-type');
            end

            % load motref
            A=load_untouch_niiz(sprintf('%s_masked.nii.gz',Basefile_pref));
            a_img=double(A.img>0);

            for im=1:numel(param_list)
                unix(sprintf('mkdir %s/meth%u',pref,im)); % create subdirectory
                unix(sprintf('align_epi_anat.py -anat2epi -anat %s/anat_procss.nii.gz -suffix _alj -epi %s_masked.nii.gz -epi_base 0  -epi_strip None  -anat_has_skull no -deoblique off -cost %s%s -volreg off -tshift off',...
                    Anatloc,Basefile_pref,costtype,param_list{im}));             
                unix(sprintf('rm anat_procss_alj+orig.BRIK anat_procss_alj+orig.BRIK.gz anat_procss_alj+orig.HEAD anat_procss_alj_e2a_only_mat.aff12.1D'));
                unix(sprintf('mv anat_procss_alj_mat.aff12.1D %s/meth%u',pref,im));
                unix(sprintf('3dAllineate -cubic -prefix %s/meth%u/anat_procss_loResAlign.nii.gz -1Dmatrix_apply %s/meth%u/anat_procss_alj_mat.aff12.1D -base %s_masked.nii.gz -input %s/anat_procss.nii.gz',...
                    pref,im, pref,im, Basefile_pref, Anatloc ));
                B=load_untouch_niiz(sprintf('%s/meth%u/anat_procss_loResAlign.nii.gz', pref,im ));
                b_img=double(B.img>0);
                % jaccard index viz motref
                jacc(im,1) = sum( a_img(:).*b_img(:) )./ sum( (a_img(:)+b_img(:))>0 );
            end
            [optval,optix] = max(jacc);
            writematrix(jacc, sprintf('%s/jacc_vals.txt',odir1));
            fprintf('\noptimal EPI-T1 approach is method %u - overlap: %.03f\n\n',optix,optval);
            % move to main directory!
            unix(sprintf('mv %s/meth%u/anat_procss_alj_mat.aff12.1D %s',pref,optix,odir1));
            unix(sprintf('mv %s/meth* %s',pref,odir1)); % for now, keep different method outputs

            cd(currPath);  % jump back to current path
        end
    
        %-- c. concatenate volreg/epi2anat/tlrc xforms
        unix(sprintf('cat_matvec -ONELINE %s/anatQQ.aff12.1D %s/anat_procss_alj_mat.aff12.1D -I %s/%s_motmat.aff12.1D > %s/%s_netaff.aff12.1D', ...
            Anatloc,odir1,odir1,prefix_set{nr},odir1,prefix_set{nr}));
    
        %-- d. apply concatenated xform: volreg/epi2anat/tlrc/NLtlrc; then apply non-linear standard-space warp
        unix(sprintf('3dNwarpApply -master %s/anat_warped.nii.gz -dxyz 3 -source %s -nwarp "%s/anatQQ_WARP.nii.gz %s/%s_netaff.aff12.1D" -prefix %s/%s_warped.nii.gz',...
            Anatloc,Funcfile_set{nr},Anatloc,odir1,prefix_set{nr},odir2,prefix_set{nr}));
    
        unix(sprintf('rm -rf %s',pref));

    else
        disp('afni-funcwarp already exists!')
    end

% %     nr=1; % just run-1 for steps beloq - QUALIRTY chex
% % 
% %     if ~exist( sprintf('%s/%s_warped_nomc.nii.gz',odir2,prefix_set{nr}) ,'file')
% %         nr=1;
% %         %-- c. concatenate volreg/epi2anat/tlrc xforms
% %         unix(sprintf('cat_matvec -ONELINE %s/anatQQ.aff12.1D %s/anat_procss_alj_mat.aff12.1D -I > %s/%s_netaff_nomc.aff12.1D', ...
% %             Anatloc,odir1,odir1,prefix_set{nr}));
% %         %-- d. apply concatenated xform: volreg/epi2anat/tlrc/NLtlrc; then apply non-linear standard-space warp
% %         unix(sprintf('3dNwarpApply -master %s/anat_warped.nii.gz -dxyz 3 -source %s -nwarp "%s/anatQQ_WARP.nii.gz %s/%s_netaff_nomc.aff12.1D" -prefix %s/%s_warped_nomc.nii.gz',...
% %             Anatloc,Funcfile_set{nr},Anatloc,odir1,prefix_set{nr},odir2,prefix_set{nr}));
% %     end
% %     if ~exist( sprintf('%s/%s_warped_noal.nii.gz',odir2,prefix_set{nr}) ,'file')
% % 
% %         %-- c. concatenate volreg/epi2anat/tlrc xforms
% %         unix(sprintf('cat_matvec -ONELINE %s/anatQQ.aff12.1D %s/%s_motmat.aff12.1D > %s/%s_netaff_noal.aff12.1D', ...
% %             Anatloc,odir1,prefix_set{nr},odir1,prefix_set{nr}));
% %         %-- d. apply concatenated xform: volreg/epi2anat/tlrc/NLtlrc; then apply non-linear standard-space warp
% %         unix(sprintf('3dNwarpApply -master %s/anat_warped.nii.gz -dxyz 3 -source %s -nwarp "%s/anatQQ_WARP.nii.gz %s/%s_netaff_noal.aff12.1D" -prefix %s/%s_warped_noal.nii.gz',...
% %             Anatloc,Funcfile_set{nr},Anatloc,odir1,prefix_set{nr},odir2,prefix_set{nr}));
% %     end
% %     if ~exist( sprintf('%s/%s_warped_prmc.nii.gz',odir2,prefix_set{nr}) ,'file')
% % 
% %         %-- d. apply concatenated xform: volreg/epi2anat/tlrc/NLtlrc; then apply non-linear standard-space warp
% %         unix(sprintf('3dNwarpApply -master %s/anat_warped.nii.gz -dxyz 3 -source %s/%s_motcor.nii.gz -nwarp "%s/anatQQ_WARP.nii.gz %s/%s_netaff_nomc.aff12.1D" -prefix %s/%s_warped_prmc.nii.gz',...
% %             Anatloc,odir1,prefix_set{nr},Anatloc,odir1,prefix_set{nr},odir2,prefix_set{nr}));
% %     end

end
