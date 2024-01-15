function pwwarp_AF2( Funcfile_set, prefix_set, odir1, odir2, Basefile_masked, Anatloc )

% func warping same as AF1, but v minor tweaks for non-t1
% compatibility... see ln.25

% > expects certain formatting in Anatloc folder!

% odir => sprintf('%s/warp',opath2f)
% prefix_set => sprintf('func%u',nr) x N_func
% Funcfile_set => sprintf('%s/prewarp/func%u_2std.nii',opath2f,nr)  x N_func
% base_set => motref_0rel x N_func
% Anatloc => sprintf('%s',opath2a);
% AnatfileSS => sprintf('%s/anat_procss.nii',opath2a)
% Anat_Afffile => sprintf('%s/anatQQ.aff12.1D',opath2a)
% Anat_Warpfile => sprintf('%s/anatQQ_WARP.nii',opath2a)
% AnatfileWW => sprintf('%s/anat_warped.nii',opath2a)

% NB: when specifying alignment matrix outputs, make sure suffix is aff12.1D
%     throughout, otherwise cat_matvec will only do 1 line!

N_func = numel(Funcfile_set);

if ~exist(sprintf('%s/netaff.aff12.1D',odir1),'file')
    
    %-- b. rigid alignment of run-1 func refbrick to run-1 t1 anat (stripped) --> ONLY FOR BASEBRICK OF FIRST RUN 
    unix(sprintf('align_epi_anat.py -anat2epi -anat %s/anat_procss.nii.gz -suffix _alj -epi %s -epi_base 0  -epi_strip None  -anat_has_skull no  -ginormous_move -deoblique off -cost lpa+ZZ -volreg off -tshift off',...
        Anatloc,Basefile_masked));
    % convert output to .nii format, push to correct directory
    %unix(sprintf('3dAFNItoNIFTI anat_procss_alj+orig -prefix anat_procss_alj.nii.gz')); 
    % % --> just delete alinged without making copy for now, also delete the e2a_only file (unused!) 
    unix(sprintf('rm anat_procss_alj+orig.BRIK anat_procss_alj+orig.BRIK.gz anat_procss_alj+orig.HEAD anat_procss_alj_e2a_only_mat.aff12.1D'));
    unix(sprintf('mv anat_procss_alj_mat.aff12.1D %s',odir1));
    
    %-- c. concatenate volreg/epi2anat/tlrc xforms
    unix(sprintf('cat_matvec -ONELINE %s/anatQQ.aff12.1D %s/anat_procss_alj_mat.aff12.1D -I > %s/netaff.aff12.1D', ...
        Anatloc,odir1,odir1));

end

for nr=1:N_func

    if ~exist( sprintf('%s/%s_warped.nii.gz',odir2,prefix_set{nr}) ,'file')
        %-- d. apply concatenated xform: volreg/epi2anat/tlrc/NLtlrc; then apply non-linear standard-space warp
        unix(sprintf('3dNwarpApply -master %s/anat_warped.nii.gz -dxyz 3 -source %s -nwarp "%s/anatQQ_WARP.nii.gz %s/netaff.aff12.1D" -prefix %s/%s_warped.nii.gz',...
            Anatloc,Funcfile_set{nr},Anatloc,odir1,odir2,prefix_set{nr}));
    else
        disp('afni-funcwarp already exists!')
    end

end
