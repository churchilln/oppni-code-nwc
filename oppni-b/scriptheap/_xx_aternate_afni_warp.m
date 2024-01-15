            disp('AFNI-style alignment');

            if ~exist( sprintf('%s/func%u_motcor.nii',opath3f,nr),'file')
                %-- a. rigid within-run alignment of all func volumes to func refbrick
                unix(sprintf('3dvolreg -zpad 1 -base %s/func%u_2std_motref.nii -1Dfile %s/func%u_mpe -prefix %s/func%u_motcor.nii -cubic -1Dmatrix_save %s/func%u_motmat.1D %s/func%u_2std.nii', opath1,nr,opath3f,nr,opath3f,nr,opath3f,nr,opath1,nr));
            else
                disp('skipping motcor...')
            end
            
            %% ---> TO-DO: motref-n to motref-1 alignment step
            unix(sprintf('3dvolreg -zpad 1 -base %s/func1_2std_motref.nii -1Dfile %s/func%u_aln2run1_mpe -prefix %s/func%u_aln2run1_motref.nii -cubic -1Dmatrix_save %s/func%u_aln2run1.1D %s/func%u_2std_motref.nii', ...
                opath1,opath3f,nr,opath3f,nr,opath3f,nr,opath1,nr));
            
            if nr==1 && ~exist( sprintf('%s/anatSS.x_alj.nii',opath3f),'file')        
                %-- b. rigid alignment of run-1 func refbrick to run-1 t1 anat (stripped) --> ONLY FOR FIRST RUN 
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
                unix(sprintf('cat_matvec -ONELINE %s/anatQQ.x.aff12.1D %s/anatSS.x_alj_mat.aff12.1D -I %s/func%u_aln2run1.1D %s/func%u_motmat.1D > %s/alg%u.x.affwarp.1D', opath3a,opath3f,opath3f,nr,opath3f,nr,opath3f,nr));
            else
                disp('skipping warp concat...')
            end

            if ~exist( sprintf('%s/ren%u.x.nlin.nii',opath3f,nr),'file')        
                %-- d. apply concatenated xform: volreg/epi2anat/tlrc/NLtlrc; then apply non-linear standard-space warp
                unix(sprintf('3dNwarpApply -master %s/anatQQ.x.nii -dxyz 3 -source %s/func%u_2std.nii -nwarp "%s/anatQQ.x_WARP.nii %s/alg%u.x.affwarp.1D" -prefix ren%u.x.nlin',opath3a,opath1,nr,opath3a,opath3f,nr,nr));
                unix(sprintf('3dAFNItoNIFTI ren%u.x.nlin+tlrc -prefix ren%u.x.nlin.nii',nr,nr));
                unix(sprintf('rm ren%u.x.nlin+tlrc.BRIK ren%u.x.nlin+tlrc.BRIK.gz ren%u.x.nlin+tlrc.HEAD',nr,nr,nr));
                unix(sprintf('mv ren%u.x.nlin.nii %s',nr,opath3f));
            else
                disp('skipping applywarp...')
            end

            if nr==1 && ~exist( sprintf('%s/ren.x.nlin_mask_clean.nii',opath3f),'file')   
                % tidied up functional mask in new space
                unix(sprintf('3dAutomask -prefix %s/ren.x.nlin_mask.nii %s/ren1.x.nlin.nii',opath3f,opath3f));
                unix(sprintf('3dresample -master %s/ren.x.nlin_mask.nii -input %s/anatQQ.x.nii -prefix %s/rm.resam.anat.nii',opath3f,opath3a,opath3f));
                unix(sprintf('3dmask_tool -dilate_input 5 -5 -fill_holes -input %s/rm.resam.anat.nii -prefix %s/rm.resam.anat_mask.nii',opath3f,opath3f))
                unix(sprintf('3dmask_tool -input %s/ren.x.nlin_mask.nii %s/rm.resam.anat_mask.nii -inter -prefix %s/ren.x.nlin_mask_clean.nii',opath3f,opath3f,opath3f))
            else
                disp('skipping newspace masking...')
            end

            tisslist = {'CSF','GM','WM'}; % tissues in increasing order of T1 intensity
            for i=1:3
                if nr==1 && ~exist( sprintf('%s/seg/rm.resam.%s.x.nlin.nii',opath3f,tisslist{i}),'file')
                    unix(sprintf('3dNwarpApply -master %s/anatQQ.x.nii -source %s/seg/anat_seg_%s.nii -nwarp "%s/anatQQ.x_WARP.nii %s/anatQQ.x.aff12.1D" -prefix %s.x.nlin',...
                        opath3a,opath3a,tisslist{i},opath3a,opath3a,tisslist{i}));
                    unix(sprintf('3dAFNItoNIFTI %s.x.nlin+tlrc -prefix %s.x.nlin.nii',tisslist{i},tisslist{i}));
                    unix(sprintf('rm %s.x.nlin+tlrc.BRIK %s.x.nlin+tlrc.BRIK.gz %s.x.nlin+tlrc.HEAD',tisslist{i},tisslist{i},tisslist{i}));
                    unix(sprintf('mv %s.x.nlin.nii %s/seg',tisslist{i},opath3f));
                    unix(sprintf('3dresample -master %s/ren.x.nlin_mask.nii -input  %s/seg/%s.x.nlin.nii -prefix %s/seg/rm.resam.%s.x.nlin.nii',opath3f,opath3f,tisslist{i},opath3f,tisslist{i}));
                else
                    disp('skipping tissue seg warping...')
                end
            end
