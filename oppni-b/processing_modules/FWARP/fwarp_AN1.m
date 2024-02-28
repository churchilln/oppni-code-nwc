function fwarp_AN1( Funcfile_set, prefix_set, odir1, odir2, base_set, Anatloc, ParamCell )

error('need to implement resampling of warped anatomic reference!!')
error('This option currently disabled -- still needs detailed QC!')

%%

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

        %-- a. rigid within-run alignment of all func volumes to func run-1 refbrick
        % disassemble volumes
        unix(sprintf('ImageMath 4 %s/%s_dis.nii.gz TimeSeriesDisassemble %s',pref,prefix_set{nr},Funcfile_set{nr}));
        e = dir(sprintf('%s/%s_dis1*.nii.gz',pref,prefix_set{nr}));
        Vtmp = load_untouch_niiz(sprintf('%s',Funcfile_set{nr}));
        nt = Vtmp.hdr.dime.dim(5);
        ts = Vtmp.hdr.dime.pixdim(5);
        if nt ~=numel(e)
            error('??');
        end
        unix(sprintf('mkdir %s/%s_motmat.aff12.1D',odir1,prefix_set{nr}))
        % go through each one
        mpemat = zeros(nt,6);
        for i=1:nt
            [i i i],
            % get the rigid alignment estimate
            nuff = sprintf('%s/%s_dis1%03u_tu',pref,prefix_set{nr},i-1);
            unix(sprintf(['antsRegistration --dimensionality 3 -f 1 ',...
              '--output [%s_,%s_Warped.nii.gz] ',...
              '--interpolation Linear --use-histogram-matching 1 ',...
              '--winsorize-image-intensities [0.005,0.995] --transform Rigid[ 0.1 ] ',...
              '--metric meansquares[%s.nii.gz,%s/%s_dis1%03u.nii.gz,1] ',...
              '--convergence [15x2,1e-6,4] --shrink-factors 2x1 ',...
              '--smoothing-sigmas 1x0vox'],...
              nuff,nuff,Basefile_pref,pref,prefix_set{nr},i-1));
            % load the affine matrix
            x=load(sprintf('%s_0GenericAffine.mat',nuff));
            % translation
	        dx=x.AffineTransform_double_3_3(10);
	        dy=x.AffineTransform_double_3_3(11);
	        dz=x.AffineTransform_double_3_3(12);
            % rotation
	        rotx = asin(x.AffineTransform_double_3_3(7));
	        roty = atan2(x.AffineTransform_double_3_3(8)/cos(rotx),x.AffineTransform_double_3_3(9)/cos(rotx));
	        rotz = atan2(x.AffineTransform_double_3_3(4)/cos(rotx),x.AffineTransform_double_3_3(1)/cos(rotx));
            % convert to degrees
	        rotx = rotx*360/(2*pi);
	        roty = roty*360/(2*pi);
	        rotz = rotz*360/(2*pi);
            % store to matrix
            mpemat(i,:) = [rotx roty rotz, dx dy dz];
            % push the affine matrices to store:
            unix(sprintf('mv %s_0GenericAffine.mat %s/%s_motmat.aff12.1D/aff1%03u.mat',nuff,odir1,prefix_set{nr},i-1))
        end
        % re-merge motion corrected files >> time_spacing=2.0s / time_origin=0
        ofile = sprintf('%s/%s_motcor.nii.gz',odir1,prefix_set{nr});
        unix(sprintf('ImageMath 4 %s TimeSeriesAssemble %u 0 %s/%s_dis1*.nii.gz',ofile,ts,pref,prefix_set{nr}))
        % also store the mpefiles
        writematrix('mpemat',sprintf('%s/%s_mpe',odir1,prefix_set{nr}));

        if nr==1
    
            %--need to test aligner of motref and t1
            filef = sprintf('%s_masked.nii.gz',Basefile_pref);
            filea = sprintf('%s/anat_procss.nii.gz',Anatloc);
    
            unix(sprintf(['antsRegistration --verbose 1 ',...
            '--dimensionality 3 ',...
            '--interpolation Linear ',...
            '--use-histogram-matching 0 ',...
            '--winsorize-image-intensities [0.005,0.995] ',...
            '--output [%s/alj,%s/alj_Warped.nii.gz,%s/alj_InverseWarped.nii.gz] ',...
            '--initial-moving-transform [%s,%s,1] ',...
            '--transform translation[0.1] ',...
            '--metric MI[%s,%s,1,32,Random,0.25] ',...
            '--convergence [50,1e-6,10] ',...
            '--shrink-factors 1 ',...
            '--smoothing-sigmas 0vox ',...
            '--transform Rigid[0.1] ',...
            '--metric MI[%s,%s,1,32,Random,0.25] ',...
            '--convergence [500x500x250x50,1e-6,10] ',...
            '--shrink-factors 6x4x2x1 ',...
            '--smoothing-sigmas 3x2x1x0vox ',...
            '--transform Affine[0.1] ',...
            '--metric MI[%s,%s,1,32] ',...
            '--convergence [500x500x250x50,1e-6,10] ',...
            '--shrink-factors 6x4x2x1 ',...
            '--smoothing-sigmas 3x2x1x0vox'],...
            odir1,odir1,odir1,filef,filea,filef,filea,filef,filea,filef,filea ));

        end

        for i=1:nt
            [i i i],
            % get the rigid alignment estimate
            nuff = sprintf('%s/%s_wrp1%03u_tu',pref,prefix_set{nr},i-1);

            %-- apply the composite transforms to warp into tempalte space
            unix(sprintf('antsApplyTransforms -d 3 -i %s/%s_dis1%03u.nii.gz -o [%s,0] -t %s/aln1Warp.nii.gz -t %s/aln0GenericAffine.mat -t [%s/alj0GenericAffine.mat,1] -t %s/%s_motmat.aff12.1D/aff1%03u.mat -r %s/anat_warped.nii.gz',...
                pref,prefix_set{nr},i-1, nuff, Anatloc, Anatloc, odir1, odir1,prefix_set{nr},i-1, Anatloc ))
        end
        ofile = sprintf('%s/%s_warped.nii.gz',odir2,prefix_set{nr});
        unix(sprintf('ImageMath 4 %s TimeSeriesAssemble %u 0 %s/%s_wrp1*.nii.gz',ofile,ts,pref,prefix_set{nr}))

        unix(sprintf('rm -rf %s',pref));

    else
        disp('ants-funcwarp already exists!')
    end


end



